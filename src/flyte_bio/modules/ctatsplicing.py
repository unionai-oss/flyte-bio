"""ctatsplicing — cancer-associated aberrant splicing with CTAT-Splicing.

Exposes:

- :func:`ctatsplicing_prepgenomelib`, which adds the cancer intron
  annotations (a ``cancer_introns.*.tsv.gz``) to a CTAT genome library
  (``ctat-splicing-lib-integration.py``), returning the updated library under
  the input archive's name.
- :func:`ctatsplicing_startocancerintrons`, which finds a sample's cancer
  introns from STAR's ``SJ.out.tab`` and ``Chimeric.out.junction`` (and,
  optionally, its coordinate-sorted BAM) against a prepared library
  (``STAR_to_cancer_introns.py``), writing ``<prefix>.cancer.introns``,
  ``<prefix>.cancer.introns.prelim`` and ``<prefix>.introns``, plus, with a
  BAM and ``--vis`` in ``args``, the cancer-intron and gene reads BAMs and an
  IGV report (``<prefix>.ctat-splicing.igv.html``, ``<prefix>.igv.tracks``,
  ``<prefix>.introns.for_IGV.bed``).

The CTAT genome library travels as a ``.tar.gz`` of its directory, as in
:mod:`flyte_bio.modules.starfusion` (upstream passes the directory, copied into
the task for prepgenomelib, which changes it in place). A BAM without an index
is indexed first, as upstream.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS
from flyte_bio.modules.starfusion import ARCHIVE_DIR, EXTRACT_GENOME_LIB

CTATSPLICING_IMAGE = "quay.io/nf-core/ctatsplicing:0.0.3"
# Upstream's process_single, with room to unpack a library.
DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory="6Gi")

ctatsplicing_prepgenomelib_cmd = shell.create(
    name="ctatsplicing_prepgenomelib",
    image=CTATSPLICING_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"genome_lib": File, "cancer_intron_tsv": File},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        {ARCHIVE_DIR}
        LIB=({{inputs.genome_lib}}); TSV=({{inputs.cancer_intron_tsv}})
        OUT={{outputs.results}}
        {EXTRACT_GENOME_LIB}
        /usr/local/src/CTAT-SPLICING/prep_genome_lib/ctat-splicing-lib-integration.py \
            --cancer_introns_tsv "${{TSV[0]}}" \
            --genome_lib_dir "$GENOME_LIB" || exit 1
        cd "$(dirname "$GENOME_LIB")"
        archive_dir "$OUT/$(basename "${{LIB[0]}}")" "$(basename "$GENOME_LIB")"
    """,
)

# As upstream: the BAM (with its index beside it, made if missing) is optional.
ctatsplicing_startocancerintrons_cmd = shell.create(
    name="ctatsplicing_startocancerintrons",
    image=CTATSPLICING_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "split_junction": File,
        "junction": File,
        "bam": File | None,
        "bai": File | None,
        "genome_lib": File,
        "prefix": str,
        "args": str,
    },
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        SJ=({{inputs.split_junction}}); JUNCTION=({{inputs.junction}}); LIB=({{inputs.genome_lib}})
        shopt -s nullglob; BAM=({{inputs.bam}}); BAI=({{inputs.bai}}); shopt -u nullglob
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        {EXTRACT_GENOME_LIB}
        W=$(mktemp -d); cd "$W"
        BAM_ARG=()
        if [ ${{#BAM[@]}} -gt 0 ]; then
            BAM_NAME=$(basename "${{BAM[0]}}")
            ln -s "${{BAM[0]}}" "$BAM_NAME"
            if [ ${{#BAI[@]}} -gt 0 ]; then
                ln -s "${{BAI[0]}}" "$BAM_NAME.bai"
            else
                samtools index "$BAM_NAME" || exit 1
            fi
            BAM_ARG=(--bam_file "$BAM_NAME")
        fi
        /usr/local/src/CTAT-SPLICING/STAR_to_cancer_introns.py \
            --SJ_tab_file "${{SJ[0]}}" \
            --chimJ_file "${{JUNCTION[0]}}" \
            "${{BAM_ARG[@]}}" \
            --output_prefix "$P" \
            --ctat_genome_lib "$GENOME_LIB" \
            $ARGS || exit 1
        shopt -s nullglob
        mv ./*.cancer_intron_reads.sorted.bam ./*.cancer_intron_reads.sorted.bam.bai \
            ./*.gene_reads.sorted.sifted.bam ./*.gene_reads.sorted.sifted.bam.bai \
            ./*.cancer.introns ./*.cancer.introns.prelim ./*"$P".introns ./*.introns.for_IGV.bed \
            ./*.ctat-splicing.igv.html ./*.igv.tracks ./*.chckpts "$OUT"/
    """,
)

env = flyte.TaskEnvironment.from_task(
    "ctatsplicing",
    ctatsplicing_prepgenomelib_cmd.as_task(),
    ctatsplicing_startocancerintrons_cmd.as_task(),
)


async def pick(results: Dir, name: str, tool: str) -> File:
    found = await results.get_file(name)
    if found is None:
        raise FileNotFoundError(f"CTAT-Splicing {tool} wrote no {name}")
    return found


async def ctatsplicing_prepgenomelib(genome_lib: File, cancer_intron_tsv: File) -> File:
    """The CTAT ``genome_lib`` archive (``.tar.gz``) with the cancer introns added, under the same name."""
    results = await ctatsplicing_prepgenomelib_cmd(genome_lib=genome_lib, cancer_intron_tsv=cancer_intron_tsv)
    return await pick(results, genome_lib.name, "ctat-splicing-lib-integration.py")


@dataclass
class CancerIntronsResult:
    cancer_introns: File  # <prefix>.cancer.introns
    cancer_introns_prelim: File  # <prefix>.cancer.introns.prelim
    introns: File  # <prefix>.introns
    results: Dir  # everything written, incl. the optional BAMs, IGV report and checkpoints


async def ctatsplicing_startocancerintrons(
    split_junction: File,
    junction: File,
    genome_lib: File,
    prefix: str,
    bam: File | None = None,
    bai: File | None = None,
    args: str = "",
) -> CancerIntronsResult:
    """A sample's cancer introns from STAR's ``SJ.out.tab`` and ``Chimeric.out.junction``.

    ``genome_lib`` is a CTAT genome library archive prepared by :func:`ctatsplicing_prepgenomelib`. The
    reads BAMs and IGV report (``--vis`` in ``args``) need ``bam``, coordinate-sorted.
    """
    results = await ctatsplicing_startocancerintrons_cmd(
        split_junction=split_junction,
        junction=junction,
        bam=bam,
        bai=bai,
        genome_lib=genome_lib,
        prefix=prefix,
        args=args,
    )
    tool = "STAR_to_cancer_introns.py"
    return CancerIntronsResult(
        cancer_introns=await pick(results, f"{prefix}.cancer.introns", tool),
        cancer_introns_prelim=await pick(results, f"{prefix}.cancer.introns.prelim", tool),
        introns=await pick(results, f"{prefix}.introns", tool),
        results=results,
    )
