"""starfusion — gene-fusion detection from RNA-seq with STAR-Fusion.

Exposes:

- :func:`starfusion_build`, which builds a CTAT genome library from a genome
  FASTA and GTF, a fusion annotation library, Pfam and Dfam HMMs and an
  annotation filter rule (``prep_genome_lib.pl``). It returns the library
  directory, ``<prefix>_genome_lib_build_dir``, as
  ``<prefix>_genome_lib_build_dir.tar.gz``.
- :func:`starfusion_detect`, which detects fusions in one sample's reads, or
  in a STAR ``Chimeric.out.junction``, against a CTAT genome library, writing
  ``<prefix>.fusion_predictions.tsv``, ``<prefix>.abridged.tsv`` and, with
  ``--examine_coding_effect`` in ``args``, ``<prefix>.abridged.coding_effect.tsv``.

A CTAT genome library is a directory of tens of GB for a human genome, and a
directory that big kills Flyte's input staging, so here it travels as one
``.tar.gz`` (unlike upstream, which passes the directory). Every tool that reads
one (STAR-Fusion, FusionInspector, CTAT-Splicing) takes the archive and
extracts it in the task: :data:`EXTRACT_GENOME_LIB`. An archive holding a single
top-level directory (as :func:`starfusion_build` writes, or as the Broad's
prebuilt ``*.plug-n-play.tar.gz`` libraries unpack) uses that directory as the
library.

Resources default to what the upstream tests run with (4 CPUs). A human
library needs far more memory than the defaults: build and detect load the
STAR index (about 32 GB), so override ``resources`` for them.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

STARFUSION_IMAGE = "community.wave.seqera.io/library/dfam_hmmer_minimap2_star-fusion:e285bb3eb373b9a7"
BUILD_RESOURCES = flyte.Resources(cpu=4, memory="16Gi")
DETECT_RESOURCES = flyte.Resources(cpu=4, memory="16Gi")

# Extracts the CTAT genome library archive in the array LIB into a scratch
# directory and sets GENOME_LIB to the library directory: the archive's single
# top-level directory, if it has one, else the extraction root. pigz, if the
# image has it, decompresses faster. Single braces: insert it with `{...}` in
# an rf-string script.
EXTRACT_GENOME_LIB = r"""
        GENOME_LIB_ROOT=$(mktemp -d)
        if command -v pigz > /dev/null; then
            tar -I pigz -xf "${LIB[0]}" -C "$GENOME_LIB_ROOT" || exit 1
        else
            tar -xzf "${LIB[0]}" -C "$GENOME_LIB_ROOT" || exit 1
        fi
        LIB_ENTRIES=("$GENOME_LIB_ROOT"/*)
        if [ ${#LIB_ENTRIES[@]} -eq 1 ] && [ -d "${LIB_ENTRIES[0]}" ]; then
            GENOME_LIB="${LIB_ENTRIES[0]}"
        else
            GENOME_LIB="$GENOME_LIB_ROOT"
        fi
""".strip()

# Archives directory $2 (relative to the working directory) into $1, gzipped
# (with pigz when the image has it).
ARCHIVE_DIR = r"""
        archive_dir() {
            if command -v pigz > /dev/null; then
                tar -I "pigz -p $TASK_CPUS" -cf "$1" "$2" || exit 1
            else
                tar -czf "$1" "$2" || exit 1
            fi
        }
""".strip()

# As upstream (stageInMode 'copy'), the inputs are copied into the working
# directory, the Pfam HMMs are gunzipped and pressed there, and the Dfam
# database is the `.hmm` among the Dfam files (its pressed `.h3*` files beside it).
starfusion_build_cmd = shell.create(
    name="starfusion_build",
    image=STARFUSION_IMAGE,
    resources=BUILD_RESOURCES,
    inputs={
        "fasta": File,
        "gtf": File,
        "fusion_annot_lib": File,
        "pfam": File,
        "dfam": list[File],
        "annot_filter_rule": File,
        "prefix": str,
        "args": str,
    },
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        {ARCHIVE_DIR}
        FA=({{inputs.fasta}}); GTF=({{inputs.gtf}}); ANNOT=({{inputs.fusion_annot_lib}})
        PFAM=({{inputs.pfam}}); DFAM=({{inputs.dfam}}); RULE=({{inputs.annot_filter_rule}})
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        W=$(mktemp -d); cd "$W"
        cp "${{FA[0]}}" "${{GTF[0]}}" "${{ANNOT[0]}}" "${{PFAM[0]}}" "${{DFAM[@]}}" "${{RULE[0]}}" .
        PFAM_GZ=$(basename "${{PFAM[0]}}")
        gunzip "$PFAM_GZ" && hmmpress "${{PFAM_GZ%.gz}}" || exit 1
        DFAM_DB=
        for f in "${{DFAM[@]}}"; do
            case "$f" in *.hmm) DFAM_DB=$(basename "$f") ;; esac
        done
        if [ -z "$DFAM_DB" ]; then echo "no .hmm among the Dfam files" >&2; exit 1; fi
        prep_genome_lib.pl \
            --genome_fa "$(basename "${{FA[0]}}")" \
            --gtf "$(basename "${{GTF[0]}}")" \
            --dfam_db "$DFAM_DB" \
            --pfam_db "${{PFAM_GZ%.gz}}" \
            --fusion_annot_lib "$(basename "${{ANNOT[0]}}")" \
            --annot_filter_rule "$(basename "${{RULE[0]}}")" \
            --CPU $TASK_CPUS \
            --output_dir "${{P}}_genome_lib_build_dir" \
            $ARGS || exit 1
        archive_dir {{outputs.results}}/"${{P}}_genome_lib_build_dir.tar.gz" "${{P}}_genome_lib_build_dir"
    """,
)

# As upstream: reads (single- or paired-end) and/or a STAR chimeric junction
# file; STAR-Fusion's outputs renamed to the prefix.
starfusion_detect_cmd = shell.create(
    name="starfusion_detect",
    image=STARFUSION_IMAGE,
    resources=DETECT_RESOURCES,
    inputs={
        "reads_1": File | None,
        "reads_2": File | None,
        "junction": File | None,
        "genome_lib": File,
        "prefix": str,
        "args": str,
    },
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        shopt -s nullglob
        R1=({{inputs.reads_1}}); R2=({{inputs.reads_2}}); JUNCTION=({{inputs.junction}})
        shopt -u nullglob
        LIB=({{inputs.genome_lib}})
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        {EXTRACT_GENOME_LIB}
        W=$(mktemp -d); cd "$W"
        INPUTS=()
        if [ ${{#R1[@]}} -gt 0 ]; then INPUTS+=(--left_fq "${{R1[0]}}"); fi
        if [ ${{#R2[@]}} -gt 0 ]; then INPUTS+=(--right_fq "${{R2[0]}}"); fi
        if [ ${{#JUNCTION[@]}} -gt 0 ]; then INPUTS+=(-J "${{JUNCTION[0]}}"); fi
        STAR-Fusion \
            --genome_lib_dir "$GENOME_LIB" \
            "${{INPUTS[@]}}" \
            --CPU $TASK_CPUS \
            --output_dir . \
            $ARGS || exit 1
        mv star-fusion.fusion_predictions.tsv "$OUT/$P.fusion_predictions.tsv"
        mv star-fusion.fusion_predictions.abridged.tsv "$OUT/$P.abridged.tsv"
        if [ -f star-fusion.fusion_predictions.abridged.coding_effect.tsv ]; then
            mv star-fusion.fusion_predictions.abridged.coding_effect.tsv "$OUT/$P.abridged.coding_effect.tsv"
        fi
    """,
)

env = flyte.TaskEnvironment.from_task(
    "starfusion",
    starfusion_build_cmd.as_task(),
    starfusion_detect_cmd.as_task(),
)


async def pick(results: Dir, name: str, tool: str) -> File:
    found = await results.get_file(name)
    if found is None:
        raise FileNotFoundError(f"{tool} wrote no {name}")
    return found


async def starfusion_build(
    fasta: File,
    gtf: File,
    fusion_annot_lib: File,
    pfam: File,
    dfam: list[File],
    annot_filter_rule: File,
    prefix: str,
    args: str = "",
) -> File:
    """A CTAT genome library: ``<prefix>_genome_lib_build_dir.tar.gz``.

    ``pfam`` is the gzipped Pfam HMMs (``Pfam-A.hmm.gz``); ``dfam`` is the Dfam HMM file ``*.hmm`` with its
    pressed ``.h3f``, ``.h3i``, ``.h3m`` and ``.h3p`` files; ``annot_filter_rule`` is ``AnnotFilterRule.pm``.
    """
    if not pfam.name.endswith(".gz"):
        raise ValueError(f"starfusion_build needs gzipped Pfam HMMs, got {pfam.name}")
    if not any(f.name.endswith(".hmm") for f in dfam):
        raise ValueError(f"starfusion_build needs a Dfam .hmm among {[f.name for f in dfam]}")
    results = await starfusion_build_cmd(
        fasta=fasta,
        gtf=gtf,
        fusion_annot_lib=fusion_annot_lib,
        pfam=pfam,
        dfam=dfam,
        annot_filter_rule=annot_filter_rule,
        prefix=prefix,
        args=args,
    )
    return await pick(results, f"{prefix}_genome_lib_build_dir.tar.gz", "STAR-Fusion prep_genome_lib.pl")


@dataclass
class StarFusionResult:
    fusions: File  # <prefix>.fusion_predictions.tsv
    abridged: File  # <prefix>.abridged.tsv
    coding_effect: File | None  # <prefix>.abridged.coding_effect.tsv, with --examine_coding_effect


async def starfusion_detect(
    genome_lib: File,
    prefix: str,
    reads_1: File | None = None,
    reads_2: File | None = None,
    junction: File | None = None,
    args: str = "",
) -> StarFusionResult:
    """Fusions in one sample against the CTAT ``genome_lib`` archive (``.tar.gz``).

    Give the reads (``reads_2`` is None for single-end), a STAR ``Chimeric.out.junction`` file, or both.
    Upstream's prefix is ``<sample>.starfusion``.
    """
    if reads_1 is None and junction is None:
        raise ValueError("starfusion_detect needs reads or a chimeric junction file")
    if reads_2 is not None and reads_1 is None:
        raise ValueError("starfusion_detect got reads_2 without reads_1")
    results = await starfusion_detect_cmd(
        reads_1=reads_1, reads_2=reads_2, junction=junction, genome_lib=genome_lib, prefix=prefix, args=args
    )
    return StarFusionResult(
        fusions=await pick(results, f"{prefix}.fusion_predictions.tsv", "STAR-Fusion"),
        abridged=await pick(results, f"{prefix}.abridged.tsv", "STAR-Fusion"),
        coding_effect=await results.get_file(f"{prefix}.abridged.coding_effect.tsv"),
    )
