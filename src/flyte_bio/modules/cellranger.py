"""cellranger — 10x Genomics Cell Ranger: V(D)J references and immune profiling.

Exposes:

- :func:`cellranger_mkvdjref`, which builds a Cell Ranger V(D)J reference
  (``cellranger mkvdjref``) from a genome FASTA and a GTF, or from
  ready-made V(D)J segment sequences (``seqs``). It returns the
  ``<reference_name>`` folder (``reference.json``, ``fasta/regions.fa``, ...)
  as ``<reference_name>.tar.gz``, plus its ``regions.fa`` and
  ``reference.json``.
- :func:`cellranger_vdj`, which assembles one sample's V(D)J transcripts and
  calls clonotypes (``cellranger vdj``) against such a reference, extracting
  it first. ``reads`` are the sample's FASTQs, named as ``bcl2fastq`` names
  them (``<sample>_S1_L001_R1_001.fastq.gz``). It returns Cell Ranger's
  ``<prefix>/outs`` folder (``airr_rearrangement.tsv``,
  ``filtered_contig_annotations.csv``, ``clonotypes.csv``,
  ``metrics_summary.csv``, ``web_summary.html``, ``vdj_reference/``, ...),
  and the same as ``<prefix>.outs.tar.gz`` (an archive of ``<prefix>/outs``)
  to pass on whole.

Both pass the task's cores and memory as ``--localcores`` and ``--localmem``.
Telemetry is disabled (``TENX_DISABLE_TELEMETRY``) and ``HOME`` is the task's
working directory.

The container is upstream's build of 10x's Cell Ranger, distributed under
10x Genomics' end-user license: it may be used only with data from 10x
Genomics products, and not redistributed.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

CELLRANGER_IMAGE = "quay.io/nf-core/cellranger:10.1.0"
NO_TELEMETRY = {"TENX_DISABLE_TELEMETRY": "1"}
# Upstream's process_high asks for 12 CPUs and 72 GB; its tests run in 6 GB.
# Real samples want 16-64 GiB (10x recommends 64 GB or more for vdj): raise
# with <cmd>.as_task().override(resources=...).
MKVDJREF_RESOURCES = flyte.Resources(cpu=4, memory="8Gi", disk="20Gi")
VDJ_RESOURCES = flyte.Resources(cpu=4, memory="16Gi", disk="50Gi")

cellranger_mkvdjref_cmd = shell.create(
    name="cellranger_mkvdjref",
    image=CELLRANGER_IMAGE,
    resources=MKVDJREF_RESOURCES,
    env_vars=NO_TELEMETRY,
    inputs={
        "fasta": File | None,
        "gtf": File | None,
        "seqs": File | None,
        "reference_name": str,
        "args": str,
    },
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        shopt -s nullglob
        FA=({{inputs.fasta}}); GTF=({{inputs.gtf}}); SEQS=({{inputs.seqs}})
        shopt -u nullglob
        NAME={{inputs.reference_name}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        export HOME="$W"

        EXTRA=()
        if [ ${{#GTF[@]}} -gt 0 ]; then EXTRA+=(--genes "${{GTF[0]}}"); fi
        if [ ${{#FA[@]}} -gt 0 ]; then EXTRA+=(--fasta "${{FA[0]}}"); fi
        if [ ${{#SEQS[@]}} -gt 0 ]; then EXTRA+=(--seqs "${{SEQS[0]}}"); fi

        cellranger \
            mkvdjref \
            --genome="$NAME" \
            "${{EXTRA[@]}}" \
            --localcores=$TASK_CPUS \
            --localmem=$(( TASK_MEMORY_MB / 1024 )) \
            $ARGS || exit 1

        tar -czf "$OUT/$NAME.tar.gz" "$NAME"
        cp "$NAME/fasta/regions.fa" "$NAME/reference.json" "$OUT"/
    """,
)

cellranger_vdj_cmd = shell.create(
    name="cellranger_vdj",
    image=CELLRANGER_IMAGE,
    resources=VDJ_RESOURCES,
    env_vars=NO_TELEMETRY,
    inputs={"reads": list[File], "reference": File, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        READS=({{inputs.reads}})
        REF_ARCHIVE=({{inputs.reference}})
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        export HOME="$W"

        mkdir reference fastqs
        tar -xzf "${{REF_ARCHIVE[0]}}" -C reference
        if [ -f reference/reference.json ]; then REF="$W/reference"; else REF=(reference/*/); REF="$W/${{REF[0]%/}}"; fi
        for R in "${{READS[@]}}"; do ln -s "$R" fastqs/; done

        cellranger \
            vdj \
            --id="$P" \
            --fastqs="$W/fastqs" \
            --reference="$REF" \
            --localcores=$TASK_CPUS \
            --localmem=$(( TASK_MEMORY_MB / 1024 )) \
            $ARGS || exit 1

        tar -czhf "$OUT/$P.outs.tar.gz" "$P/outs"
        mkdir "$OUT/outs"
        cp -rL "$P"/outs/. "$OUT/outs"/
    """,
)

env = flyte.TaskEnvironment.from_task("cellranger", cellranger_mkvdjref_cmd.as_task(), cellranger_vdj_cmd.as_task())


@dataclass
class VdjReference:
    reference: File  # <reference_name>.tar.gz, an archive of the <reference_name> folder
    regions: File  # <reference_name>/fasta/regions.fa
    reference_json: File  # <reference_name>/reference.json


async def cellranger_mkvdjref(
    reference_name: str,
    fasta: File | None = None,
    gtf: File | None = None,
    seqs: File | None = None,
    args: str = "",
) -> VdjReference:
    """A Cell Ranger V(D)J reference from ``fasta`` and ``gtf``, or from ``seqs``: ``<reference_name>.tar.gz``."""
    if seqs is None and (fasta is None or gtf is None):
        raise ValueError("cellranger_mkvdjref needs a fasta and gtf, or seqs")
    results = await cellranger_mkvdjref_cmd(fasta=fasta, gtf=gtf, seqs=seqs, reference_name=reference_name, args=args)
    files = {}
    for name in (f"{reference_name}.tar.gz", "regions.fa", "reference.json"):
        files[name] = await results.get_file(name)
        if files[name] is None:
            raise FileNotFoundError(f"cellranger mkvdjref wrote no {name}")
    return VdjReference(
        reference=files[f"{reference_name}.tar.gz"],
        regions=files["regions.fa"],
        reference_json=files["reference.json"],
    )


@dataclass
class VdjResult:
    outs: Dir  # <prefix>/outs
    archive: File  # <prefix>.outs.tar.gz, an archive of <prefix>/outs
    airr_rearrangement: File | None  # outs/airr_rearrangement.tsv
    filtered_contig_annotations: File | None  # outs/filtered_contig_annotations.csv
    metrics_summary: File | None  # outs/metrics_summary.csv
    web_summary: File | None  # outs/web_summary.html


async def cellranger_vdj(reads: list[File], reference: File, prefix: str, args: str = "") -> VdjResult:
    """Cell Ranger V(D)J assembly and clonotypes of one sample's ``reads`` against ``reference`` (a ``.tar.gz``)."""
    if not reads:
        raise ValueError("cellranger_vdj needs the sample's FASTQs")
    results = await cellranger_vdj_cmd(reads=reads, reference=reference, prefix=prefix, args=args)
    archive = await results.get_file(f"{prefix}.outs.tar.gz")
    if archive is None:
        raise FileNotFoundError(f"cellranger vdj wrote no {prefix}.outs.tar.gz")
    outs = Dir.from_existing_remote(results.path.rstrip("/") + "/outs")
    return VdjResult(
        outs=outs,
        archive=archive,
        airr_rearrangement=await outs.get_file("airr_rearrangement.tsv"),
        filtered_contig_annotations=await outs.get_file("filtered_contig_annotations.csv"),
        metrics_summary=await outs.get_file("metrics_summary.csv"),
        web_summary=await outs.get_file("web_summary.html"),
    )
