"""deepvariant — germline small-variant calling with DeepVariant.

Exposes :func:`deepvariant_rundeepvariant`, which runs DeepVariant's whole
pipeline (``run_deepvariant``: make_examples, call_variants,
postprocess_variants) on one sample's BAM or CRAM with its index, and writes,
under upstream's names:

- ``<prefix>.vcf.gz`` with its index ``<prefix>.vcf.gz.tbi`` (or ``.csi``);
- ``<prefix>.g.vcf.gz`` with its index ``<prefix>.g.vcf.gz.tbi`` (or ``.csi``);
- ``<prefix>.visual_report.html``, when DeepVariant writes one.

``fasta`` is the reference, plain or bgzipped; its ``fai`` (and, bgzipped,
its ``gzi``) are linked next to it. ``intervals`` (a BED) restricts calling
(``--regions``) and ``par_bed`` gives the pseudo-autosomal regions of the
haploid contigs (``--par_regions_bed``). ``args`` should set
``--model_type`` (``WGS``, ``WES``, ``PACBIO``, ...), which DeepVariant
requires. The number of shards is the task's CPUs.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

# Upstream's image is docker.io/google/deepvariant:1.10.0; this is Google's own
# copy of it on gcr.io (the same image: its config digest is identical), which
# the cluster can pull.
DEEPVARIANT_IMAGE = "gcr.io/deepvariant-docker/deepvariant:1.10.0"
# Upstream's label is process_high (12 CPUs, 72 GB). This is enough for the test
# data; a whole genome needs far more CPUs (and memory to match) to finish in
# reasonable time: override with deepvariant_rundeepvariant_cmd.as_task().override(resources=...).
DEFAULT_RESOURCES = flyte.Resources(cpu=4, memory="16Gi", disk="20Gi")

deepvariant_rundeepvariant_cmd = shell.create(
    name="deepvariant_rundeepvariant",
    image=DEEPVARIANT_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "reads": File,
        "index": File,
        "intervals": File | None,
        "fasta": File,
        "fai": File,
        "gzi": File | None,
        "par_bed": File | None,
        "prefix": str,
        "args": str,
    },
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        READS=({{inputs.reads}}); INDEX=({{inputs.index}})
        FA=({{inputs.fasta}}); FAI=({{inputs.fai}})
        shopt -s nullglob
        INTERVALS=({{inputs.intervals}}); GZI=({{inputs.gzi}}); PAR=({{inputs.par_bed}})
        shopt -u nullglob
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"

        # Reads and reference beside their indexes, under their own names.
        READS_NAME=$(basename "${{READS[0]}}"); ln -s "${{READS[0]}}" "$READS_NAME"
        INDEX_NAME=$(basename "${{INDEX[0]}}"); ln -s "${{INDEX[0]}}" "$INDEX_NAME"
        FA_NAME=$(basename "${{FA[0]}}"); ln -s "${{FA[0]}}" "$FA_NAME"
        ln -s "${{FAI[0]}}" "$FA_NAME.fai"
        if [ ${{#GZI[@]}} -gt 0 ]; then ln -s "${{GZI[0]}}" "$FA_NAME.gzi"; fi
        EXTRA=()
        if [ ${{#INTERVALS[@]}} -gt 0 ]; then EXTRA+=(--regions="${{INTERVALS[0]}}"); fi
        if [ ${{#PAR[@]}} -gt 0 ]; then EXTRA+=(--par_regions_bed="${{PAR[0]}}"); fi

        export MPLCONFIGDIR=$W/.matplotlib
        mkdir -p "$MPLCONFIGDIR"

        /opt/deepvariant/bin/run_deepvariant \
            --ref="$FA_NAME" \
            --reads="$READS_NAME" \
            --output_vcf="$P.vcf.gz" \
            --output_gvcf="$P.g.vcf.gz" \
            $ARGS \
            "${{EXTRA[@]}}" \
            --intermediate_results_dir=tmp \
            --num_shards=$TASK_CPUS

        shopt -s nullglob
        for f in "$P".vcf.gz* "$P".g.vcf.gz* "$P".visual_report.htm[l]; do mv "$f" "$OUT"/; done
    """,
)

env = flyte.TaskEnvironment.from_task("deepvariant", deepvariant_rundeepvariant_cmd.as_task())


@dataclass
class DeepVariantResult:
    vcf: File  # <prefix>.vcf.gz
    vcf_index: File  # its .tbi or .csi
    gvcf: File  # <prefix>.g.vcf.gz
    gvcf_index: File  # its .tbi or .csi
    report: File | None  # <prefix>.visual_report.html, when written
    results: Dir  # everything above


async def deepvariant_rundeepvariant(
    reads: File,
    index: File,
    fasta: File,
    fai: File,
    intervals: File | None = None,
    gzi: File | None = None,
    par_bed: File | None = None,
    prefix: str = "sample",
    args: str = "",
) -> DeepVariantResult:
    """Call variants in ``reads`` (a BAM or CRAM with its ``index``) against ``fasta``, as VCF and gVCF."""
    results = await deepvariant_rundeepvariant_cmd(
        reads=reads,
        index=index,
        intervals=intervals,
        fasta=fasta,
        fai=fai,
        gzi=gzi,
        par_bed=par_bed,
        prefix=prefix,
        args=args,
    )

    async def required(*names: str) -> File:
        for name in names:
            file = await results.get_file(name)
            if file is not None:
                return file
        raise FileNotFoundError(f"DeepVariant wrote none of {', '.join(names)}")

    return DeepVariantResult(
        vcf=await required(f"{prefix}.vcf.gz"),
        vcf_index=await required(f"{prefix}.vcf.gz.tbi", f"{prefix}.vcf.gz.csi"),
        gvcf=await required(f"{prefix}.g.vcf.gz"),
        gvcf_index=await required(f"{prefix}.g.vcf.gz.tbi", f"{prefix}.g.vcf.gz.csi"),
        report=await results.get_file(f"{prefix}.visual_report.html"),
        results=results,
    )
