"""hiphase — joint phasing of small and structural variants from PacBio HiFi reads with HiPhase.

Exposes :func:`hiphase`, which phases a small-variant VCF (``snv``) and/or a
structural-variant VCF (``sv``), each with its index, using one or more
aligned BAMs (with their indexes) against the reference. The VCFs may hold
several samples; ``samples`` names those to phase (else HiPhase's default).
Outputs keep upstream's names:

- ``<prefix>_snv_phased.vcf.gz`` and ``<prefix>_sv_phased.vcf.gz``, each
  with its index (``.tbi`` or ``.csi``);
- with ``output_bam``, each BAM haplotagged as
  ``<BAM name without .bam>_haplotagged.bam`` with its index;
- with ``summary_file``, ``blocks_file``, ``stats_file`` or
  ``haplotag_file``, ``<prefix>.{summary,blocks,stats,haplotag}.<file_format>``
  (``tsv`` or ``csv``).

``--threads`` follows the task's CPUs. ``args`` go to HiPhase and aren't
shell-quoted. The inputs are linked into the working directory under their
own names, indexes beside their files.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

HIPHASE_IMAGE = "community.wave.seqera.io/library/hiphase:1.6.0--f83f7bc3314ef63b"
# Upstream's label is process_medium; its tests run with 4 CPUs.
DEFAULT_RESOURCES = flyte.Resources(cpu=4, memory="8Gi")

hiphase_cmd = shell.create(
    name="hiphase",
    image=HIPHASE_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "bams": list[File],
        "bais": list[File],
        "snv": File | None,
        "snv_index": File | None,
        "sv": File | None,
        "sv_index": File | None,
        "samples": str,
        "fasta": File,
        "fai": File | None,
        "output_bam": bool,
        "summary_file": bool,
        "blocks_file": bool,
        "stats_file": bool,
        "haplotag_file": bool,
        "file_format": str,
        "prefix": str,
        "args": str,
    },
    defaults={
        "samples": "",
        "output_bam": False,
        "summary_file": False,
        "blocks_file": False,
        "stats_file": False,
        "haplotag_file": False,
        "file_format": "tsv",
        "args": "",
    },
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        BAMS=({{inputs.bams}}); BAIS=({{inputs.bais}}); FA=({{inputs.fasta}})
        shopt -s nullglob
        SNV=({{inputs.snv}}); SNV_INDEX=({{inputs.snv_index}}); SV=({{inputs.sv}}); SV_INDEX=({{inputs.sv_index}})
        FAI=({{inputs.fai}})
        shopt -u nullglob
        SAMPLE_NAMES={{inputs.samples}}
        P={{inputs.prefix}}
        FORMAT={{inputs.file_format}}
        ARGS={{inputs.args}}
        W=$(mktemp -d); cd "$W"
        # Each file under its own name, each index beside its file (the i-th BAI goes with the i-th BAM).
        link_with_index() {{ ln -s "$1" "$(basename "$1")"; ln -s "$2" "$(basename "$1").${{2##*.}}"; }}
        for I in "${{!BAMS[@]}}"; do link_with_index "${{BAMS[$I]}}" "${{BAIS[$I]}}"; done
        if [ ${{#SNV[@]}} -gt 0 ]; then link_with_index "${{SNV[0]}}" "${{SNV_INDEX[0]}}"; fi
        if [ ${{#SV[@]}} -gt 0 ]; then link_with_index "${{SV[0]}}" "${{SV_INDEX[0]}}"; fi
        FA_NAME=$(basename "${{FA[0]}}"); ln -s "${{FA[0]}}" "$FA_NAME"
        if [ ${{#FAI[@]}} -gt 0 ]; then ln -s "${{FAI[0]}}" "$FA_NAME.fai"; fi

        OPTIONS=()
        for S in $SAMPLE_NAMES; do OPTIONS+=(--sample-name "$S"); done
        for B in "${{BAMS[@]}}"; do
            NAME=$(basename "$B")
            OPTIONS+=(--bam "$NAME")
            if [ {{inputs.output_bam}} = true ]; then OPTIONS+=(--output-bam "${{NAME%.*}}_haplotagged.bam"); fi
        done
        if [ ${{#SNV[@]}} -gt 0 ]; then
            OPTIONS+=(--vcf "$(basename "${{SNV[0]}}")" --output-vcf "${{P}}_snv_phased.vcf.gz")
        fi
        if [ ${{#SV[@]}} -gt 0 ]; then
            OPTIONS+=(--vcf "$(basename "${{SV[0]}}")" --output-vcf "${{P}}_sv_phased.vcf.gz")
        fi
        if [ {{inputs.summary_file}} = true ]; then OPTIONS+=(--summary-file "$P.summary.$FORMAT"); fi
        if [ {{inputs.blocks_file}} = true ]; then OPTIONS+=(--blocks-file "$P.blocks.$FORMAT"); fi
        if [ {{inputs.stats_file}} = true ]; then OPTIONS+=(--stats-file "$P.stats.$FORMAT"); fi
        if [ {{inputs.haplotag_file}} = true ]; then OPTIONS+=(--haplotag-file "$P.haplotag.$FORMAT"); fi

        hiphase \
            $ARGS \
            --threads $TASK_CPUS \
            --reference "$FA_NAME" \
            "${{OPTIONS[@]}}" || exit 1

        shopt -s nullglob
        OUTPUTS=(./*_phased.vcf.gz ./*_phased.vcf.gz.tbi ./*_phased.vcf.gz.csi ./*_haplotagged.bam*)
        OUTPUTS+=("$P".summary.* "$P".blocks.* "$P".stats.* "$P".haplotag.*)
        mv "${{OUTPUTS[@]}}" {{outputs.results}}/
    """,
)

env = flyte.TaskEnvironment.from_task("hiphase", hiphase_cmd.as_task())


@dataclass
class HiphaseResult:
    snv_vcf: File | None  # <prefix>_snv_phased.vcf.gz
    snv_vcf_index: File | None  # its .tbi or .csi
    sv_vcf: File | None  # <prefix>_sv_phased.vcf.gz
    sv_vcf_index: File | None  # its .tbi or .csi
    summary: File | None  # <prefix>.summary.<file_format>
    blocks: File | None  # <prefix>.blocks.<file_format>
    stats: File | None  # <prefix>.stats.<file_format>
    haplotag: File | None  # <prefix>.haplotag.<file_format>
    bams: list[File]  # <BAM name without .bam>_haplotagged.bam, with output_bam
    bam_indexes: list[File]  # their .bai or .csi
    results: Dir  # everything above


async def hiphase(
    bams: list[File],
    bais: list[File],
    fasta: File,
    snv: File | None = None,
    snv_index: File | None = None,
    sv: File | None = None,
    sv_index: File | None = None,
    samples: list[str] | None = None,
    fai: File | None = None,
    output_bam: bool = False,
    summary_file: bool = False,
    blocks_file: bool = False,
    stats_file: bool = False,
    haplotag_file: bool = False,
    file_format: str = "tsv",
    prefix: str = "sample",
    args: str = "",
) -> HiphaseResult:
    """Phase ``snv`` and/or ``sv`` (VCFs with their indexes) with the reads in ``bams``."""
    results = await hiphase_cmd(
        bams=bams,
        bais=bais,
        snv=snv,
        snv_index=snv_index,
        sv=sv,
        sv_index=sv_index,
        samples=" ".join(samples or []),
        fasta=fasta,
        fai=fai,
        output_bam=output_bam,
        summary_file=summary_file,
        blocks_file=blocks_file,
        stats_file=stats_file,
        haplotag_file=haplotag_file,
        file_format=file_format,
        prefix=prefix,
        args=args,
    )
    files = sorted([f async for f in results.walk()], key=lambda f: f.name)

    def first(*names: str) -> File | None:
        return next((f for name in names for f in files if f.name == name), None)

    return HiphaseResult(
        snv_vcf=first(f"{prefix}_snv_phased.vcf.gz"),
        snv_vcf_index=first(f"{prefix}_snv_phased.vcf.gz.tbi", f"{prefix}_snv_phased.vcf.gz.csi"),
        sv_vcf=first(f"{prefix}_sv_phased.vcf.gz"),
        sv_vcf_index=first(f"{prefix}_sv_phased.vcf.gz.tbi", f"{prefix}_sv_phased.vcf.gz.csi"),
        summary=first(f"{prefix}.summary.{file_format}"),
        blocks=first(f"{prefix}.blocks.{file_format}"),
        stats=first(f"{prefix}.stats.{file_format}"),
        haplotag=first(f"{prefix}.haplotag.{file_format}"),
        bams=[f for f in files if f.name.endswith("_haplotagged.bam")],
        bam_indexes=[f for f in files if f.name.endswith(("_haplotagged.bam.bai", "_haplotagged.bam.csi"))],
        results=results,
    )
