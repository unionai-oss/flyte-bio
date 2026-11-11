"""longphase — LongPhase phasing of variants and haplotagging of long reads.

Exposes:

- :func:`longphase_phase` — phase one sample's SNVs (``--snp-file``) with its
  long-read BAMs, optionally co-phasing structural variants (``--sv-file``)
  and modifications (``--mod-file``): ``<prefix>.vcf.gz``, and
  ``<prefix>_SV.vcf.gz`` / ``<prefix>_mod.vcf.gz`` when given.
- :func:`longphase_haplotag` — tag each read of a BAM with its haplotype
  (``HP``) and phase set (``PS``) from phased variants: ``<prefix>.bam``, or
  ``<prefix>.cram`` with ``--cram`` in ``args``, and with ``--log`` the
  per-read log ``<prefix>.log``.

As upstream:

- Inputs keep their names in the working directory: LongPhase records its
  command line (with ``--threads``) in the VCF header and the BAM's ``@PG``.
- ``--threads`` follows the task's CPUs; options such as ``--ont`` or
  ``--pb`` go in ``args``. For :func:`longphase_phase`, ``args2`` go to
  ``bgzip``.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

LONGPHASE_IMAGE = "community.wave.seqera.io/library/htslib_longphase:3418df83ebc1fb9a"

# Upstream labels both process_medium; the 4 CPUs match upstream's test CPUs,
# which its snapshots record.
DEFAULT_RESOURCES = flyte.Resources(cpu=4, memory="16Gi")

# Links each input into the working directory under its own name.
LINK_INPUTS = r"""
        W=$(mktemp -d); cd "$W"
        for F in "${BAMS[@]}" "${BAIS[@]}" "${SNVS[0]}" "${SVS[@]}" "${MODS[@]}" "${FA[0]}" "${FAI[0]}"; do
            ln -s "$F" "$(basename "$F")"
        done
        EXTRA=()
        if [ ${#SVS[@]} -gt 0 ]; then EXTRA+=(--sv-file "$(basename "${SVS[0]}")"); fi
        if [ ${#MODS[@]} -gt 0 ]; then EXTRA+=(--mod-file "$(basename "${MODS[0]}")"); fi
"""

INPUTS = {
    "snvs": File,
    "svs": File | None,
    "mods": File | None,
    "fasta": File,
    "fai": File,
    "prefix": str,
    "args": str,
}

longphase_phase_cmd = shell.create(
    name="longphase_phase",
    image=LONGPHASE_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"bams": list[File], "bais": list[File], **INPUTS, "args2": str},
    defaults={"args": "", "args2": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        BAMS=({{inputs.bams}}); SNVS=({{inputs.snvs}}); FA=({{inputs.fasta}}); FAI=({{inputs.fai}})
        shopt -s nullglob; BAIS=({{inputs.bais}}); SVS=({{inputs.svs}}); MODS=({{inputs.mods}}); shopt -u nullglob
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        ARGS2={{inputs.args2}}
        {LINK_INPUTS}
        READS=()
        for F in "${{BAMS[@]}}"; do READS+=(-b "$(basename "$F")"); done
        longphase \
            phase \
            $ARGS \
            --threads $TASK_CPUS \
            -o "$P" \
            --reference "$(basename "${{FA[0]}}")" \
            --snp-file "$(basename "${{SNVS[0]}}")" \
            "${{READS[@]}}" \
            "${{EXTRA[@]}}" || exit 1
        bgzip --threads $TASK_CPUS $ARGS2 "$P"*.vcf || exit 1
        for F in "$P.vcf.gz" "${{P}}_SV.vcf.gz" "${{P}}_mod.vcf.gz"; do
            if [ -f "$F" ]; then mv "$F" {{outputs.results}}/; fi
        done
    """,
)

longphase_haplotag_cmd = shell.create(
    name="longphase_haplotag",
    image=LONGPHASE_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"bam": File, "bai": File | None, **INPUTS},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        BAMS=({{inputs.bam}}); SNVS=({{inputs.snvs}}); FA=({{inputs.fasta}}); FAI=({{inputs.fai}})
        shopt -s nullglob; BAIS=({{inputs.bai}}); SVS=({{inputs.svs}}); MODS=({{inputs.mods}}); shopt -u nullglob
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        {LINK_INPUTS}
        longphase \
            haplotag \
            $ARGS \
            --threads $TASK_CPUS \
            -o "$P" \
            --reference "$(basename "${{FA[0]}}")" \
            --snp-file "$(basename "${{SNVS[0]}}")" \
            --bam "$(basename "${{BAMS[0]}}")" \
            "${{EXTRA[@]}}" || exit 1
        if [ -f "$P.out" ]; then mv "$P.out" "$P.log"; fi
        for F in "$P.bam" "$P.cram" "$P.log"; do
            if [ -f "$F" ] && [ ! -L "$F" ]; then mv "$F" {{outputs.results}}/; fi
        done
    """,
)

env = flyte.TaskEnvironment.from_task("longphase", longphase_phase_cmd.as_task(), longphase_haplotag_cmd.as_task())


@dataclass
class LongphasePhaseResult:
    snv_vcf: File  # <prefix>.vcf.gz
    sv_vcf: File | None  # <prefix>_SV.vcf.gz, with svs
    mod_vcf: File | None  # <prefix>_mod.vcf.gz, with mods
    results: Dir


@dataclass
class LongphaseHaplotagResult:
    alignment: File  # <prefix>.bam, or <prefix>.cram with --cram
    log: File | None  # <prefix>.log, with --log
    results: Dir


async def longphase_phase(
    bams: list[File],
    bais: list[File],
    snvs: File,
    fasta: File,
    fai: File,
    prefix: str,
    svs: File | None = None,
    mods: File | None = None,
    args: str = "",
    args2: str = "",
) -> LongphasePhaseResult:
    """Phase ``snvs`` (and ``svs``, ``mods``) with the reads in ``bams``."""
    results = await longphase_phase_cmd(
        bams=bams,
        bais=bais,
        snvs=snvs,
        svs=svs,
        mods=mods,
        fasta=fasta,
        fai=fai,
        prefix=prefix,
        args=args,
        args2=args2,
    )
    snv_vcf = await results.get_file(f"{prefix}.vcf.gz")
    if snv_vcf is None:
        raise FileNotFoundError(f"longphase phase wrote no {prefix}.vcf.gz")
    return LongphasePhaseResult(
        snv_vcf=snv_vcf,
        sv_vcf=await results.get_file(f"{prefix}_SV.vcf.gz"),
        mod_vcf=await results.get_file(f"{prefix}_mod.vcf.gz"),
        results=results,
    )


async def longphase_haplotag(
    bam: File,
    snvs: File,
    fasta: File,
    fai: File,
    prefix: str,
    bai: File | None = None,
    svs: File | None = None,
    mods: File | None = None,
    args: str = "",
) -> LongphaseHaplotagResult:
    """Tag the reads in ``bam`` with their haplotypes from phased ``snvs`` (and ``svs``, ``mods``)."""
    name = bam.name or bam.path.rsplit("/", 1)[-1]
    if name in (f"{prefix}.bam", f"{prefix}.cram"):
        raise ValueError("longphase_haplotag: input and output names are the same; pass another prefix")
    results = await longphase_haplotag_cmd(
        bam=bam, bai=bai, snvs=snvs, svs=svs, mods=mods, fasta=fasta, fai=fai, prefix=prefix, args=args
    )
    alignment = await results.get_file(f"{prefix}.cram") or await results.get_file(f"{prefix}.bam")
    if alignment is None:
        raise FileNotFoundError(f"longphase haplotag wrote no {prefix}.bam or {prefix}.cram")
    return LongphaseHaplotagResult(alignment=alignment, log=await results.get_file(f"{prefix}.log"), results=results)
