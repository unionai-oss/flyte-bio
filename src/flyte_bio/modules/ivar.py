"""ivar — consensus sequences and variant calls for viral amplicon sequencing.

Both tools read a pileup that, as upstream, ``samtools mpileup`` makes from
a BAM against the reference FASTA and pipes in (samtools is in iVar's
image). ``mpileup_args`` carries samtools mpileup's options (upstream's
``ext.args2``, e.g. ``-aa -A -d 0 -Q 0``), ``args`` iVar's (upstream's
``ext.args``). With ``save_mpileup`` the pileup is also kept, as
``<prefix>.mpileup``. As upstream (whose shell doesn't set ``pipefail``),
the pipeline's status is iVar's: a samtools error, such as on an unsorted
BAM, ends the pileup early rather than failing the task, and iVar works on
what it got.

Exposes:

- :func:`ivar_consensus` — the consensus sequence ``<prefix>.fa`` and its
  base qualities ``<prefix>.qual.txt``.
- :func:`ivar_variants` — the variant calls ``<prefix>.tsv``; with a GFF of
  the reference's features, each call's amino-acid change too.

The reference (and its ``.fai`` for ivar_variants) is staged in a working
directory by name, so samtools can index the FASTA beside it if needed.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

IVAR_IMAGE = "quay.io/biocontainers/ivar:1.4.4--h077b44d_0"
# Upstream: process_medium; the tools are single-threaded.
DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory="4Gi")

ivar_consensus_cmd = shell.create(
    name="ivar_consensus",
    image=IVAR_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "bam": File,
        "fasta": File,
        "save_mpileup": bool,
        "prefix": str,
        "mpileup_args": str,
        "args": str,
    },
    defaults={"save_mpileup": False, "mpileup_args": "", "args": ""},
    outputs={"results": Dir},
    script=r"""
        BAM=({inputs.bam}); FA=({inputs.fasta})
        P={inputs.prefix}
        MPILEUP_ARGS={inputs.mpileup_args}
        ARGS={inputs.args}
        OUT={outputs.results}
        W=$(mktemp -d); cd "$W"
        REF=$(basename "${FA[0]}")
        ln -s "${FA[0]}" "$REF"
        SAVE=(cat)
        if [ {inputs.save_mpileup} = true ]; then SAVE=(tee "$P.mpileup"); fi
        set +o pipefail
        samtools \
            mpileup \
            --reference "$REF" \
            $MPILEUP_ARGS \
            "${BAM[0]}" \
            | "${SAVE[@]}" \
            | ivar \
                consensus \
                $ARGS \
                -p "$P" \
            || exit 1
        shopt -s nullglob
        mv "$P.fa" "$P.qual.txt" ./*.mpileup "$OUT"/
    """,
)

ivar_variants_cmd = shell.create(
    name="ivar_variants",
    image=IVAR_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "bam": File,
        "fasta": File,
        "fai": File,
        "gff": File | None,
        "save_mpileup": bool,
        "prefix": str,
        "mpileup_args": str,
        "args": str,
    },
    defaults={"save_mpileup": False, "mpileup_args": "", "args": ""},
    outputs={"results": Dir},
    script=r"""
        BAM=({inputs.bam}); FA=({inputs.fasta}); FAI=({inputs.fai})
        shopt -s nullglob; GFF=({inputs.gff}); shopt -u nullglob
        P={inputs.prefix}
        MPILEUP_ARGS={inputs.mpileup_args}
        ARGS={inputs.args}
        OUT={outputs.results}
        W=$(mktemp -d); cd "$W"
        REF=$(basename "${FA[0]}")
        ln -s "${FA[0]}" "$REF"
        ln -s "${FAI[0]}" "$REF.fai"
        FEATURES=()
        if [ ${#GFF[@]} -gt 0 ]; then FEATURES=(-g "${GFF[0]}"); fi
        SAVE=(cat)
        if [ {inputs.save_mpileup} = true ]; then SAVE=(tee "$P.mpileup"); fi
        set +o pipefail
        samtools \
            mpileup \
            $MPILEUP_ARGS \
            --reference "$REF" \
            "${BAM[0]}" \
            | "${SAVE[@]}" \
            | ivar \
                variants \
                $ARGS \
                "${FEATURES[@]}" \
                -r "$REF" \
                -p "$P" \
            || exit 1
        shopt -s nullglob
        mv "$P.tsv" ./*.mpileup "$OUT"/
    """,
)

env = flyte.TaskEnvironment.from_task("ivar", ivar_consensus_cmd.as_task(), ivar_variants_cmd.as_task())


async def pick_file(results: Dir, name: str, tool: str) -> File:
    found = await results.get_file(name)
    if found is None:
        raise FileNotFoundError(f"{tool} wrote no {name}")
    return found


@dataclass
class IvarConsensusResult:
    fasta: File  # <prefix>.fa
    qual: File  # <prefix>.qual.txt
    mpileup: File | None  # with save_mpileup


@dataclass
class IvarVariantsResult:
    tsv: File  # <prefix>.tsv
    mpileup: File | None  # with save_mpileup


async def ivar_consensus(
    bam: File,
    fasta: File,
    save_mpileup: bool = False,
    prefix: str = "sample",
    mpileup_args: str = "",
    args: str = "",
) -> IvarConsensusResult:
    """The consensus sequence of ``bam``'s reads against the reference ``fasta``."""
    results = await ivar_consensus_cmd(
        bam=bam, fasta=fasta, save_mpileup=save_mpileup, prefix=prefix, mpileup_args=mpileup_args, args=args
    )
    return IvarConsensusResult(
        fasta=await pick_file(results, f"{prefix}.fa", "ivar consensus"),
        qual=await pick_file(results, f"{prefix}.qual.txt", "ivar consensus"),
        mpileup=await pick_file(results, f"{prefix}.mpileup", "samtools mpileup") if save_mpileup else None,
    )


async def ivar_variants(
    bam: File,
    fasta: File,
    fai: File,
    gff: File | None = None,
    save_mpileup: bool = False,
    prefix: str = "sample",
    mpileup_args: str = "",
    args: str = "",
) -> IvarVariantsResult:
    """Variant calls of ``bam``'s reads against the reference ``fasta`` (indexed by ``fai``)."""
    results = await ivar_variants_cmd(
        bam=bam,
        fasta=fasta,
        fai=fai,
        gff=gff,
        save_mpileup=save_mpileup,
        prefix=prefix,
        mpileup_args=mpileup_args,
        args=args,
    )
    return IvarVariantsResult(
        tsv=await pick_file(results, f"{prefix}.tsv", "ivar variants"),
        mpileup=await pick_file(results, f"{prefix}.mpileup", "samtools mpileup") if save_mpileup else None,
    )
