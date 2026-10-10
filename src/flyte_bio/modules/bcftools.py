"""bcftools — VCF/BCF utilities.

Exposes :func:`bcftools_stats`, which writes ``bcftools stats`` for a VCF to
``<prefix>.bcftools_stats.txt``. The VCF (and a reference FASTA, for indel
context) are linked into the working directory under their own names, so the
command line bcftools records matches upstream's.
"""

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

BCFTOOLS_IMAGE = "community.wave.seqera.io/library/bcftools_htslib:1.23.1--9f08ec665533d64a"

DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory="2Gi")

bcftools_stats_cmd = shell.create(
    name="bcftools_stats",
    image=BCFTOOLS_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"vcf": File, "tbi": File | None, "fasta": File | None, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        V=({inputs.vcf})
        shopt -s nullglob; T=({inputs.tbi}); FA=({inputs.fasta}); shopt -u nullglob
        ARGS={inputs.args}
        W=$(mktemp -d); cd "$W"
        NAME=$(basename "${V[0]}"); ln -s "${V[0]}" "$NAME"
        if [ ${#T[@]} -gt 0 ]; then ln -s "${T[0]}" "$NAME.tbi"; fi
        REF=()
        if [ ${#FA[@]} -gt 0 ]; then
            FA_NAME=$(basename "${FA[0]}"); ln -s "${FA[0]}" "$FA_NAME"; REF=(--fasta-ref "$FA_NAME")
        fi
        bcftools stats $ARGS "${REF[@]}" "$NAME" > {outputs.results}/{inputs.prefix}.bcftools_stats.txt
    """,
)


env = flyte.TaskEnvironment.from_task(
    "bcftools",
    bcftools_stats_cmd.as_task(),
)


async def bcftools_stats(
    vcf: File, prefix: str, tbi: File | None = None, fasta: File | None = None, args: str = ""
) -> File:
    """``<prefix>.bcftools_stats.txt``: bcftools' statistics of ``vcf``."""
    results = await bcftools_stats_cmd(vcf=vcf, tbi=tbi, fasta=fasta, prefix=prefix, args=args)
    stats = await results.get_file(f"{prefix}.bcftools_stats.txt")
    if stats is None:
        raise FileNotFoundError(f"bcftools stats wrote no {prefix}.bcftools_stats.txt")
    return stats
