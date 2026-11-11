"""svim — SVIM structural-variant calling from long-read alignments.

Exposes :func:`svim_alignment`, which calls deletions, insertions,
inversions, tandem and interspersed duplications and breakends from one
sample's sorted, indexed long-read BAM against its reference FASTA
(``svim alignment``). As upstream, SVIM writes into a ``<prefix>/`` working
directory and the output is its ``variants.vcf``, as ``<prefix>.vcf``.
Options such as ``--min_mapq`` go in ``args``. The VCF header records the
date.
"""

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

SVIM_IMAGE = "quay.io/biocontainers/svim:2.0.0--pyhdfd78af_0"

# Upstream labels it process_high; SVIM is single-threaded.
DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory="8Gi")

svim_alignment_cmd = shell.create(
    name="svim_alignment",
    image=SVIM_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"bam": File, "bai": File, "fasta": File, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        BAM=({inputs.bam}); BAI=({inputs.bai}); FA=({inputs.fasta})
        P={inputs.prefix}
        ARGS={inputs.args}
        W=$(mktemp -d); cd "$W"
        for F in "${BAM[0]}" "${BAI[0]}" "${FA[0]}"; do ln -s "$F" "$(basename "$F")"; done
        export MPLCONFIGDIR=$(mktemp -d)
        svim alignment \
            $ARGS \
            "$P" \
            "$(basename "${BAM[0]}")" \
            "$(basename "${FA[0]}")" || exit 1
        mv "$P/variants.vcf" {outputs.results}/"$P.vcf"
    """,
)

env = flyte.TaskEnvironment.from_task("svim", svim_alignment_cmd.as_task())


async def svim_alignment(bam: File, bai: File, fasta: File, prefix: str, args: str = "") -> File:
    """Call structural variants in ``bam`` (with its ``bai``) as ``<prefix>.vcf``."""
    results = await svim_alignment_cmd(bam=bam, bai=bai, fasta=fasta, prefix=prefix, args=args)
    vcf = await results.get_file(f"{prefix}.vcf")
    if vcf is None:
        raise FileNotFoundError(f"svim wrote no {prefix}.vcf")
    return vcf
