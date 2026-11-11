"""cutesv — cuteSV structural-variant calling from long-read alignments.

Exposes :func:`cutesv`, which calls insertions, deletions, duplications,
inversions and translocations from one sample's sorted, indexed long-read
BAM against its reference FASTA, optionally only in the regions of a BED
file (``-include_bed``). The output is ``<prefix>.vcf``.

As upstream, cuteSV keeps its signature files in the working directory and
gets the task's CPUs as threads; platform-specific clustering parameters
(``--max_cluster_bias_INS``, ``--diff_ratio_merging_DEL``, …) and
``--genotype`` go in ``args``. The VCF header records the date and command
line.
"""

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

CUTESV_IMAGE = "community.wave.seqera.io/library/cutesv:2.1.4--f0ec240b872b2a98"

# Upstream labels it process_high; these fit test-sized and targeted runs.
DEFAULT_RESOURCES = flyte.Resources(cpu=4, memory="16Gi")

cutesv_cmd = shell.create(
    name="cutesv",
    image=CUTESV_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"bam": File, "bai": File, "fasta": File, "regions": File | None, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        BAM=({{inputs.bam}}); BAI=({{inputs.bai}}); FA=({{inputs.fasta}})
        shopt -s nullglob; REGIONS=({{inputs.regions}}); shopt -u nullglob
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        W=$(mktemp -d); cd "$W"
        for F in "${{BAM[0]}}" "${{BAI[0]}}" "${{FA[0]}}" "${{REGIONS[@]}}"; do ln -s "$F" "$(basename "$F")"; done
        BED=()
        if [ ${{#REGIONS[@]}} -gt 0 ]; then BED=(-include_bed "$(basename "${{REGIONS[0]}}")"); fi
        cuteSV \
            "$(basename "${{BAM[0]}}")" \
            "$(basename "${{FA[0]}}")" \
            "$P.vcf" \
            . \
            --threads $TASK_CPUS \
            "${{BED[@]}}" \
            $ARGS || exit 1
        mv "$P.vcf" {{outputs.results}}/
    """,
)

env = flyte.TaskEnvironment.from_task("cutesv", cutesv_cmd.as_task())


async def cutesv(bam: File, bai: File, fasta: File, prefix: str, regions: File | None = None, args: str = "") -> File:
    """Call structural variants in ``bam`` (with its ``bai``) as ``<prefix>.vcf``."""
    results = await cutesv_cmd(bam=bam, bai=bai, fasta=fasta, regions=regions, prefix=prefix, args=args)
    vcf = await results.get_file(f"{prefix}.vcf")
    if vcf is None:
        raise FileNotFoundError(f"cuteSV wrote no {prefix}.vcf")
    return vcf
