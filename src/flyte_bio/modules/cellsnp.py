"""cellsnp — per-cell allele counts at known SNPs with cellsnp-lite (mode 1a).

Exposes :func:`cellsnp_modea`, which pileups a droplet-barcoded BAM at the
SNPs of ``region_vcf`` (or genome-wide without it) for the listed cell
barcodes. Outputs keep upstream's names:

- ``<prefix>.base.vcf.gz``: the SNPs found, aggregated over cells;
- ``<prefix>.cells.vcf.gz``: per-cell genotypes (only with ``--genotype`` in
  ``args``);
- ``<prefix>.samples.tsv``: the cell barcodes, in matrix column order;
- ``<prefix>.tag.{AD,DP,OTH}.mtx``: per-cell alternative-allele, total and
  other-allele depths (Matrix Market).

The BAM is linked beside its index. cellsnp-lite gets the task's CPUs
(``--nproc``); ``args`` go to it too (e.g. ``--genotype --minMAF 0.1``).
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

CELLSNP_IMAGE = "community.wave.seqera.io/library/cellsnp-lite:1.2.3--777dfc22dc268b85"
# Upstream: process_medium.
DEFAULT_RESOURCES = flyte.Resources(cpu=4, memory="8Gi")

cellsnp_modea_cmd = shell.create(
    name="cellsnp_modea",
    image=CELLSNP_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"bam": File, "bai": File, "region_vcf": File | None, "barcodes": File, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        BAM=({{inputs.bam}}); BAI=({{inputs.bai}}); BARCODES=({{inputs.barcodes}})
        shopt -s nullglob; REGIONS=({{inputs.region_vcf}}); shopt -u nullglob
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        BAM_NAME=$(basename "${{BAM[0]}}")
        ln -s "${{BAM[0]}}" "$BAM_NAME"; ln -s "${{BAI[0]}}" "$BAM_NAME.bai"
        REGION_OPT=()
        if [ ${{#REGIONS[@]}} -gt 0 ]; then REGION_OPT=(-R "${{REGIONS[0]}}"); fi
        cellsnp-lite -s "$BAM_NAME" \
            -b "${{BARCODES[0]}}" \
            "${{REGION_OPT[@]}}" \
            -O . \
            --gzip \
            --nproc $TASK_CPUS \
            $ARGS || exit 1
        mv cellSNP.base.vcf.gz "$OUT/$P.base.vcf.gz"
        if [[ "$ARGS" == *"--genotype"* ]]; then
            mv cellSNP.cells.vcf.gz "$OUT/$P.cells.vcf.gz"
        fi
        mv cellSNP.tag.AD.mtx "$OUT/$P.tag.AD.mtx"
        mv cellSNP.tag.DP.mtx "$OUT/$P.tag.DP.mtx"
        mv cellSNP.tag.OTH.mtx "$OUT/$P.tag.OTH.mtx"
        mv cellSNP.samples.tsv "$OUT/$P.samples.tsv"
    """,
)

env = flyte.TaskEnvironment.from_task("cellsnp", cellsnp_modea_cmd.as_task())


@dataclass
class CellsnpResult:
    base: File  # <prefix>.base.vcf.gz
    cells: File | None  # <prefix>.cells.vcf.gz, with --genotype
    samples: File  # <prefix>.samples.tsv
    allele_depth: File  # <prefix>.tag.AD.mtx
    depth_coverage: File  # <prefix>.tag.DP.mtx
    depth_other: File  # <prefix>.tag.OTH.mtx


async def cellsnp_modea(
    bam: File,
    bai: File,
    barcodes: File,
    prefix: str,
    region_vcf: File | None = None,
    args: str = "",
) -> CellsnpResult:
    """cellsnp-lite mode 1a: per-cell allele depths at ``region_vcf``'s SNPs for the cells in ``barcodes``."""
    results = await cellsnp_modea_cmd(
        bam=bam, bai=bai, region_vcf=region_vcf, barcodes=barcodes, prefix=prefix, args=args
    )

    async def pick(name: str) -> File:
        found = await results.get_file(name)
        if found is None:
            raise FileNotFoundError(f"cellsnp_modea wrote no {name}")
        return found

    return CellsnpResult(
        base=await pick(f"{prefix}.base.vcf.gz"),
        cells=await results.get_file(f"{prefix}.cells.vcf.gz"),
        samples=await pick(f"{prefix}.samples.tsv"),
        allele_depth=await pick(f"{prefix}.tag.AD.mtx"),
        depth_coverage=await pick(f"{prefix}.tag.DP.mtx"),
        depth_other=await pick(f"{prefix}.tag.OTH.mtx"),
    )
