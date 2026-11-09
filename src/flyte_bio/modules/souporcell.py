"""souporcell — genotype-free demultiplexing of pooled single-cell data.

Exposes :func:`souporcell`, which runs ``souporcell_pipeline.py`` on a
droplet-barcoded BAM, its cell barcodes and the reference FASTA: it remaps
the reads, calls variants, clusters the cells into ``clusters`` donors and
estimates ambient RNA. Outputs keep upstream's names (upstream emits them from
its ``<prefix>/`` output directory):

- ``clusters.tsv``: each cell's assignment and cluster log-probabilities;
- ``cluster_genotypes.vcf``: the clusters' genotypes;
- ``ambient_rna.txt``: the estimated ambient RNA fraction.

The pipeline gets the task's CPUs (``-t``); ``args`` go to it too (e.g.
``--skip_remap SKIP_REMAP``). The FASTA is linked into the working
directory, where the pipeline indexes it.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

SOUPORCELL_IMAGE = "community.wave.seqera.io/library/souporcell_gxx:f648658dde2cdd53"
# Upstream: process_medium.
DEFAULT_RESOURCES = flyte.Resources(cpu=4, memory="16Gi")

souporcell_cmd = shell.create(
    name="souporcell",
    image=SOUPORCELL_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"bam": File, "barcodes": File, "fasta": File, "clusters": int, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        BAM=({{inputs.bam}}); BARCODES=({{inputs.barcodes}}); FA=({{inputs.fasta}})
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        BAM_NAME=$(basename "${{BAM[0]}}"); ln -s "${{BAM[0]}}" "$BAM_NAME"
        BARCODES_NAME=$(basename "${{BARCODES[0]}}"); ln -s "${{BARCODES[0]}}" "$BARCODES_NAME"
        FA_NAME=$(basename "${{FA[0]}}"); ln -s "${{FA[0]}}" "$FA_NAME"
        mkdir -p temp
        export TMPDIR=./temp
        souporcell_pipeline.py \
            -i "$BAM_NAME" \
            -b "$BARCODES_NAME" \
            -f "$FA_NAME" \
            -t $TASK_CPUS \
            -o "$P" \
            -k {{inputs.clusters}} \
            $ARGS || exit 1
        mv "$P"/clusters.tsv "$P"/cluster_genotypes.vcf "$P"/ambient_rna.txt "$OUT"/
    """,
)

env = flyte.TaskEnvironment.from_task("souporcell", souporcell_cmd.as_task())


@dataclass
class SouporcellResult:
    clusters: File  # clusters.tsv
    vcf: File  # cluster_genotypes.vcf
    ambient_rna: File  # ambient_rna.txt


async def souporcell(
    bam: File, barcodes: File, fasta: File, clusters: int, prefix: str, args: str = ""
) -> SouporcellResult:
    """souporcell: cluster the cells in ``barcodes`` into ``clusters`` donors by their variants."""
    results = await souporcell_cmd(bam=bam, barcodes=barcodes, fasta=fasta, clusters=clusters, prefix=prefix, args=args)

    async def pick(name: str) -> File:
        found = await results.get_file(name)
        if found is None:
            raise FileNotFoundError(f"souporcell wrote no {name}")
        return found

    return SouporcellResult(
        clusters=await pick("clusters.tsv"),
        vcf=await pick("cluster_genotypes.vcf"),
        ambient_rna=await pick("ambient_rna.txt"),
    )
