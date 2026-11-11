"""ucsc — UCSC genome browser utilities.

Exposes :func:`bedclip`, which drops bedGraph records outside the
chromosome bounds in a chrom.sizes file, and :func:`bedgraphtobigwig`, which
converts a (clipped, sorted) bedGraph to bigWig, and :func:`wigtobigwig`,
which converts a wiggle file (plain or gzipped) to bigWig. Outputs keep
upstream's names: ``<prefix>.bedGraph``, ``<prefix>.bigWig`` and
``<prefix>.bw``.

:func:`bigwigaverageoverbed` averages a bigWig over each region of a BED,
writing the table ``<prefix>.tab`` (name, size, covered bases, sum, mean over
the region, mean over covered bases).

It also exposes :func:`gtftogenepred`, which converts a GTF to genePred,
``<prefix>.genepred``. As upstream, when ``args`` has both ``-genePredExt``
and ``-geneNameAsName2`` it also writes ``<prefix>.refflat``, the refFlat
table (gene name first) Picard's CollectRnaSeqMetrics reads.

:func:`liftover` converts a BED's coordinates to another assembly with a
chain file (``liftOver``), writing ``<prefix>.lifted.bed`` and the records it
couldn't convert, ``<prefix>.unlifted.bed``.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

# Pinned biocontainer URIs. Upstream's main.nf names the bare
# `biocontainers/...` (Docker Hub) — use the quay.io mirrors.
BEDCLIP_IMAGE = "quay.io/biocontainers/ucsc-bedclip:377--h0b8a92a_2"
BEDGRAPHTOBIGWIG_IMAGE = "quay.io/biocontainers/ucsc-bedgraphtobigwig:469--h9b8f530_0"
WIGTOBIGWIG_IMAGE = "community.wave.seqera.io/library/ucsc-wigtobigwig:482--4eb74e15915ee878"
GTFTOGENEPRED_IMAGE = "community.wave.seqera.io/library/ucsc-gtftogenepred:482--ace98be9f7cb069f"
BIGWIGAVERAGEOVERBED_IMAGE = "community.wave.seqera.io/library/ucsc-bigwigaverageoverbed:482--09b4685dee7c8b33"
LIFTOVER_IMAGE = "community.wave.seqera.io/library/ucsc-liftover:482--f8639164c0f3b5bf"

DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory="4Gi")

bedclip_cmd = shell.create(
    name="ucsc_bedclip",
    image=BEDCLIP_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"bedgraph": File, "sizes": File, "prefix": str},
    outputs={"results": Dir},
    script=r"""
        bedClip {inputs.bedgraph} {inputs.sizes} {outputs.results}/{inputs.prefix}.bedGraph
    """,
)

bedgraphtobigwig_cmd = shell.create(
    name="ucsc_bedgraphtobigwig",
    image=BEDGRAPHTOBIGWIG_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"bedgraph": File, "sizes": File, "prefix": str},
    outputs={"results": Dir},
    script=r"""
        bedGraphToBigWig {inputs.bedgraph} {inputs.sizes} {outputs.results}/{inputs.prefix}.bigWig
    """,
)

wigtobigwig_cmd = shell.create(
    name="ucsc_wigtobigwig",
    image=WIGTOBIGWIG_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"wig": File, "sizes": File, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        ARGS={inputs.args}
        wigToBigWig $ARGS {inputs.wig} {inputs.sizes} {outputs.results}/{inputs.prefix}.bw
    """,
)

gtftogenepred_cmd = shell.create(
    name="ucsc_gtftogenepred",
    image=GTFTOGENEPRED_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"gtf": File, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        GTF=({inputs.gtf})
        P={inputs.prefix}
        ARGS={inputs.args}
        OUT={outputs.results}
        gtfToGenePred $ARGS "${GTF[0]}" "$OUT/$P.genepred" || exit 1
        if [[ " $ARGS " == *" -genePredExt "* ]] && [[ " $ARGS " == *" -geneNameAsName2 "* ]]; then
            awk 'BEGIN { OFS="\t" } { print $12, $1, $2, $3, $4, $5, $6, $7, $8, $9, $10 }' \
                "$OUT/$P.genepred" > "$OUT/$P.refflat" || exit 1
        fi
    """,
)

# As upstream, the BED's chromosome names must match the bigWig's (UCSC style).
bigwigaverageoverbed_cmd = shell.create(
    name="ucsc_bigwigaverageoverbed",
    image=BIGWIGAVERAGEOVERBED_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"bed": File, "bigwig": File, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        BED=({inputs.bed}); BW=({inputs.bigwig})
        ARGS={inputs.args}
        bigWigAverageOverBed $ARGS "${BW[0]}" "${BED[0]}" {outputs.results}/{inputs.prefix}.tab || exit 1
    """,
)

liftover_cmd = shell.create(
    name="ucsc_liftover",
    image=LIFTOVER_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"bed": File, "chain": File, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        BED=({inputs.bed}); CHAIN=({inputs.chain})
        P={inputs.prefix}
        ARGS={inputs.args}
        OUT={outputs.results}
        liftOver $ARGS "${BED[0]}" "${CHAIN[0]}" "$OUT/$P.lifted.bed" "$OUT/$P.unlifted.bed" || exit 1
    """,
)

# One image per tool, so one env each; `env` aggregates them.
bedclip_env = flyte.TaskEnvironment.from_task("ucsc_bedclip", bedclip_cmd.as_task())
bedgraphtobigwig_env = flyte.TaskEnvironment.from_task("ucsc_bedgraphtobigwig", bedgraphtobigwig_cmd.as_task())
wigtobigwig_env = flyte.TaskEnvironment.from_task("ucsc_wigtobigwig", wigtobigwig_cmd.as_task())
gtftogenepred_env = flyte.TaskEnvironment.from_task("ucsc_gtftogenepred", gtftogenepred_cmd.as_task())
bigwigaverageoverbed_env = flyte.TaskEnvironment.from_task(
    "ucsc_bigwigaverageoverbed", bigwigaverageoverbed_cmd.as_task()
)
liftover_env = flyte.TaskEnvironment.from_task("ucsc_liftover", liftover_cmd.as_task())
env = flyte.TaskEnvironment(
    name="ucsc",
    depends_on=[
        bedclip_env,
        bedgraphtobigwig_env,
        wigtobigwig_env,
        gtftogenepred_env,
        bigwigaverageoverbed_env,
        liftover_env,
    ],
)


async def pick(results: Dir, name: str) -> File:
    found = await results.get_file(name)
    if found is None:
        raise FileNotFoundError(f"no {name}")
    return found


async def bedclip(bedgraph: File, sizes: File, prefix: str) -> File:
    """Clip ``bedgraph`` to the chromosome sizes; returns ``<prefix>.bedGraph``."""
    return await pick(await bedclip_cmd(bedgraph=bedgraph, sizes=sizes, prefix=prefix), f"{prefix}.bedGraph")


async def bedgraphtobigwig(bedgraph: File, sizes: File, prefix: str) -> File:
    """Convert ``bedgraph`` to ``<prefix>.bigWig``."""
    return await pick(await bedgraphtobigwig_cmd(bedgraph=bedgraph, sizes=sizes, prefix=prefix), f"{prefix}.bigWig")


async def wigtobigwig(wig: File, sizes: File, prefix: str, args: str = "") -> File:
    """Convert ``wig`` (plain or gzipped) to ``<prefix>.bw``."""
    return await pick(await wigtobigwig_cmd(wig=wig, sizes=sizes, prefix=prefix, args=args), f"{prefix}.bw")


async def bigwigaverageoverbed(bed: File, bigwig: File, prefix: str, args: str = "") -> File:
    """``<prefix>.tab``: ``bigwig``'s coverage averaged over each region of ``bed``."""
    results = await bigwigaverageoverbed_cmd(bed=bed, bigwig=bigwig, prefix=prefix, args=args)
    return await pick(results, f"{prefix}.tab")


@dataclass
class GtfToGenePredResult:
    genepred: File  # <prefix>.genepred
    refflat: File | None  # <prefix>.refflat, with -genePredExt -geneNameAsName2
    results: Dir


async def gtftogenepred(gtf: File, prefix: str, args: str = "") -> GtfToGenePredResult:
    """Convert ``gtf`` to ``<prefix>.genepred`` (and ``<prefix>.refflat`` when ``args`` asks for gene names)."""
    results = await gtftogenepred_cmd(gtf=gtf, prefix=prefix, args=args)
    flags = args.split()
    refflat = None
    if "-genePredExt" in flags and "-geneNameAsName2" in flags:
        refflat = await pick(results, f"{prefix}.refflat")
    return GtfToGenePredResult(genepred=await pick(results, f"{prefix}.genepred"), refflat=refflat, results=results)


@dataclass
class LiftoverResult:
    lifted: File  # <prefix>.lifted.bed
    unlifted: File  # <prefix>.unlifted.bed: the records liftOver couldn't convert, with the reason


async def liftover(bed: File, chain: File, prefix: str, args: str = "") -> LiftoverResult:
    """Convert ``bed``'s coordinates to another assembly with ``chain`` (plain or gzipped)."""
    results = await liftover_cmd(bed=bed, chain=chain, prefix=prefix, args=args)
    return LiftoverResult(
        lifted=await pick(results, f"{prefix}.lifted.bed"), unlifted=await pick(results, f"{prefix}.unlifted.bed")
    )
