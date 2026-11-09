"""autocycler — Autocycler consensus assembly of long-read bacterial genomes.

Autocycler builds a consensus from several assemblies of the same genome. Its
subcommands chain: subsample the reads, assemble each subsample with other
tools, then compress → cluster → trim (per cluster) → resolve (per cluster)
→ combine. Exposes:

- :func:`autocycler_subsample` — split long reads into subsets for separate
  assemblies: ``sample_01.fastq.gz``, ``sample_02.fastq.gz``, …. As upstream,
  ``--seed 42`` is added unless ``args`` sets a seed.
- :func:`autocycler_compress` — compress the input assemblies into a unitig
  graph: ``input_assemblies.gfa`` and ``input_assemblies.yaml``.
- :func:`autocycler_cluster` — cluster the graph's contigs:
  ``clustering.{newick,tsv,yaml}``, ``pairwise_distances.phylip`` and, per
  passing cluster, ``qc_pass/cluster_NNN/1_untrimmed.{gfa,yaml}``.
- :func:`autocycler_trim` — trim one cluster's overlaps:
  ``2_trimmed.{gfa,yaml}``.
- :func:`autocycler_resolve` — resolve one trimmed cluster to a consensus:
  ``<prefix>_3_bridged.gfa``, ``<prefix>_4_merged.gfa`` and
  ``<prefix>_5_final.gfa``.
- :func:`autocycler_combine` — combine resolved clusters into the final
  assembly: ``consensus_assembly.{fasta,gfa,yaml}``.

``compress`` and ``trim`` get the task's CPUs as threads. Autocycler exits
255 on error, which the platform would retry forever, so its failures exit 1.
"""

import asyncio
from dataclasses import dataclass
from pathlib import PurePosixPath

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

AUTOCYCLER_IMAGE = "community.wave.seqera.io/library/autocycler:0.5.2--91b58d5b8fccdfb9"
# Upstream labels subsample, compress, cluster and trim process_low and
# resolve and combine process_single.
LOW_RESOURCES = flyte.Resources(cpu=2, memory="8Gi")
SINGLE_RESOURCES = flyte.Resources(cpu=1, memory="4Gi")

autocycler_subsample_cmd = shell.create(
    name="autocycler_subsample",
    image=AUTOCYCLER_IMAGE,
    resources=LOW_RESOURCES,
    inputs={"reads": File, "genome_size": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        READS=({inputs.reads})
        GENOME_SIZE={inputs.genome_size}
        ARGS={inputs.args}
        OUT={outputs.results}
        if [[ " $ARGS " != *--seed* ]]; then ARGS="$ARGS --seed 42"; fi
        W=$(mktemp -d); cd "$W"
        autocycler subsample $ARGS --reads "${READS[0]}" --out_dir subsampled --genome_size "$GENOME_SIZE" || exit 1
        for F in subsampled/*.fastq; do
            gzip -n "$F"
            mv "$F.gz" "$OUT"/
        done
    """,
)

autocycler_compress_cmd = shell.create(
    name="autocycler_compress",
    image=AUTOCYCLER_IMAGE,
    resources=LOW_RESOURCES,
    inputs={"assemblies": list[File], "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        ASSEMBLIES=({{inputs.assemblies}})
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        # As upstream, the assemblies sit in the working directory under their own names.
        for F in "${{ASSEMBLIES[@]}}"; do ln -s "$F" "$(basename "$F")"; done
        autocycler compress $ARGS -t $TASK_CPUS -i . -a compress || exit 1
        mv compress/*.gfa compress/*.yaml "$OUT"/
    """,
)

autocycler_cluster_cmd = shell.create(
    name="autocycler_cluster",
    image=AUTOCYCLER_IMAGE,
    resources=LOW_RESOURCES,
    inputs={"gfa": File, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        GFA=({inputs.gfa})
        ARGS={inputs.args}
        OUT={outputs.results}
        W=$(mktemp -d); cd "$W"
        ln -s "${GFA[0]}" input_assemblies.gfa
        autocycler cluster $ARGS -a . || exit 1
        cp clustering/*.newick clustering/*.tsv clustering/*.phylip clustering/*.yaml "$OUT"/
        for C in clustering/qc_pass/*/; do
            NAME=$(basename "$C")
            mkdir -p "$OUT/qc_pass/$NAME"
            cp "$C"*.gfa "$C"*.yaml "$OUT/qc_pass/$NAME"/
        done
    """,
)

autocycler_trim_cmd = shell.create(
    name="autocycler_trim",
    image=AUTOCYCLER_IMAGE,
    resources=LOW_RESOURCES,
    inputs={"gfa": File, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        GFA=({{inputs.gfa}})
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        ln -s "${{GFA[0]}}" 1_untrimmed.gfa
        autocycler trim $ARGS --threads $TASK_CPUS -c . || exit 1
        mv 2_trimmed.gfa 2_trimmed.yaml "$OUT"/
    """,
)

autocycler_resolve_cmd = shell.create(
    name="autocycler_resolve",
    image=AUTOCYCLER_IMAGE,
    resources=SINGLE_RESOURCES,
    inputs={"gfa": File, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        GFA=({inputs.gfa})
        P={inputs.prefix}
        ARGS={inputs.args}
        OUT={outputs.results}
        W=$(mktemp -d); cd "$W"
        ln -s "${GFA[0]}" 2_trimmed.gfa
        autocycler resolve $ARGS -c . || exit 1
        mv 3_bridged.gfa "$OUT/${P}_3_bridged.gfa"
        mv 4_merged.gfa "$OUT/${P}_4_merged.gfa"
        mv 5_final.gfa "$OUT/${P}_5_final.gfa"
    """,
)

autocycler_combine_cmd = shell.create(
    name="autocycler_combine",
    image=AUTOCYCLER_IMAGE,
    resources=SINGLE_RESOURCES,
    inputs={"clusters": list[File], "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        CLUSTERS=({inputs.clusters})
        ARGS={inputs.args}
        OUT={outputs.results}
        W=$(mktemp -d); cd "$W"
        NAMES=()
        for F in "${CLUSTERS[@]}"; do ln -s "$F" "$(basename "$F")"; NAMES+=("$(basename "$F")"); done
        autocycler combine $ARGS -i "${NAMES[@]}" -a combine || exit 1
        mv combine/consensus_assembly.fasta combine/consensus_assembly.gfa combine/consensus_assembly.yaml "$OUT"/
    """,
)

env = flyte.TaskEnvironment.from_task(
    "autocycler",
    autocycler_subsample_cmd.as_task(),
    autocycler_compress_cmd.as_task(),
    autocycler_cluster_cmd.as_task(),
    autocycler_trim_cmd.as_task(),
    autocycler_resolve_cmd.as_task(),
    autocycler_combine_cmd.as_task(),
)


async def pick(results: Dir, name: str, tool: str) -> File:
    found = await results.get_file(name)
    if found is None:
        raise FileNotFoundError(f"autocycler {tool} wrote no {name}")
    return found


async def autocycler_subsample(reads: File, genome_size: str | int, args: str = "") -> list[File]:
    """Subsample ``reads`` for ``genome_size`` (e.g. ``5m``): ``sample_NN.fastq.gz``, sorted by name."""
    results = await autocycler_subsample_cmd(reads=reads, genome_size=str(genome_size), args=args)
    files = [f async for f in results.walk() if f.path.endswith(".fastq.gz")]
    if not files:
        raise FileNotFoundError("autocycler subsample wrote no reads")
    return sorted(files, key=lambda f: PurePosixPath(f.path).name)


@dataclass
class AutocyclerCompressResult:
    gfa: File  # input_assemblies.gfa
    stats: File  # input_assemblies.yaml
    results: Dir


async def autocycler_compress(assemblies: list[File], args: str = "") -> AutocyclerCompressResult:
    """Compress ``assemblies`` (FASTAs of the same genome) into Autocycler's unitig graph."""
    if not assemblies:
        raise ValueError("autocycler compress needs assemblies")
    results = await autocycler_compress_cmd(assemblies=assemblies, args=args)
    return AutocyclerCompressResult(
        gfa=await pick(results, "input_assemblies.gfa", "compress"),
        stats=await pick(results, "input_assemblies.yaml", "compress"),
        results=results,
    )


@dataclass
class AutocyclerClusterResult:
    clusters: dict[str, File]  # cluster name (``cluster_001``, …) → 1_untrimmed.gfa, QC-passing clusters only
    cluster_stats: dict[str, File]  # cluster name → 1_untrimmed.yaml
    newick: File  # clustering.newick
    tsv: File  # clustering.tsv
    pairwise_distances: File  # pairwise_distances.phylip
    stats: File  # clustering.yaml
    results: Dir


async def autocycler_cluster(gfa: File, args: str = "") -> AutocyclerClusterResult:
    """Cluster the contigs of ``gfa`` (``autocycler compress``'s graph)."""
    results = await autocycler_cluster_cmd(gfa=gfa, args=args)
    clusters: dict[str, File] = {}
    cluster_stats: dict[str, File] = {}
    async for f in results.walk():
        path = PurePosixPath(f.path)
        if path.parent.parent.name != "qc_pass":
            continue
        if path.suffix == ".gfa":
            clusters[path.parent.name] = f
        elif path.suffix == ".yaml":
            cluster_stats[path.parent.name] = f
    newick, tsv, phylip, stats = await asyncio.gather(
        pick(results, "clustering.newick", "cluster"),
        pick(results, "clustering.tsv", "cluster"),
        pick(results, "pairwise_distances.phylip", "cluster"),
        pick(results, "clustering.yaml", "cluster"),
    )
    return AutocyclerClusterResult(
        clusters=dict(sorted(clusters.items())),
        cluster_stats=dict(sorted(cluster_stats.items())),
        newick=newick,
        tsv=tsv,
        pairwise_distances=phylip,
        stats=stats,
        results=results,
    )


@dataclass
class AutocyclerTrimResult:
    gfa: File  # 2_trimmed.gfa
    stats: File  # 2_trimmed.yaml
    results: Dir


async def autocycler_trim(gfa: File, args: str = "") -> AutocyclerTrimResult:
    """Trim one cluster's ``1_untrimmed.gfa``."""
    results = await autocycler_trim_cmd(gfa=gfa, args=args)
    return AutocyclerTrimResult(
        gfa=await pick(results, "2_trimmed.gfa", "trim"),
        stats=await pick(results, "2_trimmed.yaml", "trim"),
        results=results,
    )


@dataclass
class AutocyclerResolveResult:
    bridged: File  # <prefix>_3_bridged.gfa
    merged: File  # <prefix>_4_merged.gfa
    resolved: File  # <prefix>_5_final.gfa
    results: Dir


async def autocycler_resolve(gfa: File, prefix: str = "sample", args: str = "") -> AutocyclerResolveResult:
    """Resolve one cluster's ``2_trimmed.gfa`` to a consensus."""
    results = await autocycler_resolve_cmd(gfa=gfa, prefix=prefix, args=args)
    bridged, merged, resolved = await asyncio.gather(
        pick(results, f"{prefix}_3_bridged.gfa", "resolve"),
        pick(results, f"{prefix}_4_merged.gfa", "resolve"),
        pick(results, f"{prefix}_5_final.gfa", "resolve"),
    )
    return AutocyclerResolveResult(bridged=bridged, merged=merged, resolved=resolved, results=results)


@dataclass
class AutocyclerCombineResult:
    fasta: File  # consensus_assembly.fasta
    gfa: File  # consensus_assembly.gfa
    stats: File  # consensus_assembly.yaml
    results: Dir


async def autocycler_combine(clusters: File | list[File], args: str = "") -> AutocyclerCombineResult:
    """Combine resolved clusters (``autocycler resolve``'s final graphs) into one assembly."""
    clusters = clusters if isinstance(clusters, list) else [clusters]
    if not clusters:
        raise ValueError("autocycler combine needs clusters")
    results = await autocycler_combine_cmd(clusters=clusters, args=args)
    fasta, gfa, stats = await asyncio.gather(
        pick(results, "consensus_assembly.fasta", "combine"),
        pick(results, "consensus_assembly.gfa", "combine"),
        pick(results, "consensus_assembly.yaml", "combine"),
    )
    return AutocyclerCombineResult(fasta=fasta, gfa=gfa, stats=stats, results=results)
