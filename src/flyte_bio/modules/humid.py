"""humid — UMI-based deduplication of FASTQ reads, without alignment.

Exposes :func:`humid`, which clusters one sample's reads (single- or
paired-end) by their UMIs and sequence with HUMID and returns, as upstream:

- its log, ``<prefix>.log``;
- the deduplicated reads, named by HUMID after the input files with
  ``_dedup`` before the first dot (``test_R1.fastp.fastq.gz`` ->
  ``test_R1_dedup.fastp.fastq.gz``), unless ``args`` turns them off;
- with ``-a`` in ``args``, the reads annotated with their cluster
  (``*_annotated*.fastq.gz``);
- with ``-s`` in ``args``, the statistics folder ``<prefix>/``
  (``stats.dat``, ``neigh.dat``, ``counts.dat``, ``clusters.dat``).

The UMIs are read from the read names (as fastp's ``--umi`` puts them), or
from a separate ``umi_file`` FASTQ.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

HUMID_IMAGE = "quay.io/biocontainers/humid:1.0.4--hadf994f_0"
# Upstream: process_single.
DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory="4Gi")

humid_cmd = shell.create(
    name="humid",
    image=HUMID_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"reads_1": File, "reads_2": File | None, "umi_file": File | None, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        R1=({inputs.reads_1})
        shopt -s nullglob; R2=({inputs.reads_2}); UMIS=({inputs.umi_file}); shopt -u nullglob
        P={inputs.prefix}
        ARGS={inputs.args}
        OUT={outputs.results}
        W=$(mktemp -d); cd "$W"
        humid \
            $ARGS \
            -d "$P" \
            -l "$P.log" \
            "${R1[0]}" "${R2[@]}" \
            "${UMIS[@]}" \
            || exit 1
        shopt -s nullglob
        mv "$P"/*.fastq* "$P.log" "$OUT"/
        STATS=("$P"/*)
        if [ ${#STATS[@]} -gt 0 ]; then mv "$P" "$OUT"/; fi
        true
    """,
)

env = flyte.TaskEnvironment.from_task("humid", humid_cmd.as_task())


@dataclass
class HumidResult:
    log: File  # <prefix>.log
    dedup: list[File]  # the deduplicated reads: one file, or mates 1 and 2
    annotated: list[File]  # with -a: the reads annotated with their cluster
    stats: Dir | None  # with -s: <prefix>/ of stats.dat, neigh.dat, counts.dat, clusters.dat


async def humid(
    reads_1: File,
    reads_2: File | None = None,
    umi_file: File | None = None,
    prefix: str = "sample",
    args: str = "",
) -> HumidResult:
    """Deduplicate one sample's reads by UMI; ``reads_2`` is None for single-end reads."""
    results = await humid_cmd(reads_1=reads_1, reads_2=reads_2, umi_file=umi_file, prefix=prefix, args=args)
    top = sorted(await results.list_files(), key=lambda f: f.name)
    logs = [f for f in top if f.name == f"{prefix}.log"]
    if not logs:
        raise FileNotFoundError(f"humid wrote no {prefix}.log")
    stats = None
    async for file in results.walk():
        if file.name.endswith(".dat"):
            stats = Dir.from_existing_remote(results.path.rstrip("/") + f"/{prefix}")
            break
    return HumidResult(
        log=logs[0],
        dedup=[f for f in top if "_dedup" in f.name and ".fastq" in f.name],
        annotated=[f for f in top if "_annotated" in f.name and ".fastq" in f.name],
        stats=stats,
    )
