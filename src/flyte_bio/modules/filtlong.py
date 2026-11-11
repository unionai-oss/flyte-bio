"""filtlong — quality filtering of long reads with Filtlong.

Exposes :func:`filtlong`, which filters one sample's long reads, optionally
scoring them against Illumina reads of the same sample (one single-end file,
``-1``, or an R1/R2 pair, ``-1``/``-2``). Pass the filters (``--min_length``,
``--keep_percent``, ``--target_bases``, ...) in ``args``.

As upstream, every input is staged into the working directory under its
original name, and the outputs keep upstream's names:

- ``<prefix>.fastq.gz``: the kept reads, gzipped with ``gzip -n``;
- ``<prefix>.log``: Filtlong's report (its stderr).

Filtlong is single-threaded.
"""

import os
from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

FILTLONG_IMAGE = "quay.io/biocontainers/filtlong:0.2.1--h9a82719_0"

# Filtlong holds the long reads' scores in memory and, with Illumina reads,
# their k-mers; raise memory with filtlong_cmd.override(resources=...)
# for large read sets.
DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory="8Gi")

filtlong_cmd = shell.create(
    name="filtlong",
    image=FILTLONG_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "long_reads": File,
        "short_reads_1": File | None,
        "short_reads_2": File | None,
        "prefix": str,
        "args": str,
    },
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        LONG=({inputs.long_reads})
        shopt -s nullglob; SHORT=({inputs.short_reads_1} {inputs.short_reads_2}); shopt -u nullglob
        P={inputs.prefix}
        ARGS={inputs.args}
        OUT={outputs.results}
        W=$(mktemp -d); cd "$W"
        for F in "${LONG[@]}" "${SHORT[@]}"; do ln -s "$F" "$(basename "$F")"; done
        SHORT_ARGS=()
        if [ ${#SHORT[@]} -ge 1 ]; then SHORT_ARGS+=(-1 "$(basename "${SHORT[0]}")"); fi
        if [ ${#SHORT[@]} -ge 2 ]; then SHORT_ARGS+=(-2 "$(basename "${SHORT[1]}")"); fi
        set -o pipefail
        filtlong \
            "${SHORT_ARGS[@]}" \
            $ARGS \
            "$(basename "${LONG[0]}")" \
            2> "$P.log" \
            | gzip -n > "$OUT/$P.fastq.gz"
        cat "$P.log" >&2
        mv "$P.log" "$OUT"/
    """,
)

env = flyte.TaskEnvironment.from_task("filtlong", filtlong_cmd.as_task())


@dataclass
class FiltlongResult:
    reads: File  # <prefix>.fastq.gz, the kept long reads
    log: File  # <prefix>.log, Filtlong's report
    results: Dir  # both


async def filtlong(
    long_reads: File,
    short_reads: list[File] | None = None,
    prefix: str = "sample",
    args: str = "",
) -> FiltlongResult:
    """Filter ``long_reads``, scored against ``short_reads`` (one file, or an R1/R2 pair) when given."""
    short_reads = short_reads or []
    if len(short_reads) > 2:
        raise ValueError(f"filtlong takes one or two Illumina files, got {len(short_reads)}")
    if (long_reads.name or os.path.basename(long_reads.path)) == f"{prefix}.fastq.gz":
        raise ValueError(
            "Longread FASTQ input and output names are the same, set prefix in module configuration to disambiguate!"
        )
    results = await filtlong_cmd(
        long_reads=long_reads,
        short_reads_1=short_reads[0] if short_reads else None,
        short_reads_2=short_reads[1] if len(short_reads) == 2 else None,
        prefix=prefix,
        args=args,
    )

    async def pick(name: str) -> File:
        found = await results.get_file(name)
        if found is None:
            raise FileNotFoundError(f"filtlong wrote no {name}")
        return found

    return FiltlongResult(reads=await pick(f"{prefix}.fastq.gz"), log=await pick(f"{prefix}.log"), results=results)
