"""wipertools — repair malformed FASTQ files, in parallel chunks.

Exposes four steps, in the order a repair runs them:

- :func:`wipertools_fastqscatter` splits a FASTQ into ``num_splits`` chunks,
  ``<prefix>_000.fastq.gz``, ``<prefix>_001.fastq.gz``…, returned in order.
- :func:`wipertools_fastqwiper` repairs one FASTQ (or chunk): it drops
  malformed lines and reads (e.g. a quality string whose length differs from
  its sequence), writing ``<prefix>.fastq.gz`` and a report of what it
  dropped, ``<prefix>.report``. Upstream's prefix is ``<sample>_wiped``.
- :func:`wipertools_fastqgather` joins chunks, in name order as upstream,
  into ``<prefix>.fastq.gz`` (upstream's prefix: ``<sample>_gather``).
- :func:`wipertools_reportgather` sums fastqwiper reports into
  ``<prefix>.report`` (upstream's prefix: ``<sample>_gather``).

As upstream, fastqscatter writes into the folder ``args`` names with
``-o``/``--out_folder``, else ``chunks`` (added to ``args``); the chunks are
returned from it.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

WIPERTOOLS_IMAGE = "community.wave.seqera.io/library/wipertools:1.1.4--042a9b3f81640596"
# Upstream labels every step process_single (1 CPU, 6 GB).
DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory="6Gi")

wipertools_fastqscatter_cmd = shell.create(
    name="wipertools_fastqscatter",
    image=WIPERTOOLS_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"fastq": File, "num_splits": int, "out_folder": str, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        FQ=({inputs.fastq})
        N={inputs.num_splits}
        FOLDER={inputs.out_folder}
        P={inputs.prefix}
        ARGS={inputs.args}
        OUT={outputs.results}
        W=$(mktemp -d); cd "$W"
        wipertools fastqscatter -f "${FQ[0]}" -n "$N" -p "$P" $ARGS || exit 1
        mv "$FOLDER"/* "$OUT"/
    """,
)

wipertools_fastqwiper_cmd = shell.create(
    name="wipertools_fastqwiper",
    image=WIPERTOOLS_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"fastq": File, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        FQ=({inputs.fastq})
        P={inputs.prefix}
        ARGS={inputs.args}
        OUT={outputs.results}
        W=$(mktemp -d); cd "$W"
        wipertools fastqwiper -i "${FQ[0]}" -o "$OUT/$P.fastq.gz" -r "$OUT/$P.report" $ARGS || exit 1
    """,
)

# A list[File] is staged under the files' names; as upstream, they're passed
# sorted by name (byte order).
wipertools_fastqgather_cmd = shell.create(
    name="wipertools_fastqgather",
    image=WIPERTOOLS_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"fastqs": list[File], "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        FQ=({inputs.fastqs})
        P={inputs.prefix}
        ARGS={inputs.args}
        OUT={outputs.results}
        W=$(mktemp -d); cd "$W"
        for F in "${FQ[@]}"; do ln -s "$F" "$(basename "$F")"; done
        mapfile -t NAMES < <(for F in "${FQ[@]}"; do basename "$F"; done | LC_ALL=C sort)
        wipertools fastqgather -i "${NAMES[@]}" -o "$P.fastq.gz" $ARGS || exit 1
        mv "$P.fastq.gz" "$OUT"/
    """,
)

wipertools_reportgather_cmd = shell.create(
    name="wipertools_reportgather",
    image=WIPERTOOLS_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"reports": list[File], "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        REPORTS=({inputs.reports})
        P={inputs.prefix}
        ARGS={inputs.args}
        OUT={outputs.results}
        W=$(mktemp -d); cd "$W"
        for F in "${REPORTS[@]}"; do ln -s "$F" "$(basename "$F")"; done
        mapfile -t NAMES < <(for F in "${REPORTS[@]}"; do basename "$F"; done | LC_ALL=C sort)
        wipertools reportgather -r "${NAMES[@]}" -f "$P.report" $ARGS || exit 1
        mv "$P.report" "$OUT"/
    """,
)

env = flyte.TaskEnvironment.from_task(
    "wipertools",
    wipertools_fastqscatter_cmd.as_task(),
    wipertools_fastqwiper_cmd.as_task(),
    wipertools_fastqgather_cmd.as_task(),
    wipertools_reportgather_cmd.as_task(),
)


def file_name(file: File) -> str:
    return file.path.rstrip("/").rsplit("/", 1)[-1]


async def pick(results: Dir, name: str) -> File:
    found = await results.get_file(name)
    if found is None:
        raise FileNotFoundError(f"wipertools wrote no {name}")
    return found


@dataclass
class ScatterResult:
    chunks: list[File]  # <prefix>_000.fastq.gz, <prefix>_001.fastq.gz, … in order
    results: Dir


async def wipertools_fastqscatter(
    fastq: File, num_splits: int, prefix: str = "sample", args: str = ""
) -> ScatterResult:
    """Split ``fastq`` into ``num_splits`` chunks."""
    if num_splits < 1:
        raise ValueError(f"fastqscatter needs at least one split, got {num_splits}")
    tokens = args.split()
    out_folder = "chunks"
    for flag in ("--out_folder", "-o"):
        if flag in tokens:
            index = tokens.index(flag)
            if index + 1 >= len(tokens):
                raise ValueError(f"fastqscatter args end with {flag} but name no folder")
            out_folder = tokens[index + 1]
            break
    else:
        args = f"{args} -o {out_folder}".strip()
    results = await wipertools_fastqscatter_cmd(
        fastq=fastq, num_splits=num_splits, out_folder=out_folder, prefix=prefix, args=args
    )
    chunks = sorted([f async for f in results.walk()], key=file_name)
    if not chunks:
        raise FileNotFoundError("fastqscatter wrote no chunks")
    return ScatterResult(chunks=chunks, results=results)


@dataclass
class WiperResult:
    fastq: File  # <prefix>.fastq.gz
    report: File  # <prefix>.report
    results: Dir


async def wipertools_fastqwiper(fastq: File, prefix: str = "sample_wiped", args: str = "") -> WiperResult:
    """Repair ``fastq``, dropping malformed reads; returns the repaired reads and a report."""
    if file_name(fastq) == f"{prefix}.fastq.gz":
        raise ValueError(f"fastqwiper input {file_name(fastq)} has the output's name; choose another prefix")
    results = await wipertools_fastqwiper_cmd(fastq=fastq, prefix=prefix, args=args)
    return WiperResult(
        fastq=await pick(results, f"{prefix}.fastq.gz"),
        report=await pick(results, f"{prefix}.report"),
        results=results,
    )


def check_gather_inputs(files: list[File], output: str, label: str) -> None:
    if not files:
        raise ValueError(f"{label} needs at least one file")
    names = [file_name(f) for f in files]
    if len(set(names)) != len(names):
        raise ValueError(f"{label} inputs need distinct names, got {names}")
    if output in names:
        raise ValueError(f"{label} output {output} matches an input's name; choose another prefix")


async def wipertools_fastqgather(fastqs: list[File], prefix: str = "sample_gather", args: str = "") -> File:
    """Join FASTQ chunks, in name order, into ``<prefix>.fastq.gz``."""
    check_gather_inputs(fastqs, f"{prefix}.fastq.gz", "fastqgather")
    results = await wipertools_fastqgather_cmd(fastqs=fastqs, prefix=prefix, args=args)
    return await pick(results, f"{prefix}.fastq.gz")


async def wipertools_reportgather(reports: list[File], prefix: str = "sample_gather", args: str = "") -> File:
    """Sum fastqwiper ``reports`` into ``<prefix>.report``."""
    check_gather_inputs(reports, f"{prefix}.report", "reportgather")
    results = await wipertools_reportgather_cmd(reports=reports, prefix=prefix, args=args)
    return await pick(results, f"{prefix}.report")
