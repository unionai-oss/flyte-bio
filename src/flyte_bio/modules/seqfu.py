"""seqfu — sequence-length statistics for FASTA/FASTQ files.

Exposes :func:`seqfu_stats`, ``seqfu stats`` over one or more FASTA/FASTQ
files (gzipped or not), one row per file in the caller's order. As upstream,
the files are passed under their own names, which seqfu records in the
tables' first column. Outputs keep upstream's names:

- the statistics table: ``<prefix>.tsv``;
- the same statistics as a MultiQC custom-content table:
  ``<prefix>_mqc.txt``.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

SEQFU_IMAGE = "quay.io/biocontainers/seqfu:1.22.3--hfd12232_2"
# Upstream's process_single.
DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory="2Gi")

# A list[File] is staged under the files' names and globbed alphabetically, so
# ``names`` carries the caller's order (one name per line) into the script.
seqfu_stats_cmd = shell.create(
    name="seqfu_stats",
    image=SEQFU_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"files": list[File], "names": str, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        IN=({inputs.files})
        NAMES={inputs.names}
        P={inputs.prefix}
        ARGS={inputs.args}
        OUT={outputs.results}
        set -f
        W=$(mktemp -d); cd "$W"
        D=$(dirname "${IN[0]}")
        FILES=()
        while IFS= read -r NAME; do
            ln -s "$D/$NAME" "$NAME"
            FILES+=("$NAME")
        done <<< "$NAMES"
        seqfu stats $ARGS --multiqc "$OUT/${P}_mqc.txt" "${FILES[@]}" > "$OUT/$P.tsv" || exit 1
    """,
)

env = flyte.TaskEnvironment.from_task("seqfu", seqfu_stats_cmd.as_task())


@dataclass
class SeqfuStatsResult:
    stats: File  # <prefix>.tsv
    multiqc: File  # <prefix>_mqc.txt


def file_name(file: File) -> str:
    return file.path.rstrip("/").rsplit("/", 1)[-1]


async def seqfu_stats(files: list[File], prefix: str = "sample", args: str = "") -> SeqfuStatsResult:
    """``seqfu stats`` of ``files`` (one row per file, in list order)."""
    if not files:
        raise ValueError("seqfu_stats needs at least one file")
    names = [file_name(f) for f in files]
    if len(set(names)) != len(names):
        raise ValueError(f"seqfu_stats inputs need distinct names, got {names}")
    results = await seqfu_stats_cmd(files=files, names="\n".join(names), prefix=prefix, args=args)

    async def pick(name: str) -> File:
        found = await results.get_file(name)
        if found is None:
            raise FileNotFoundError(f"seqfu stats wrote no {name}")
        return found

    return SeqfuStatsResult(stats=await pick(f"{prefix}.tsv"), multiqc=await pick(f"{prefix}_mqc.txt"))
