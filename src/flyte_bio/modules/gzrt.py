"""gzrt — recover the readable data from corrupt gzip files.

Exposes :func:`gzrt`, which runs ``gzrecover`` over one sample's gzipped
FASTQ (single-end) or two (paired-end), recompressing what it recovers. As
upstream, it writes ``<prefix>.fastq.gz`` or ``<prefix>_1.fastq.gz`` /
``<prefix>_2.fastq.gz``, and an output with nothing recovered holds one empty
line. Upstream's prefix is ``<sample>_recovered``.

gzrecover's exit status is ignored, as upstream (its pipe into gzip keeps only
gzip's): it reports the corruption it skips over.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

GZRT_IMAGE = "quay.io/biocontainers/gzrt:0.9.1--h577a1d6_1"
# Upstream labels gzrt process_single (1 CPU, 6 GB).
DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory="6Gi")

gzrt_cmd = shell.create(
    name="gzrt",
    image=GZRT_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"reads_1": File, "reads_2": File | None, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        R1=({inputs.reads_1})
        shopt -s nullglob; R2=({inputs.reads_2}); shopt -u nullglob
        P={inputs.prefix}
        ARGS={inputs.args}
        OUT={outputs.results}
        W=$(mktemp -d); cd "$W"
        recover() {
            (set +o pipefail; gzrecover $ARGS -p "$1" | gzip > "$2") || exit 1
            if [ ! -s "$2" ]; then echo "" | gzip > "$2" || exit 1; fi
        }
        if [ ${#R2[@]} -eq 0 ]; then
            recover "${R1[0]}" "$OUT/$P.fastq.gz"
        else
            recover "${R1[0]}" "$OUT/${P}_1.fastq.gz"
            recover "${R2[0]}" "$OUT/${P}_2.fastq.gz"
        fi
    """,
)

env = flyte.TaskEnvironment.from_task("gzrt", gzrt_cmd.as_task())


@dataclass
class GzrtResult:
    reads: list[File]  # <prefix>.fastq.gz, or <prefix>_1.fastq.gz and <prefix>_2.fastq.gz
    results: Dir


def file_name(file: File) -> str:
    return file.path.rstrip("/").rsplit("/", 1)[-1]


async def gzrt(
    reads_1: File, reads_2: File | None = None, prefix: str = "sample_recovered", args: str = ""
) -> GzrtResult:
    """Recover what's readable of one sample's gzipped FASTQs; ``reads_2`` is None for single-end."""
    inputs = [reads_1] if reads_2 is None else [reads_1, reads_2]
    names = [f"{prefix}.fastq.gz"] if reads_2 is None else [f"{prefix}_1.fastq.gz", f"{prefix}_2.fastq.gz"]
    for read in inputs:
        if not file_name(read).endswith(".gz"):
            raise ValueError(f"gzrt works with .gz files only, got {file_name(read)}")
        if file_name(read) in names:
            raise ValueError(f"gzrt input {file_name(read)} has an output's name; choose another prefix")
    results = await gzrt_cmd(reads_1=reads_1, reads_2=reads_2, prefix=prefix, args=args)
    reads = []
    for name in names:
        found = await results.get_file(name)
        if found is None:
            raise FileNotFoundError(f"gzrt wrote no {name}")
        reads.append(found)
    return GzrtResult(reads=reads, results=results)
