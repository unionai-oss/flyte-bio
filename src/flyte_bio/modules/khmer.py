"""khmer — k-mer counting with khmer.

Exposes :func:`khmer_uniquekmers`, which estimates the number of unique
k-mers of length ``kmer_size`` in a FASTA (with ``unique-kmers.py``, a
HyperLogLog counter) and returns, under upstream's names:

- khmer's report, ``<prefix>.report.txt``;
- the count alone, ``<prefix>.kmers.txt``, as pipelines use it for MACS's
  effective genome size.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

KHMER_IMAGE = "quay.io/biocontainers/khmer:3.0.0a3--py37haa7609a_2"
DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory="4Gi")

# As upstream: the count is the "number of unique k-mers" line of the report.
khmer_uniquekmers_cmd = shell.create(
    name="khmer_uniquekmers",
    image=KHMER_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"fasta": File, "kmer_size": int, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        FA=({inputs.fasta})
        P={inputs.prefix}
        ARGS={inputs.args}
        OUT={outputs.results}
        W=$(mktemp -d); cd "$W"
        ln -s "${FA[0]}" "$(basename "${FA[0]}")"
        unique-kmers.py -k {inputs.kmer_size} -R "$OUT/$P.report.txt" $ARGS "$(basename "${FA[0]}")"
        grep ^number "$OUT/$P.report.txt" | sed 's/^.*:.[[:blank:]]//g' > "$OUT/$P.kmers.txt"
    """,
)

env = flyte.TaskEnvironment.from_task("khmer", khmer_uniquekmers_cmd.as_task())


@dataclass
class KhmerUniquekmersResult:
    report: File  # <prefix>.report.txt
    kmers: File  # <prefix>.kmers.txt: the number of unique k-mers
    results: Dir  # everything the task wrote


async def khmer_uniquekmers(
    fasta: File, kmer_size: int, prefix: str = "sample", args: str = ""
) -> KhmerUniquekmersResult:
    """Estimate the number of unique ``kmer_size``-mers in ``fasta``."""
    results = await khmer_uniquekmers_cmd(fasta=fasta, kmer_size=kmer_size, prefix=prefix, args=args)

    async def pick(name: str) -> File:
        file = await results.get_file(name)
        if file is None:
            raise FileNotFoundError(f"khmer wrote no {name}")
        return file

    return KhmerUniquekmersResult(
        report=await pick(f"{prefix}.report.txt"), kmers=await pick(f"{prefix}.kmers.txt"), results=results
    )
