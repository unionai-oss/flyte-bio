"""nanolyse — NanoLyse: remove lambda phage reads from Oxford Nanopore FASTQ.

Exposes :func:`nanolyse`, which streams one sample's gzipped reads through
NanoLyse against a reference FASTA (usually the lambda phage DNA control
strand) and, as upstream, writes the reads that don't map to
``<prefix>.fastq.gz`` and NanoLyse's log to ``<prefix>.nanolyse.log``.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

NANOLYSE_IMAGE = "quay.io/biocontainers/nanolyse:1.2.0--py_0"
# Upstream labels NanoLyse process_low (2 CPUs, 12 GB); it maps with one
# minimap2 thread, and test-sized runs need little memory.
NANOLYSE_RESOURCES = flyte.Resources(cpu=2, memory="4Gi")

nanolyse_cmd = shell.create(
    name="nanolyse",
    image=NANOLYSE_IMAGE,
    resources=NANOLYSE_RESOURCES,
    inputs={"fastq": File, "fasta": File, "prefix": str},
    outputs={"results": Dir},
    script=r"""
        FQ=({inputs.fastq}); FA=({inputs.fasta})
        P={inputs.prefix}
        OUT={outputs.results}
        W=$(mktemp -d); cd "$W"
        FA_NAME=$(basename "${FA[0]}"); ln -s "${FA[0]}" "$FA_NAME"
        gunzip -c "${FQ[0]}" | NanoLyse -r "$FA_NAME" | gzip > "$OUT/$P.fastq.gz" || exit 1
        mv NanoLyse.log "$OUT/$P.nanolyse.log"
    """,
)

env = flyte.TaskEnvironment.from_task("nanolyse", nanolyse_cmd.as_task())


@dataclass
class NanoLyseResult:
    reads: File  # <prefix>.fastq.gz: the reads that don't map to the reference
    log: File  # <prefix>.nanolyse.log


async def nanolyse(fastq: File, fasta: File, prefix: str) -> NanoLyseResult:
    """Remove the reads of gzipped ``fastq`` that map to ``fasta`` (e.g. lambda phage)."""
    results = await nanolyse_cmd(fastq=fastq, fasta=fasta, prefix=prefix)

    async def pick(name: str) -> File:
        found = await results.get_file(name)
        if found is None:
            raise FileNotFoundError(f"nanolyse wrote no {name}")
        return found

    return NanoLyseResult(reads=await pick(f"{prefix}.fastq.gz"), log=await pick(f"{prefix}.nanolyse.log"))
