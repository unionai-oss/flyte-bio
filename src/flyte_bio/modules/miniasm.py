"""miniasm — Miniasm, de novo assembly of long reads from their all-vs-all overlaps.

Exposes :func:`miniasm`, which lays out an assembly graph from long reads
and their pairwise overlaps (a PAF, e.g. from ``minimap2 -x ava-ont``).
Miniasm does no consensus, so the contigs keep the reads' error rate.

The outputs keep upstream's names, gzipped with ``gzip -n``:

- ``<prefix>.gfa.gz``: the assembly graph;
- ``<prefix>.fasta.gz``: its segments (``S`` lines) as FASTA, folded to 80
  columns.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

MINIASM_IMAGE = "quay.io/biocontainers/miniasm:0.3_r179--h5bf99c6_2"
# Miniasm is single-threaded.
DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory="8Gi")

miniasm_cmd = shell.create(
    name="miniasm",
    image=MINIASM_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"reads": File, "paf": File, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        READS=({inputs.reads}); PAF=({inputs.paf})
        P={inputs.prefix}
        ARGS={inputs.args}
        OUT={outputs.results}
        W=$(mktemp -d); cd "$W"
        miniasm $ARGS -f "${READS[0]}" "${PAF[0]}" > "$P.gfa"
        awk '/^S/{print ">"$2"\n"$3}' "$P.gfa" | fold > "$P.fasta"
        gzip -n "$P.gfa" "$P.fasta"
        mv "$P.gfa.gz" "$P.fasta.gz" "$OUT"/
    """,
)

env = flyte.TaskEnvironment.from_task("miniasm", miniasm_cmd.as_task())


@dataclass
class MiniasmResult:
    gfa: File  # <prefix>.gfa.gz
    assembly: File  # <prefix>.fasta.gz
    results: Dir


async def miniasm(reads: File, paf: File, prefix: str = "sample", args: str = "") -> MiniasmResult:
    """Assemble ``reads`` from their all-vs-all overlaps ``paf``."""
    results = await miniasm_cmd(reads=reads, paf=paf, prefix=prefix, args=args)
    gfa, fasta = await results.get_file(f"{prefix}.gfa.gz"), await results.get_file(f"{prefix}.fasta.gz")
    if gfa is None or fasta is None:
        raise FileNotFoundError(f"miniasm wrote no {prefix}.gfa.gz or {prefix}.fasta.gz")
    return MiniasmResult(gfa=gfa, assembly=fasta, results=results)
