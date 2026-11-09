"""raven — Raven, de novo assembly of long uncorrected reads.

Exposes :func:`raven`, which assembles one sample's long reads (PacBio or
Oxford Nanopore). Further Raven options go in ``args``.

The outputs keep upstream's names:

- ``<prefix>.fasta.gz``: the assembly (polished with Racon, two rounds by
  default);
- ``<prefix>.gfa.gz``: the assembly graph.

``-t`` is the task's CPUs.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

RAVEN_IMAGE = "community.wave.seqera.io/library/raven-assembler:1.6.1--2388e7e21e19cd95"
# Upstream's label asks for 6 CPUs and 36 GB; its tests got 4 CPUs.
DEFAULT_RESOURCES = flyte.Resources(cpu=4, memory="8Gi")

raven_cmd = shell.create(
    name="raven",
    image=RAVEN_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"reads": list[File], "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        set -o pipefail
        READS=({{inputs.reads}})
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        raven -t $TASK_CPUS --graphical-fragment-assembly "$P.gfa" $ARGS "${{READS[@]}}" | gzip -c > "$P.fasta.gz"
        gzip -c "$P.gfa" > "$P.gfa.gz"
        mv "$P.fasta.gz" "$P.gfa.gz" "$OUT"/
    """,
)

env = flyte.TaskEnvironment.from_task("raven", raven_cmd.as_task())


@dataclass
class RavenResult:
    fasta: File  # <prefix>.fasta.gz
    gfa: File  # <prefix>.gfa.gz
    results: Dir


async def raven(reads: File | list[File], prefix: str = "sample", args: str = "") -> RavenResult:
    """Assemble ``reads`` (one file or several)."""
    reads = reads if isinstance(reads, list) else [reads]
    if not reads:
        raise ValueError("raven needs reads")
    results = await raven_cmd(reads=reads, prefix=prefix, args=args)
    fasta, gfa = await results.get_file(f"{prefix}.fasta.gz"), await results.get_file(f"{prefix}.gfa.gz")
    if fasta is None or gfa is None:
        raise FileNotFoundError(f"raven wrote no {prefix}.fasta.gz or {prefix}.gfa.gz")
    return RavenResult(fasta=fasta, gfa=gfa, results=results)
