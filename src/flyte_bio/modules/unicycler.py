"""unicycler — Unicycler, hybrid or single-technology assembly of bacterial genomes.

Exposes :func:`unicycler`, which assembles one sample from short reads
(single-end, or an R1/R2 pair), long reads, or both. Further Unicycler
options (``--mode``, …) go in ``args``.

The outputs keep upstream's names:

- ``<prefix>.scaffolds.fa.gz``: the assembly;
- ``<prefix>.assembly.gfa.gz``: the assembly graph;
- ``<prefix>.unicycler.log``.

Both assemblies are gzipped with ``gzip -n``. ``--threads`` is the task's
CPUs.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

UNICYCLER_IMAGE = "community.wave.seqera.io/library/unicycler:0.5.1--b9d21c454db1e56b"
# Upstream's label asks for 12 CPUs and 72 GB; these defaults fit bacterial
# genomes. Its tests got 4 CPUs.
DEFAULT_RESOURCES = flyte.Resources(cpu=4, memory="16Gi")

unicycler_cmd = shell.create(
    name="unicycler",
    image=UNICYCLER_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "shortreads_1": File | None,
        "shortreads_2": File | None,
        "longreads": File | None,
        "prefix": str,
        "args": str,
    },
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        shopt -s nullglob
        SHORT=({{inputs.shortreads_1}} {{inputs.shortreads_2}}); LONG=({{inputs.longreads}})
        shopt -u nullglob
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        READS=()
        if [ ${{#SHORT[@]}} -eq 1 ]; then READS+=(-s "${{SHORT[0]}}"); fi
        if [ ${{#SHORT[@]}} -eq 2 ]; then READS+=(-1 "${{SHORT[0]}}" -2 "${{SHORT[1]}}"); fi
        if [ ${{#LONG[@]}} -gt 0 ]; then READS+=(-l "${{LONG[0]}}"); fi
        unicycler \
            --threads $TASK_CPUS \
            $ARGS \
            "${{READS[@]}}" \
            --out ./

        mv assembly.fasta "$P.scaffolds.fa"
        mv assembly.gfa "$P.assembly.gfa"
        gzip -n "$P.scaffolds.fa" "$P.assembly.gfa"
        mv unicycler.log "$P.unicycler.log"
        mv "$P.scaffolds.fa.gz" "$P.assembly.gfa.gz" "$P.unicycler.log" "$OUT"/
    """,
)

env = flyte.TaskEnvironment.from_task("unicycler", unicycler_cmd.as_task())


@dataclass
class UnicyclerResult:
    scaffolds: File  # <prefix>.scaffolds.fa.gz
    gfa: File  # <prefix>.assembly.gfa.gz
    log: File  # <prefix>.unicycler.log
    results: Dir


async def unicycler(
    shortreads: list[File] | None = None, longreads: File | None = None, prefix: str = "sample", args: str = ""
) -> UnicyclerResult:
    """Assemble ``shortreads`` (one file, single-end, or an R1/R2 pair) and/or ``longreads``."""
    shortreads = shortreads or []
    if len(shortreads) > 2:
        raise ValueError(f"unicycler takes one or two short-read files, got {len(shortreads)}")
    if not shortreads and longreads is None:
        raise ValueError("unicycler needs reads")
    results = await unicycler_cmd(
        shortreads_1=shortreads[0] if shortreads else None,
        shortreads_2=shortreads[1] if len(shortreads) == 2 else None,
        longreads=longreads,
        prefix=prefix,
        args=args,
    )
    found = [await results.get_file(f"{prefix}.{s}") for s in ("scaffolds.fa.gz", "assembly.gfa.gz", "unicycler.log")]
    if any(f is None for f in found):
        raise FileNotFoundError(f"unicycler wrote no scaffolds, graph or log for {prefix}")
    scaffolds, gfa, log = found
    return UnicyclerResult(scaffolds=scaffolds, gfa=gfa, log=log, results=results)
