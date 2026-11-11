"""flye — Flye de novo assembly of long reads (PacBio or Oxford Nanopore).

Exposes :func:`flye`, which assembles one sample's long reads in one of
Flye's read modes (``--pacbio-raw``, ``--pacbio-corr``, ``--pacbio-hifi``,
``--nano-raw``, ``--nano-corr``, ``--nano-hq``). Further Flye options
(``--genome-size``, ``--meta``, ...) go in ``args``.

The outputs keep upstream's names:

- ``<prefix>.assembly.fasta.gz``: the assembly;
- ``<prefix>.assembly_graph.gfa.gz`` and ``<prefix>.assembly_graph.gv.gz``:
  the assembly graph, as GFA and Graphviz;
- ``<prefix>.assembly_info.txt``: per-contig length, coverage and
  circularity;
- ``<prefix>.flye.log`` and ``<prefix>.params.json``: Flye's log and the
  parameters it ran with.

``--threads`` is the task's CPUs.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

FLYE_IMAGE = "community.wave.seqera.io/library/flye:2.9.5--d577924c8416ccd8"

# Upstream's label asks for 12 CPUs and 72 GB. These defaults fit bacterial
# genomes and the test data; raise them with flye_cmd.override(resources=...)
# for larger genomes (a human assembly needs several hundred GB). The 4 CPUs
# match upstream's tests.
DEFAULT_RESOURCES = flyte.Resources(cpu=4, memory="16Gi")

MODES = ("--pacbio-raw", "--pacbio-corr", "--pacbio-hifi", "--nano-raw", "--nano-corr", "--nano-hq")

flye_cmd = shell.create(
    name="flye",
    image=FLYE_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"reads": list[File], "mode": str, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        READS=({{inputs.reads}})
        MODE={{inputs.mode}}
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        flye \
            $MODE \
            "${{READS[@]}}" \
            --out-dir . \
            --threads $TASK_CPUS \
            $ARGS

        gzip -c assembly.fasta > "$OUT/$P.assembly.fasta.gz"
        gzip -c assembly_graph.gfa > "$OUT/$P.assembly_graph.gfa.gz"
        gzip -c assembly_graph.gv > "$OUT/$P.assembly_graph.gv.gz"
        mv assembly_info.txt "$OUT/$P.assembly_info.txt"
        mv flye.log "$OUT/$P.flye.log"
        mv params.json "$OUT/$P.params.json"
    """,
)

env = flyte.TaskEnvironment.from_task("flye", flye_cmd.as_task())


@dataclass
class FlyeResult:
    fasta: File  # the assembly
    gfa: File  # the assembly graph, GFA
    gv: File  # the assembly graph, Graphviz
    txt: File  # assembly_info.txt
    log: File
    json: File  # params.json
    results: Dir  # everything above


async def flye(reads: File | list[File], mode: str, prefix: str = "sample", args: str = "") -> FlyeResult:
    """Assemble ``reads`` (one file or several) in Flye read ``mode``, e.g. ``--nano-hq``."""
    if mode not in MODES:
        raise ValueError(f"Unrecognised mode to run Flye. Options: {', '.join(MODES)}")
    reads = reads if isinstance(reads, list) else [reads]
    if not reads:
        raise ValueError("flye needs reads")
    results = await flye_cmd(reads=reads, mode=mode, prefix=prefix, args=args)

    async def pick(suffix: str) -> File:
        found = await results.get_file(f"{prefix}.{suffix}")
        if found is None:
            raise FileNotFoundError(f"flye wrote no {prefix}.{suffix}")
        return found

    return FlyeResult(
        fasta=await pick("assembly.fasta.gz"),
        gfa=await pick("assembly_graph.gfa.gz"),
        gv=await pick("assembly_graph.gv.gz"),
        txt=await pick("assembly_info.txt"),
        log=await pick("flye.log"),
        json=await pick("params.json"),
        results=results,
    )
