"""dragonflye — Dragonflye, long-read bacterial genome assembly and polishing.

Exposes :func:`dragonflye`, which assembles one sample's long reads with
Flye (the default), Miniasm or Raven (``--assembler`` in ``args``), polishes
the assembly with the long reads and, given an R1/R2 pair of short reads,
with those too. The genome size (``--gsize``) and other options go in
``args``.

The outputs keep upstream's names:

- ``<prefix>.fa`` (and ``<prefix>.reoriented.fa`` when Dragonflye
  reorients the contigs): the final contigs;
- ``{flye,miniasm,raven}.fasta``: the assembler's raw contigs;
- ``{flye,miniasm,raven}-unpolished.gfa``: its assembly graph, when written;
- ``flye-info.txt``: Flye's contig info (Flye only);
- ``dragonflye.log``.

``--cpus`` and ``--ram`` (whole GiB) follow the task's CPUs and memory.
"""

from dataclasses import dataclass
from pathlib import PurePosixPath

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

DRAGONFLYE_IMAGE = "quay.io/biocontainers/dragonflye:1.2.1--hdfd78af_0"
# Upstream's label asks for 6 CPUs and 36 GB; its tests got 4 CPUs.
DEFAULT_RESOURCES = flyte.Resources(cpu=4, memory="16Gi")

dragonflye_cmd = shell.create(
    name="dragonflye",
    image=DRAGONFLYE_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"longreads": File, "shortreads_1": File | None, "shortreads_2": File | None, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        LONG=({{inputs.longreads}})
        shopt -s nullglob; SHORT=({{inputs.shortreads_1}} {{inputs.shortreads_2}}); shopt -u nullglob
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        # As upstream, the reads sit in the working directory under their own names.
        for F in "${{LONG[@]}}" "${{SHORT[@]}}"; do ln -s "$F" "$(basename "$F")"; done
        POLISH=()
        if [ ${{#SHORT[@]}} -eq 2 ]; then
            POLISH=(--R1 "$(basename "${{SHORT[0]}}")" --R2 "$(basename "${{SHORT[1]}}")")
        fi
        dragonflye \
            --reads "$(basename "${{LONG[0]}}")" \
            "${{POLISH[@]}}" \
            $ARGS \
            --prefix "$P" \
            --cpus $TASK_CPUS \
            --ram $(( TASK_MEMORY_MB / 1024 )) \
            --outdir ./ \
            --force

        shopt -s nullglob
        for F in *.fa dragonflye.log {{flye,miniasm,raven}}.fasta {{flye,miniasm,raven}}-unpolished.gfa \
                flye-info.txt; do
            if [ -e "$F" ]; then mv "$F" "$OUT"/; fi
        done
    """,
)

env = flyte.TaskEnvironment.from_task("dragonflye", dragonflye_cmd.as_task())


@dataclass
class DragonflyeResult:
    contigs: list[File]  # *.fa: <prefix>.fa and, when written, <prefix>.reoriented.fa
    log: File  # dragonflye.log
    raw_contigs: File  # {flye,miniasm,raven}.fasta
    gfa: File | None  # {flye,miniasm,raven}-unpolished.gfa
    txt: File | None  # flye-info.txt
    results: Dir  # everything above


async def dragonflye(
    longreads: File, shortreads: list[File] | None = None, prefix: str = "sample", args: str = ""
) -> DragonflyeResult:
    """Assemble ``longreads``, polishing with an optional R1/R2 pair of ``shortreads``."""
    shortreads = shortreads or []
    if len(shortreads) not in (0, 2):
        raise ValueError(f"dragonflye polishes with an R1/R2 pair of short reads, got {len(shortreads)} files")
    results = await dragonflye_cmd(
        longreads=longreads,
        shortreads_1=shortreads[0] if shortreads else None,
        shortreads_2=shortreads[1] if shortreads else None,
        prefix=prefix,
        args=args,
    )
    return await read_dragonflye_results(results)


async def read_dragonflye_results(results: Dir) -> DragonflyeResult:
    """Pick the outputs out of ``dragonflye_cmd``'s results Dir (e.g. from an overridden task)."""
    files = {PurePosixPath(f.path).name: f async for f in results.walk()}
    contigs = [files[name] for name in sorted(files) if name.endswith(".fa")]
    raw = [files[f"{a}.fasta"] for a in ("flye", "miniasm", "raven") if f"{a}.fasta" in files]
    gfa = [files[f"{a}-unpolished.gfa"] for a in ("flye", "miniasm", "raven") if f"{a}-unpolished.gfa" in files]
    if not contigs or len(raw) != 1 or "dragonflye.log" not in files:
        raise FileNotFoundError(f"dragonflye wrote {sorted(files)}: missing contigs, raw contigs or log")
    return DragonflyeResult(
        contigs=contigs,
        log=files["dragonflye.log"],
        raw_contigs=raw[0],
        gfa=gfa[0] if gfa else None,
        txt=files.get("flye-info.txt"),
        results=results,
    )
