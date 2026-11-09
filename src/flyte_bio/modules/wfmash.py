"""wfmash — all-to-all mapping and base-level alignment of genomes, the first step of pangenome graph construction.

Exposes :func:`wfmash`, which maps (and aligns) the sequences of a FASTA to
themselves or to a list of queries: ``<prefix>.paf``.

As upstream:

- The FASTA may be bgzipped, with its ``.gzi`` and ``.fai`` linked beside it.
- With ``query_self`` the FASTA is also passed as the query, for an
  all-to-all mapping; ``query_file_list`` (a file naming the query sequences)
  adds ``--query-file-list``.
- ``paf`` (existing mappings) adds ``--input-paf``, so wfmash only computes
  the alignments.
- wfmash gets the task's CPUs as threads; options in ``args`` come later.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

WFMASH_IMAGE = "community.wave.seqera.io/library/wfmash:0.13.0--aad0fd683181ae57"

# Upstream labels wfmash process_medium.
WFMASH_RESOURCES = flyte.Resources(cpu=4, memory="8Gi")

wfmash_cmd = shell.create(
    name="wfmash",
    image=WFMASH_IMAGE,
    resources=WFMASH_RESOURCES,
    inputs={
        "fasta": File,
        "gzi": File | None,
        "fai": File | None,
        "paf": File | None,
        "query_file_list": File | None,
        "query_self": bool,
        "prefix": str,
        "args": str,
    },
    defaults={"query_self": True, "args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        set -e
        FA=({{inputs.fasta}})
        shopt -s nullglob
        GZI=({{inputs.gzi}})
        FAI=({{inputs.fai}})
        PAF=({{inputs.paf}})
        QUERY_LIST=({{inputs.query_file_list}})
        shopt -u nullglob
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        for f in "${{FA[@]}}" "${{GZI[@]}}" "${{FAI[@]}}" "${{PAF[@]}}" "${{QUERY_LIST[@]}}"; do
            ln -s "$f" "$(basename "$f")"
        done
        NAME=$(basename "${{FA[0]}}")
        QUERY=()
        if [ {{inputs.query_self}} = true ]; then QUERY=("$NAME"); fi
        if [ ${{#QUERY_LIST[@]}} -gt 0 ]; then QUERY+=(--query-file-list "$(basename "${{QUERY_LIST[0]}}")"); fi
        MAPPINGS=()
        if [ ${{#PAF[@]}} -gt 0 ]; then MAPPINGS=(--input-paf "$(basename "${{PAF[0]}}")"); fi
        wfmash "$NAME" "${{QUERY[@]}}" --threads $TASK_CPUS "${{MAPPINGS[@]}}" $ARGS > "$OUT/$P.paf" || exit 1
    """,
)

env = flyte.TaskEnvironment.from_task("wfmash", wfmash_cmd.as_task())


@dataclass
class WfmashResult:
    paf: File  # <prefix>.paf
    results: Dir


async def wfmash(
    fasta: File,
    gzi: File | None = None,
    fai: File | None = None,
    query_self: bool = True,
    paf: File | None = None,
    query_file_list: File | None = None,
    prefix: str = "sample",
    args: str = "",
) -> WfmashResult:
    """Map (and align) ``fasta``'s sequences, all-to-all with ``query_self``: ``<prefix>.paf``."""
    results = await wfmash_cmd(
        fasta=fasta,
        gzi=gzi,
        fai=fai,
        paf=paf,
        query_file_list=query_file_list,
        query_self=query_self,
        prefix=prefix,
        args=args,
    )
    found = await results.get_file(f"{prefix}.paf")
    if found is None:
        raise FileNotFoundError(f"wfmash wrote no {prefix}.paf")
    return WfmashResult(paf=found, results=results)
