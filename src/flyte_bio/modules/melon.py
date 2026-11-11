"""melon — taxonomic profiling of long metagenomic reads with marker genes.

Exposes :func:`melon`, which profiles one sample's long reads against a
melon database Dir (protein ``prot.dmnd``, nucleotide ``nucl.*.mmi`` indexes
and ``metadata.tsv``), optionally first filtering out reads a Kraken2
database ``k2_db`` assigns (``--db-kraken``). Outputs keep upstream's names:
the profile ``<prefix>.tsv`` and the per-read assignments ``<prefix>.json``,
which melon writes into ``<prefix>/``, and the log ``<prefix>.log``.

``--threads`` is the task's CPUs; options go in ``args``.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

MELON_IMAGE = "quay.io/biocontainers/melon:0.2.5--pyhdfd78af_0"
# Upstream labels it process_medium.
DEFAULT_RESOURCES = flyte.Resources(cpu=4, memory="8Gi")

melon_cmd = shell.create(
    name="melon",
    image=MELON_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"reads": File, "database": Dir, "k2_db": Dir | None, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        READS=({{inputs.reads}})
        DB={{inputs.database}}
        K2_DB={{inputs.k2_db}}
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        # An unset Dir input renders as a path that doesn't exist (or is empty).
        KRAKEN=()
        if [ -d "$K2_DB" ] && [ -n "$(ls -A "$K2_DB")" ]; then KRAKEN=(--db-kraken "$K2_DB"); fi
        melon \
            "${{READS[0]}}" \
            --db "$DB" \
            --output "$P" \
            --threads $TASK_CPUS \
            "${{KRAKEN[@]}}" \
            $ARGS \
            2> "$P.log" \
            || {{ tail -n 20 "$P.log" >&2; exit 1; }}
        cat "$P.log" >&2
        mv "$P".log "$P"/*.tsv "$P"/*.json "$OUT"/
    """,
)

env = flyte.TaskEnvironment.from_task("melon", melon_cmd.as_task())


@dataclass
class MelonResult:
    tsv: File  # <prefix>.tsv, the profile
    json: File  # <prefix>.json, the per-read assignments
    log: File  # <prefix>.log
    results: Dir  # all of the above


async def melon(
    reads: File, database: Dir, k2_db: Dir | None = None, prefix: str = "sample", args: str = ""
) -> MelonResult:
    """Profile one sample's long ``reads`` against the melon ``database`` (after Kraken2 filtering with ``k2_db``)."""
    results = await melon_cmd(reads=reads, database=database, k2_db=k2_db, prefix=prefix, args=args)

    async def pick(name: str) -> File:
        found = await results.get_file(name)
        if found is None:
            raise FileNotFoundError(f"melon wrote no {name}")
        return found

    return MelonResult(
        tsv=await pick(f"{prefix}.tsv"),
        json=await pick(f"{prefix}.json"),
        log=await pick(f"{prefix}.log"),
        results=results,
    )
