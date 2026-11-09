"""mifaser — mi-faser, functional annotation of metagenomic reads with enzyme (EC) numbers.

Exposes :func:`mifaser`, which maps one sample's reads (one FASTQ, ``-f``, or
several, ``-l``) to a mi-faser database with DIAMOND and writes, under
upstream's names:

- ``<prefix>_multi_ec.tsv``: reads matching several EC numbers;
- ``<prefix>_analysis.tsv``: the analysis summary;
- ``<prefix>_ec_count.tsv``: the read count per EC number.

The database (e.g. ``GS-24-all``: ``database.dmnd``, ``ec_mapping.tsv``,
``ec_annotation.tsv``, ``sequences.fasta``) is a tar archive (any compression),
extracted in the task. When its entries share one top-level directory, that
directory is the database, as upstream's untar strips it. As upstream,
mi-faser runs one worker thread (``--threads 1``) with the task's CPUs
(``--cpu``).
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

MIFASER_IMAGE = "quay.io/biocontainers/mifaser:1.64--pyh106432d_0"
# Upstream labels it process_medium.
DEFAULT_RESOURCES = flyte.Resources(cpu=4, memory="8Gi")

mifaser_cmd = shell.create(
    name="mifaser",
    image=MIFASER_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"reads": list[File], "database": File, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        READS=({{inputs.reads}})
        DATABASE=({{inputs.database}})
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        mkdir database
        python -c 'import sys, tarfile; tarfile.open(sys.argv[1]).extractall(sys.argv[2])' \
            "${{DATABASE[0]}}" database || exit 1
        ENTRIES=(database/*)
        DB_DIR=database
        if [ ${{#ENTRIES[@]}} -eq 1 ] && [ -d "${{ENTRIES[0]}}" ]; then DB_DIR="${{ENTRIES[0]}}"; fi
        INPUT_FLAG=-l
        if [ ${{#READS[@]}} -eq 1 ]; then INPUT_FLAG=-f; fi

        mifaser \
            $ARGS \
            $INPUT_FLAG "${{READS[@]}}" \
            --threads 1 \
            --cpu $TASK_CPUS \
            --databasefolder "$W/$DB_DIR" \
            --outputfolder "mifaser-$P/" || exit 1

        for SUFFIX in multi_ec.tsv analysis.tsv ec_count.tsv; do
            mv "mifaser-$P/$SUFFIX" "$OUT/${{P}}_$SUFFIX"
        done
    """,
)

env = flyte.TaskEnvironment.from_task("mifaser", mifaser_cmd.as_task())


@dataclass
class MifaserResult:
    multi_ec: File  # <prefix>_multi_ec.tsv
    analysis: File  # <prefix>_analysis.tsv
    ec_counts: File  # <prefix>_ec_count.tsv
    results: Dir  # everything above


async def mifaser(reads: list[File], database: File, prefix: str = "sample", args: str = "") -> MifaserResult:
    """EC-number annotation of one sample's ``reads`` against a mi-faser ``database`` archive."""
    results = await mifaser_cmd(reads=reads, database=database, prefix=prefix, args=args)

    async def pick(suffix: str) -> File:
        found = await results.get_file(f"{prefix}_{suffix}")
        if found is None:
            raise FileNotFoundError(f"mi-faser wrote no {prefix}_{suffix}")
        return found

    return MifaserResult(
        multi_ec=await pick("multi_ec.tsv"),
        analysis=await pick("analysis.tsv"),
        ec_counts=await pick("ec_count.tsv"),
        results=results,
    )
