"""rundbcan — carbohydrate-active enzyme (CAZyme) annotation with run_dbcan.

Exposes:

- :func:`rundbcan_database`, which downloads the dbCAN database
  (``run_dbcan database --aws_s3``) and returns it as ``dbcan_db.tar.gz``, an
  archive of the ``dbcan_db/`` folder (``CAZy.dmnd``, ``dbCAN.hmm``,
  ``dbCAN-sub.hmm``, ``dbCAN-PUL/``, ...). It's a task of its own, so the
  database downloads once and is cached. It's an archive rather than a Dir
  because ``dbCAN-PUL/`` holds hundreds of files, more than a task's input
  staging copes with.
- :func:`rundbcan_cazymeannotation`, ``run_dbcan CAZyme_annotation`` of a
  protein FASTA (may be gzipped) against that archive (or any ``.tar.gz`` of
  a dbCAN database folder), extracted in the task. Outputs keep upstream's
  names: ``<prefix>_overview.tsv`` (the per-protein summary of the three
  methods), ``<prefix>_dbCAN_hmm_results.tsv``,
  ``<prefix>_dbCANsub_hmm_results.tsv`` and ``<prefix>_diamond.out``.

The annotation runs in protein mode with ``--threads`` of the task's CPUs.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

DBCAN_IMAGE = "community.wave.seqera.io/library/dbcan:5.2.9--fb040f67fa175a9b"
# Upstream: database process_low, annotation process_medium. The database is a
# few GiB; both tasks hold it packed and unpacked.
DATABASE_RESOURCES = flyte.Resources(cpu=1, memory="4Gi", disk="40Gi")
DEFAULT_RESOURCES = flyte.Resources(cpu=4, memory="16Gi", disk="40Gi")

rundbcan_database_cmd = shell.create(
    name="rundbcan_database",
    image=DBCAN_IMAGE,
    resources=DATABASE_RESOURCES,
    inputs={"args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        ARGS={inputs.args}
        OUT={outputs.results}
        W=$(mktemp -d); cd "$W"
        run_dbcan database \
            --db_dir dbcan_db \
            --aws_s3 \
            $ARGS || exit 1
        tar -czf "$OUT/dbcan_db.tar.gz" dbcan_db
    """,
)

rundbcan_cazymeannotation_cmd = shell.create(
    name="rundbcan_cazymeannotation",
    image=DBCAN_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"fasta": File, "db": File, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        FA=({{inputs.fasta}})
        DB_ARCHIVE=({{inputs.db}})
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        mkdir db
        tar -xzf "${{DB_ARCHIVE[0]}}" -C db || exit 1
        MARKER=$(find db -name CAZy.dmnd | head -n 1)
        if [ -z "$MARKER" ]; then echo "the dbCAN archive holds no CAZy.dmnd" >&2; exit 1; fi
        DB=$(dirname "$MARKER")
        mkdir out
        run_dbcan CAZyme_annotation \
            --mode protein \
            --db_dir "$DB" \
            --input_raw_data "${{FA[0]}}" \
            --output_dir out \
            --threads $TASK_CPUS \
            $ARGS || exit 1
        mv out/overview.tsv "$OUT/${{P}}_overview.tsv"
        mv out/dbCAN_hmm_results.tsv "$OUT/${{P}}_dbCAN_hmm_results.tsv"
        mv out/dbCANsub_hmm_results.tsv "$OUT/${{P}}_dbCANsub_hmm_results.tsv"
        mv out/diamond.out "$OUT/${{P}}_diamond.out"
    """,
)

env = flyte.TaskEnvironment.from_task(
    "rundbcan", rundbcan_database_cmd.as_task(), rundbcan_cazymeannotation_cmd.as_task()
)


async def rundbcan_database(args: str = "") -> File:
    """Download the dbCAN database: ``dbcan_db.tar.gz``, to pass to :func:`rundbcan_cazymeannotation`."""
    results = await rundbcan_database_cmd(args=args)
    archive = await results.get_file("dbcan_db.tar.gz")
    if archive is None:
        raise FileNotFoundError("run_dbcan database wrote no dbcan_db.tar.gz")
    return archive


@dataclass
class RundbcanResult:
    cazyme_annotation: File  # <prefix>_overview.tsv
    dbcanhmm_results: File  # <prefix>_dbCAN_hmm_results.tsv
    dbcansub_results: File  # <prefix>_dbCANsub_hmm_results.tsv
    dbcandiamond_results: File  # <prefix>_diamond.out


async def rundbcan_cazymeannotation(fasta: File, db: File, prefix: str = "sample", args: str = "") -> RundbcanResult:
    """CAZymes in the proteins of ``fasta``, against the dbCAN database archive ``db``."""
    results = await rundbcan_cazymeannotation_cmd(fasta=fasta, db=db, prefix=prefix, args=args)

    async def picked(name: str) -> File:
        found = await results.get_file(name)
        if found is None:
            raise FileNotFoundError(f"run_dbcan CAZyme_annotation wrote no {name}")
        return found

    return RundbcanResult(
        cazyme_annotation=await picked(f"{prefix}_overview.tsv"),
        dbcanhmm_results=await picked(f"{prefix}_dbCAN_hmm_results.tsv"),
        dbcansub_results=await picked(f"{prefix}_dbCANsub_hmm_results.tsv"),
        dbcandiamond_results=await picked(f"{prefix}_diamond.out"),
    )
