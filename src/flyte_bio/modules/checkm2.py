"""checkm2 — CheckM2, completeness and contamination of genomes and metagenome-assembled genomes.

Exposes:

- :func:`checkm2_databasedownload` — download CheckM2's DIAMOND database
  from Zenodo (record ``14897628``, version 3, by default):
  ``checkm2_db_v<version>.dmnd``. As upstream, the record's version and md5
  come from Zenodo's API, and ``aria2c`` checks the download against the md5.
- :func:`checkm2_predict` — predict the quality of one or more genomes
  (plain or gzipped FASTA) with a CheckM2 database: CheckM2's output
  directory ``<prefix>/`` and its report ``<prefix>_checkm2_report.tsv``.

As upstream, the database goes to ``checkm2 predict`` through ``CHECKM2DB``
rather than ``--database_path``, which skips CheckM2's checksum check, so
reduced databases (like the test one) work. ``--threads`` follows the task's
CPUs; ``args`` come last.
"""

import asyncio
import json
import urllib.request
from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

CHECKM2_DOWNLOAD_IMAGE = "community.wave.seqera.io/library/aria2:1.37.0--3a9ec328469995dd"
CHECKM2_IMAGE = "community.wave.seqera.io/library/checkm2:1.1.0--60f287bc25d7a10d"

# Upstream: process_single. The archive is about 1.7 GB, the database about 3 GB.
DOWNLOAD_RESOURCES = flyte.Resources(cpu=1, memory="2Gi", disk="16Gi")
# Upstream: process_medium. The 4 CPUs match upstream's test CPUs; the full database wants about 16 GiB
# (less with --lowmem).
PREDICT_RESOURCES = flyte.Resources(cpu=4, memory="16Gi")

# Version 3 of CheckM2's database, upstream's default.
DEFAULT_ZENODO_ID = 14897628

checkm2_databasedownload_cmd = shell.create(
    name="checkm2_databasedownload",
    image=CHECKM2_DOWNLOAD_IMAGE,
    resources=DOWNLOAD_RESOURCES,
    inputs={"zenodo_id": int, "db_version": str, "checksum": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        ID={inputs.zenodo_id}
        VERSION={inputs.db_version}
        CHECKSUM={inputs.checksum}
        ARGS={inputs.args}
        W=$(mktemp -d); cd "$W"
        # Automatic download is broken when using singularity/apptainer (https://github.com/chklovski/CheckM2/issues/73)
        # So it's necessary to download the database manually
        aria2c \
            $ARGS \
            --checksum "$CHECKSUM" \
            "https://zenodo.org/records/$ID/files/checkm2_database.tar.gz" || exit 1

        tar -xzf checkm2_database.tar.gz || exit 1
        rm -f checkm2_database.tar.gz
        DB_PATH=$(find . -name '*.dmnd')
        mv "$DB_PATH" {outputs.results}/"checkm2_db_v$VERSION.dmnd"
    """,
)

checkm2_predict_cmd = shell.create(
    name="checkm2_predict",
    image=CHECKM2_IMAGE,
    resources=PREDICT_RESOURCES,
    inputs={"fasta": list[File], "db": File, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        FA=({{inputs.fasta}}); DB=({{inputs.db}})
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        mkdir input_bins
        INPUTS=()
        for F in "${{FA[@]}}"; do
            ln -s "$F" input_bins/"$(basename "$F")"
            INPUTS+=(input_bins/"$(basename "$F")")
        done

        # Passing the database via CHECKM2DB instead of --database_path skips CheckM2's
        # hardcoded checksum check, allowing custom/reduced databases (e.g. for testing)
        export CHECKM2DB="${{DB[0]}}"

        checkm2 \
            predict \
            --input "${{INPUTS[@]}}" \
            --output-directory "$P" \
            --threads $TASK_CPUS \
            $ARGS || exit 1

        cp "$P"/quality_report.tsv "${{P}}_checkm2_report.tsv"
        mv "$P" "${{P}}_checkm2_report.tsv" "$OUT"/
    """,
)

# One image per environment: the download runs in an aria2 image, prediction in CheckM2's.
checkm2_download_env = flyte.TaskEnvironment.from_task(
    "checkm2_databasedownload", checkm2_databasedownload_cmd.as_task()
)
checkm2_predict_env = flyte.TaskEnvironment.from_task("checkm2_predict", checkm2_predict_cmd.as_task())
env = flyte.TaskEnvironment(name="checkm2", depends_on=[checkm2_download_env, checkm2_predict_env])


def zenodo_record(zenodo_id: int) -> dict:
    with urllib.request.urlopen(f"https://zenodo.org/api/records/{zenodo_id}", timeout=60) as response:
        return json.load(response)


@dataclass
class Checkm2Database:
    database: File  # checkm2_db_v<version>.dmnd
    version: str


async def checkm2_databasedownload(zenodo_id: int | None = None, args: str = "") -> Checkm2Database:
    """Download CheckM2's database from Zenodo record ``zenodo_id`` (default: version 3)."""
    record_id = zenodo_id or DEFAULT_ZENODO_ID
    record = await asyncio.to_thread(zenodo_record, record_id)
    version = str(record["metadata"]["version"])
    checksum = record["files"][0]["checksum"].replace("md5:", "md5=", 1)
    # As upstream, a user agent is added unless args set one.
    if "--user-agent" not in args:
        args = f"{args} --user-agent=Wget/1.21.4".strip()
    results = await checkm2_databasedownload_cmd(zenodo_id=record_id, db_version=version, checksum=checksum, args=args)
    name = f"checkm2_db_v{version}.dmnd"
    database = await results.get_file(name)
    if database is None:
        raise FileNotFoundError(f"CheckM2 database download wrote no {name}")
    return Checkm2Database(database=database, version=version)


@dataclass
class Checkm2Result:
    report: File  # <prefix>_checkm2_report.tsv, a copy of quality_report.tsv
    results: Dir  # the report and CheckM2's output directory <prefix>/


async def checkm2_predict(fasta: File | list[File], db: File, prefix: str = "sample", args: str = "") -> Checkm2Result:
    """Predict completeness and contamination of each genome in ``fasta`` with the CheckM2 database ``db``."""
    results = await checkm2_predict_cmd(
        fasta=fasta if isinstance(fasta, list) else [fasta], db=db, prefix=prefix, args=args
    )
    name = f"{prefix}_checkm2_report.tsv"
    report = await results.get_file(name)
    if report is None:
        raise FileNotFoundError(f"CheckM2 wrote no {name}")
    return Checkm2Result(report=report, results=results)
