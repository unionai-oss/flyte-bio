"""checkv — quality and completeness of viral genomes with CheckV.

Exposes:

- :func:`checkv_downloaddatabase`, which downloads CheckV's database
  (``checkv download_database``; about 1.5 GB from NERSC at run time) and
  returns it as ``checkv_db.tar.gz``, an archive of the database's
  ``README.txt``, ``genome_db/`` and ``hmm_db/`` (the database holds many
  files, so it's passed as one archive rather than a Dir). It's a task of
  its own, cached on its ``args``, so the database downloads once.
- :func:`checkv_endtoend`, which runs CheckV's whole pipeline
  (``checkv end_to_end``) on a FASTA of viral contigs against a database
  archive, and returns upstream's outputs: ``quality_summary.tsv``,
  ``completeness.tsv``, ``contamination.tsv``, ``complete_genomes.tsv``,
  ``proviruses.fna`` and ``viruses.fna``. ``db`` is a tar archive (any
  compression) of a CheckV database, its files at the top or in a single
  folder, such as :func:`checkv_downloaddatabase` returns.

As upstream, ``-t`` is the task's CPUs.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

CHECKV_IMAGE = "quay.io/biocontainers/checkv:1.0.3--pyhdfd78af_0"
# The download is about 1.5 GB, extracted to some 6 GB, then archived again.
DOWNLOAD_RESOURCES = flyte.Resources(cpu=2, memory="4Gi", disk="40Gi")
# Upstream: process_medium. The full database is some 6 GB extracted.
DEFAULT_RESOURCES = flyte.Resources(cpu=4, memory="16Gi", disk="30Gi")
OUTPUT_NAMES = (
    "quality_summary.tsv",
    "completeness.tsv",
    "contamination.tsv",
    "complete_genomes.tsv",
    "proviruses.fna",
    "viruses.fna",
)
OUTPUT_LIST = " ".join(OUTPUT_NAMES)

checkv_downloaddatabase_cmd = shell.create(
    name="checkv_downloaddatabase",
    image=CHECKV_IMAGE,
    resources=DOWNLOAD_RESOURCES,
    inputs={"args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        ARGS={inputs.args}
        OUT={outputs.results}
        W=$(mktemp -d); cd "$W"
        checkv download_database $ARGS ./checkv_db/ || exit 1
        # The database is a folder (checkv-db-v<version>) inside the destination.
        README=$(find checkv_db -name README.txt | head -n 1)
        if [ -z "$README" ]; then echo "checkv download_database wrote no README.txt" >&2; exit 1; fi
        cd "$(dirname "$README")"
        # The top-level entries first, so the archive's start lists them.
        { echo README.txt; echo genome_db; echo hmm_db; find genome_db hmm_db -mindepth 1; } > "$W/members.txt"
        tar -czf "$OUT/checkv_db.tar.gz" --no-recursion -T "$W/members.txt" || exit 1
    """,
)

checkv_endtoend_cmd = shell.create(
    name="checkv_endtoend",
    image=CHECKV_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"fasta": File, "db": File, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        FA=({{inputs.fasta}}); DB_ARCHIVE=({{inputs.db}})
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        mkdir db
        tar -xf "${{DB_ARCHIVE[0]}}" -C db || exit 1
        DB=db
        # An archive of the database's folder holds it one level down.
        ENTRIES=(db/*)
        if [ ! -e db/genome_db ] && [ ${{#ENTRIES[@]}} -eq 1 ] && [ -d "${{ENTRIES[0]}}" ]; then DB=${{ENTRIES[0]}}; fi
        checkv \
            end_to_end \
            $ARGS \
            -t $TASK_CPUS \
            -d "$DB" \
            "${{FA[0]}}" \
            "$P" \
            || exit 1
        for NAME in {OUTPUT_LIST}; do mv "$P/$NAME" "$OUT"/ || exit 1; done
    """,
)

env = flyte.TaskEnvironment.from_task("checkv", checkv_downloaddatabase_cmd.as_task(), checkv_endtoend_cmd.as_task())


async def checkv_downloaddatabase(args: str = "") -> File:
    """Download CheckV's database: ``checkv_db.tar.gz``, to pass to :func:`checkv_endtoend` as ``db``."""
    results = await checkv_downloaddatabase_cmd(args=args)
    archive = await results.get_file("checkv_db.tar.gz")
    if archive is None:
        raise FileNotFoundError("checkv_downloaddatabase wrote no checkv_db.tar.gz")
    return archive


@dataclass
class CheckvResult:
    quality_summary: File
    completeness: File
    contamination: File
    complete_genomes: File
    proviruses: File
    viruses: File


async def checkv_endtoend(fasta: File, db: File, prefix: str = "sample", args: str = "") -> CheckvResult:
    """Assess the viral contigs in ``fasta`` against the CheckV database archive ``db``."""
    results = await checkv_endtoend_cmd(fasta=fasta, db=db, prefix=prefix, args=args)
    files = {}
    for name in OUTPUT_NAMES:
        found = await results.get_file(name)
        if found is None:
            raise FileNotFoundError(f"checkv end_to_end wrote no {name}")
        files[name.split(".")[0]] = found
    return CheckvResult(**files)
