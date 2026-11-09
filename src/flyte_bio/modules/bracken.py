"""bracken — species (or other rank) abundance re-estimation from a Kraken2 report.

Exposes :func:`bracken`, which runs Bracken on one sample's Kraken2 report
(``--report`` of :func:`flyte_bio.modules.kraken2.kraken2`) against ``db``, a
Dir holding a Kraken2 database with Bracken's ``database<read length>mers.kmer_distrib``
files. Outputs keep upstream's names:

- ``<prefix>.tsv``: Bracken's abundance table (``-o``);
- ``<prefix>.kraken2.report_bracken.txt``: the Kraken2-style report with the
  re-estimated counts (``-w``).

The rank (``-l``), read length (``-r``) and threshold (``-t``) go in
``args``, as upstream's ``ext.args``; Bracken's defaults are ``-l S -t 10
-r 100``.

:func:`bracken_build` adds Bracken's k-mer distribution files to a Kraken2
database (``bracken-build``). It takes the uncleaned database archive of
:func:`flyte_bio.modules.kraken2.kraken2_build` (``clean=False``) and, as
upstream, keeps its ``*.k2d`` and ``*.map`` files, ``library/added/`` and
``taxonomy/``. It returns the built database's top-level files
(``database.kraken``, ``database<read length>mers.kmer_distrib``,
``database<read length>mers.kraken``, the ``.k2d`` files and
``seqid2taxid.map``): a Dir :func:`bracken` and
:func:`flyte_bio.modules.kraken2.kraken2` read. The whole database
directory, library and taxonomy included, also comes back as
``bracken-database.tar.gz``. The k-mer length (``-k``) and read length
(``-l``) go in ``args``; ``-t`` is the task's CPUs.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

BRACKEN_IMAGE = "community.wave.seqera.io/library/bracken:3.1--22a4e66ce04c5e01"
DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory="2Gi")

bracken_cmd = shell.create(
    name="bracken",
    image=BRACKEN_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"kraken_report": File, "db": Dir, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        REPORT=({inputs.kraken_report})
        DB={inputs.db}
        P={inputs.prefix}
        ARGS={inputs.args}
        OUT={outputs.results}
        W=$(mktemp -d); cd "$W"
        # Bracken's estimation step is Python and may exit 255 on errors, which
        # the cluster takes for an interruption; any failure exits 1 instead.
        bracken \
            $ARGS \
            -d "$DB" \
            -i "${REPORT[0]}" \
            -o "$P.tsv" \
            -w "$P.kraken2.report_bracken.txt" \
            || exit 1
        mv "$P.tsv" "$P.kraken2.report_bracken.txt" "$OUT"/
    """,
)

# Upstream labels it process_high.
BUILD_RESOURCES = flyte.Resources(cpu=4, memory="16Gi")

bracken_build_cmd = shell.create(
    name="bracken_build",
    image=BRACKEN_IMAGE,
    resources=BUILD_RESOURCES,
    inputs={"kraken2_db": File, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        DB=({{inputs.kraken2_db}})
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"

        mkdir kraken2-database
        tar -xzf "${{DB[0]}}" -C kraken2-database --strip-components 1
        mkdir -p bracken-database/library/added bracken-database/taxonomy
        shopt -s nullglob
        for F in kraken2-database/*k2d kraken2-database/*map; do mv "$F" bracken-database/; done
        for F in kraken2-database/library/added/*; do mv "$F" bracken-database/library/added/; done
        for F in kraken2-database/taxonomy/*; do mv "$F" bracken-database/taxonomy/; done
        shopt -u nullglob
        rm -rf kraken2-database

        bracken-build \
            $ARGS \
            -t $TASK_CPUS \
            -d bracken-database/ \
            || exit 1
        tar -czf "$OUT/bracken-database.tar.gz" bracken-database
        mkdir "$OUT/bracken-database"
        find bracken-database -maxdepth 1 -type f -exec cp {{}} "$OUT/bracken-database/" \;
    """,
)

env = flyte.TaskEnvironment.from_task("bracken", bracken_cmd.as_task(), bracken_build_cmd.as_task())


@dataclass
class BrackenResult:
    report: File  # <prefix>.tsv, the abundance table
    kraken_style_report: File  # <prefix>.kraken2.report_bracken.txt
    results: Dir  # both of the above


async def bracken(kraken_report: File, db: Dir, prefix: str = "sample", args: str = "") -> BrackenResult:
    """Re-estimate the abundances in ``kraken_report`` with Bracken against ``db``."""
    results = await bracken_cmd(kraken_report=kraken_report, db=db, prefix=prefix, args=args)

    async def pick(name: str) -> File:
        found = await results.get_file(name)
        if found is None:
            raise FileNotFoundError(f"bracken wrote no {name}")
        return found

    return BrackenResult(
        report=await pick(f"{prefix}.tsv"),
        kraken_style_report=await pick(f"{prefix}.kraken2.report_bracken.txt"),
        results=results,
    )


@dataclass
class BrackenBuildResult:
    db: Dir  # the database's top-level files: Kraken2's and Bracken's
    archive: File  # the whole database directory, bracken-database.tar.gz


async def bracken_build(kraken2_db: File, args: str = "") -> BrackenBuildResult:
    """Add Bracken's files to an uncleaned Kraken2 database archive (``kraken2_build(clean=False)``)."""
    results = await bracken_build_cmd(kraken2_db=kraken2_db, args=args)
    archive = await results.get_file("bracken-database.tar.gz")
    if archive is None:
        raise FileNotFoundError("bracken_build wrote no bracken-database.tar.gz")
    return BrackenBuildResult(
        db=Dir.from_existing_remote(results.path.rstrip("/") + "/bracken-database"), archive=archive
    )
