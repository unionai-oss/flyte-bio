"""fcsgx — NCBI Foreign Contamination Screen: cross-species contamination (FCS-GX).

Exposes:

- :func:`fcsgx_fetchdb` — download a GX database by its manifest URL
  (``sync_files.py get``), checking each file's md5, as
  ``<prefix>.tar.gz`` of the directory ``<prefix>/`` (default prefix
  ``gxdb_<manifest name>``). The full database is about 470 GB (its manifest
  is https://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/FCS/database/latest/all.manifest);
  the test database (``test-only``) is 4.6 GB.
- :func:`fcsgx_rungx` — screen an assembly against a database from
  :func:`fcsgx_fetchdb` for sequences from organisms other than the given
  NCBI taxid (``run_gx.py``): ``<prefix>…fcs_gx_report.txt`` (the action
  report), ``<prefix>…taxonomy.rpt``, ``<prefix>…summary.txt`` (the log) and,
  when there were any, the hits ``<prefix>…hits.tsv.gz``.
- :func:`fcsgx_cleangenome` — apply an action report to the assembly
  (``gx clean-genome``): ``<prefix>.cleaned.fasta`` and the removed
  ``<prefix>.contaminants.fasta``.

Databases travel as one ``.tar.gz`` (a few very large files), extracted on the
task's disk; :func:`fcsgx_rungx` reads it from there with ``GX_NUM_CORES`` set
to the task's CPUs. Upstream can first copy the database to a RAM disk, which
makes the run much faster on a real database; that option isn't ported. A
real database needs the task's memory and disk raised well beyond the
defaults (NCBI recommends 512 GB of memory), with
``fcsgx_rungx_cmd.override(resources=...)``.

FCS-GX and its databases are US government works, in the public domain.
"""

from dataclasses import dataclass
from pathlib import PurePosixPath
from urllib.parse import urlsplit

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

FCSGX_IMAGE = "quay.io/biocontainers/ncbi-fcs-gx:0.5.5--h9948957_0"
# Upstream: fetchdb process_single, rungx process_high, cleangenome process_low.
# These fit the 4.6 GB test database.
FETCHDB_RESOURCES = flyte.Resources(cpu=1, memory="4Gi", disk="30Gi")
RUNGX_RESOURCES = flyte.Resources(cpu=8, memory="24Gi", disk="30Gi")
CLEANGENOME_RESOURCES = flyte.Resources(cpu=1, memory="4Gi")

fcsgx_fetchdb_cmd = shell.create(
    name="fcsgx_fetchdb",
    image=FCSGX_IMAGE,
    resources=FETCHDB_RESOURCES,
    inputs={"manifest": str, "prefix": str},
    outputs={"results": Dir},
    script=r"""
        set -o pipefail
        MANIFEST={inputs.manifest}
        P={inputs.prefix}
        W=$(mktemp -d); cd "$W"
        sync_files.py get --mft "$MANIFEST" --dir "$P" || exit 1
        tar -cf - "$P" | gzip -c > {outputs.results}/"$P.tar.gz" || exit 1
    """,
)

fcsgx_rungx_cmd = shell.create(
    name="fcsgx_rungx",
    image=FCSGX_IMAGE,
    resources=RUNGX_RESOURCES,
    inputs={"fasta": File, "taxid": str, "gxdb": File, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        FA=({{inputs.fasta}})
        DB=({{inputs.gxdb}})
        TAXID={{inputs.taxid}}
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"

        # The database: the archive's single top-level directory, or its root.
        mkdir gxdb
        tar -xzf "${{DB[0]}}" -C gxdb || exit 1
        DATABASE=gxdb
        TOP=(gxdb/*)
        if [ ${{#TOP[@]}} -eq 1 ] && [ -d "${{TOP[0]}}" ]; then DATABASE=${{TOP[0]}}; fi

        mkdir run; cd run
        ln -s "${{FA[0]}}" "$(basename "${{FA[0]}}")"
        export GX_NUM_CORES=$TASK_CPUS
        run_gx.py \
            --fasta "$(basename "${{FA[0]}}")" \
            --gx-db "$W/$DATABASE" \
            --tax-id "$TAXID" \
            --generate-logfile true \
            --out-basename "$P" \
            --out-dir . \
            $ARGS || exit 1
        shopt -s nullglob
        mv ./*.fcs_gx_report.txt ./*.taxonomy.rpt ./*.summary.txt ./*.hits.tsv.gz "$OUT"/
    """,
)

fcsgx_cleangenome_cmd = shell.create(
    name="fcsgx_cleangenome",
    image=FCSGX_IMAGE,
    resources=CLEANGENOME_RESOURCES,
    inputs={"fasta": File, "fcsgx_report": File, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        FA=({inputs.fasta})
        REPORT=({inputs.fcsgx_report})
        P={inputs.prefix}
        ARGS={inputs.args}
        OUT={outputs.results}
        gx clean-genome \
            --input "${FA[0]}" \
            --action-report "${REPORT[0]}" \
            --output "$OUT/$P.cleaned.fasta" \
            --contam-fasta-out "$OUT/$P.contaminants.fasta" \
            $ARGS || exit 1
    """,
)

env = flyte.TaskEnvironment.from_task(
    "fcsgx", fcsgx_fetchdb_cmd.as_task(), fcsgx_rungx_cmd.as_task(), fcsgx_cleangenome_cmd.as_task()
)


def file_name(file: File) -> str:
    return file.path.rstrip("/").rsplit("/", 1)[-1]


async def fcsgx_fetchdb(manifest: str, prefix: str = "") -> File:
    """The GX database of the manifest at URL ``manifest``: ``<prefix>.tar.gz`` of ``<prefix>/``.

    ``prefix`` defaults, as upstream, to ``gxdb_<manifest name without extension>``.
    """
    if not prefix:
        prefix = "gxdb_" + PurePosixPath(urlsplit(manifest).path).stem
    results = await fcsgx_fetchdb_cmd(manifest=manifest, prefix=prefix)
    found = await results.get_file(f"{prefix}.tar.gz")
    if found is None:
        raise FileNotFoundError(f"fcsgx fetchdb wrote no {prefix}.tar.gz")
    return found


@dataclass
class FcsgxRungxResult:
    fcsgx_report: File  # *.fcs_gx_report.txt: the action report
    taxonomy_report: File  # *.taxonomy.rpt
    log: File  # *.summary.txt
    hits: File | None  # *.hits.tsv.gz, when there were any
    results: Dir  # everything run_gx.py kept


async def fcsgx_rungx(fasta: File, taxid: str, gxdb: File, prefix: str = "sample", args: str = "") -> FcsgxRungxResult:
    """Screen ``fasta`` (from NCBI taxid ``taxid``) against the GX database archive ``gxdb``."""
    results = await fcsgx_rungx_cmd(fasta=fasta, taxid=taxid, gxdb=gxdb, prefix=prefix, args=args)
    files = await results.list_files()

    def ending(suffix: str) -> File | None:
        found = [f for f in files if file_name(f).endswith(suffix)]
        return found[0] if found else None

    picked = {}
    for suffix in (".fcs_gx_report.txt", ".taxonomy.rpt", ".summary.txt"):
        picked[suffix] = ending(suffix)
        if picked[suffix] is None:
            raise FileNotFoundError(f"fcsgx run_gx.py wrote no *{suffix}")
    return FcsgxRungxResult(
        fcsgx_report=picked[".fcs_gx_report.txt"],
        taxonomy_report=picked[".taxonomy.rpt"],
        log=picked[".summary.txt"],
        hits=ending(".hits.tsv.gz"),
        results=results,
    )


@dataclass
class FcsgxCleangenomeResult:
    cleaned: File  # <prefix>.cleaned.fasta
    contaminants: File  # <prefix>.contaminants.fasta


async def fcsgx_cleangenome(
    fasta: File, fcsgx_report: File, prefix: str = "sample", args: str = ""
) -> FcsgxCleangenomeResult:
    """``fasta`` with the actions of ``fcsgx_report`` applied, and the sequences it removed."""
    results = await fcsgx_cleangenome_cmd(fasta=fasta, fcsgx_report=fcsgx_report, prefix=prefix, args=args)

    async def pick(name: str) -> File:
        found = await results.get_file(name)
        if found is None:
            raise FileNotFoundError(f"gx clean-genome wrote no {name}")
        return found

    return FcsgxCleangenomeResult(
        cleaned=await pick(f"{prefix}.cleaned.fasta"), contaminants=await pick(f"{prefix}.contaminants.fasta")
    )
