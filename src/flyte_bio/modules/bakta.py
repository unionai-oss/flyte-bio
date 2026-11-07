"""bakta — annotation of bacterial genomes, plasmids and MAGs, with its database download.

Exposes:

- :func:`bakta_db_download`, which downloads Bakta's database (``bakta_db
  download``) and returns it as ``bakta_db.tar.gz``, an archive of the
  database folder (``version.json`` and ``bakta.db`` first, then the rest).
  It's a task of its own, cached on its ``args``, so each database downloads
  once. It's an archive rather than a Dir because the database holds hundreds
  of files, too many for a Dir input to stage reliably. The default is the full database (tens of GiB); ``args="--type
  light"`` gets the light one (a few GiB), as upstream's tests do.
- :func:`bakta`, which annotates one assembly FASTA (plain or gzipped)
  against that database, extracting the archive into its working directory
  first. Outputs keep upstream's names: ``<prefix>.embl``,
  ``.faa``, ``.ffn``, ``.fna``, ``.gbff``, ``.gff3``, ``.hypotheticals.tsv``,
  ``.hypotheticals.faa``, ``.tsv``, ``.txt`` and ``.json``. Bakta's other
  outputs (such as its genome plots) are in ``results`` too.

``proteins`` is a trusted protein FASTA or GenBank file (``--proteins``),
``prodigal_tf`` a Prodigal training file (``--prodigal-tf``), ``regions`` a
file of pre-annotated regions (``--regions``) and ``hmms`` a file of trusted
HMMs (``--hmms``).

As upstream, both run with ``HOME`` in the working directory (fontconfig
needs a writable cache directory). The archive keeps the database's
symlinks (its AMRFinderPlus part links ``latest`` to a dated folder).
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

BAKTA_IMAGE = "community.wave.seqera.io/library/bakta_diamond:7830b94718da4f96"
# Both tasks hold the database twice, packed and unpacked: the light database
# (a few GiB unpacked) fits in 20Gi of disk. The full one is tens of GiB, so
# raise both tasks' disk with <cmd>.as_task().override(resources=...).
DOWNLOAD_RESOURCES = flyte.Resources(cpu=1, memory="4Gi", disk="20Gi")
# Upstream's label asks for 6 CPUs and 36 GB; its test runs bakta in 7 GB with
# the light database. Real genomes with the full database need more memory.
DEFAULT_RESOURCES = flyte.Resources(cpu=4, memory="8Gi", disk="20Gi")

HOME_SETUP = r"""
        mkdir -p "$W/home" "$W/home/.cache" "$W/home/.matplotlib"
        export HOME="$W/home" XDG_CACHE_HOME="$W/home/.cache" MPLCONFIGDIR="$W/home/.matplotlib"
""".strip()

bakta_db_download_cmd = shell.create(
    name="bakta_db_download",
    image=BAKTA_IMAGE,
    resources=DOWNLOAD_RESOURCES,
    inputs={"args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        {HOME_SETUP}

        bakta_db download $ARGS

        DB=$(find . -maxdepth 1 -type d -name 'db*' | head -n 1)
        if [ -z "$DB" ]; then echo "bakta_db wrote no db* folder" >&2; exit 1; fi
        cd "$DB"
        REST=()
        for F in * .[!.]*; do
            if [ -e "$F" ] && [ "$F" != version.json ] && [ "$F" != bakta.db ]; then REST+=("$F"); fi
        done
        tar -czf "$OUT/bakta_db.tar.gz" version.json bakta.db "${{REST[@]}}"
    """,
)

bakta_cmd = shell.create(
    name="bakta",
    image=BAKTA_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "fasta": File,
        "db": File,
        "proteins": File | None,
        "prodigal_tf": File | None,
        "regions": File | None,
        "hmms": File | None,
        "prefix": str,
        "args": str,
    },
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        FA=({{inputs.fasta}})
        shopt -s nullglob
        PROTEINS=({{inputs.proteins}}); TF=({{inputs.prodigal_tf}})
        REGIONS=({{inputs.regions}}); HMMS=({{inputs.hmms}})
        shopt -u nullglob
        DB_ARCHIVE=({{inputs.db}})
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        {HOME_SETUP}
        DB="$W/db"; mkdir "$DB"
        tar -xzf "${{DB_ARCHIVE[0]}}" -C "$DB"

        EXTRA=()
        if [ ${{#PROTEINS[@]}} -gt 0 ]; then EXTRA+=(--proteins "${{PROTEINS[0]}}"); fi
        if [ ${{#TF[@]}} -gt 0 ]; then EXTRA+=(--prodigal-tf "${{TF[0]}}"); fi
        if [ ${{#REGIONS[@]}} -gt 0 ]; then EXTRA+=(--regions "${{REGIONS[0]}}"); fi
        if [ ${{#HMMS[@]}} -gt 0 ]; then EXTRA+=(--hmms "${{HMMS[0]}}"); fi

        mkdir out; cd out
        bakta \
            "${{FA[0]}}" \
            $ARGS \
            --threads $TASK_CPUS \
            --prefix "$P" \
            "${{EXTRA[@]}}" \
            --db "$DB"

        mv ./* "$OUT"/
    """,
)

env = flyte.TaskEnvironment.from_task("bakta", bakta_db_download_cmd.as_task(), bakta_cmd.as_task())


async def bakta_db_download(args: str = "") -> File:
    """Download Bakta's database: ``bakta_db.tar.gz``, to pass to :func:`bakta` as ``db``."""
    results = await bakta_db_download_cmd(args=args)
    archive = await results.get_file("bakta_db.tar.gz")
    if archive is None:
        raise FileNotFoundError("bakta_db_download wrote no bakta_db.tar.gz")
    return archive


@dataclass
class BaktaResult:
    embl: File
    faa: File
    ffn: File
    fna: File
    gbff: File
    gff: File  # <prefix>.gff3
    hypotheticals_tsv: File
    hypotheticals_faa: File
    tsv: File
    txt: File
    json: File
    results: Dir  # everything bakta wrote


async def bakta(
    fasta: File,
    db: File,
    proteins: File | None = None,
    prodigal_tf: File | None = None,
    regions: File | None = None,
    hmms: File | None = None,
    prefix: str = "sample",
    args: str = "",
) -> BaktaResult:
    """Annotate the assembly ``fasta`` against ``db``, an archive from :func:`bakta_db_download`."""
    results = await bakta_cmd(
        fasta=fasta,
        db=db,
        proteins=proteins,
        prodigal_tf=prodigal_tf,
        regions=regions,
        hmms=hmms,
        prefix=prefix,
        args=args,
    )

    async def pick(suffix: str) -> File:
        found = await results.get_file(f"{prefix}.{suffix}")
        if found is None:
            raise FileNotFoundError(f"bakta wrote no {prefix}.{suffix}")
        return found

    return BaktaResult(
        embl=await pick("embl"),
        faa=await pick("faa"),
        ffn=await pick("ffn"),
        fna=await pick("fna"),
        gbff=await pick("gbff"),
        gff=await pick("gff3"),
        hypotheticals_tsv=await pick("hypotheticals.tsv"),
        hypotheticals_faa=await pick("hypotheticals.faa"),
        tsv=await pick("tsv"),
        txt=await pick("txt"),
        json=await pick("json"),
        results=results,
    )
