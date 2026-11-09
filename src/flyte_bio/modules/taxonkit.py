"""taxonkit — NCBI taxonomy lookups with TaxonKit.

Exposes :func:`taxonkit_lineage`, which writes the full lineage of taxon IDs
(``taxonkit lineage``) as ``<prefix>.tsv``: one line per ID, the ID then its
lineage, plus whatever columns ``args`` adds (``-n`` names, ``-r`` ranks,
``-R`` the lineage's ranks, ...). The IDs are either one ``taxid`` or a file of
IDs, one per line (``taxidfile``), never both.

``taxdb`` is the taxonomy: NCBI's ``taxdump.tar.gz`` (or any ``.tar.gz`` of a
directory holding ``nodes.dmp``, ``names.dmp``, ``delnodes.dmp`` and
``merged.dmp``), extracted in the task, or a Dir of those files. TaxonKit
runs with ``--threads`` of the task's CPUs.
"""

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

TAXONKIT_IMAGE = "community.wave.seqera.io/library/taxonkit:0.18.0--36a4cff30b2fab4e"
# Upstream's process_low; NCBI's taxdump is a few hundred MiB unpacked.
DEFAULT_RESOURCES = flyte.Resources(cpu=2, memory="4Gi", disk="10Gi")

taxonkit_lineage_cmd = shell.create(
    name="taxonkit_lineage",
    image=TAXONKIT_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "taxid": str,
        "taxidfile": File | None,
        "taxdb_archive": File | None,
        "taxdb_dir": Dir | None,
        "prefix": str,
        "args": str,
    },
    defaults={"taxid": "", "args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        shopt -s nullglob
        TAXIDFILE=({{inputs.taxidfile}}); ARCHIVE=({{inputs.taxdb_archive}})
        shopt -u nullglob
        TAXID={{inputs.taxid}}
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        if [ ${{#ARCHIVE[@]}} -gt 0 ]; then
            mkdir taxdb
            tar -xzf "${{ARCHIVE[0]}}" -C taxdb || exit 1
            NODES=$(find taxdb -name nodes.dmp | head -n 1)
            if [ -z "$NODES" ]; then echo "the taxonomy archive holds no nodes.dmp" >&2; exit 1; fi
            TAXDB=$(dirname "$NODES")
        else
            TAXDB={{inputs.taxdb_dir}}
        fi
        if [ ${{#TAXIDFILE[@]}} -gt 0 ]; then
            IDS="${{TAXIDFILE[0]}}"
        else
            printf '%s\n' "$TAXID" > taxid.txt
            IDS=taxid.txt
        fi
        taxonkit \
            lineage \
            $ARGS \
            --data-dir "$TAXDB" \
            --threads $TASK_CPUS \
            --out-file "$OUT/$P.tsv" \
            "$IDS" || exit 1
    """,
)

env = flyte.TaskEnvironment.from_task("taxonkit", taxonkit_lineage_cmd.as_task())


async def taxonkit_lineage(
    taxdb: File | Dir,
    taxid: str = "",
    taxidfile: File | None = None,
    prefix: str = "sample",
    args: str = "",
) -> File:
    """The lineage of ``taxid`` or of each ID in ``taxidfile`` (exactly one), against ``taxdb``: ``<prefix>.tsv``."""
    if bool(taxid) == (taxidfile is not None):
        raise ValueError("taxonkit_lineage needs exactly one of taxid and taxidfile")
    archive = taxdb if isinstance(taxdb, File) else None
    directory = taxdb if isinstance(taxdb, Dir) else None
    results = await taxonkit_lineage_cmd(
        taxid=taxid,
        taxidfile=taxidfile,
        taxdb_archive=archive,
        taxdb_dir=directory,
        prefix=prefix,
        args=args,
    )
    out = await results.get_file(f"{prefix}.tsv")
    if out is None:
        raise FileNotFoundError(f"taxonkit lineage wrote no {prefix}.tsv")
    return out
