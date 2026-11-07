"""sam2lca — lowest-common-ancestor taxonomic assignment of aligned reads.

Exposes :func:`sam2lca_analyze`, which assigns each read of an indexed BAM
(the ``.bai`` is linked next to it) to the lowest common ancestor of the
references it aligns to, and returns upstream's outputs:

- ``<prefix>.sam2lca.csv``: reads per taxon;
- ``<prefix>.sam2lca.json``: the same, with each taxon's lineage;
- with ``-b`` in ``args``: ``<prefix>.sam2lca.bam`` (reads tagged ``XT`` with
  their TAXID), plus per-taxon BAMs with ``-r``.

``database`` is a sam2lca database directory (as ``sam2lca update-db``
writes it) packed as a ``.tar.gz``: the contents, or a single top-level
directory holding them. It's extracted into the task's working directory,
which sam2lca must be able to write. Without one, sam2lca starts from an empty
database directory and builds the taxonomy and accession-to-taxon databases
that ``args`` select (``-t``, ``-a``), downloading them; upstream's test uses
its small ``test`` databases this way. The real NCBI databases are tens of
GiB: prepare them once and pass them as ``database``, with the task's memory
and disk raised.

As upstream, sam2lca's ``-p`` is left to its default (2 processes) unless
``args`` sets it.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

SAM2LCA_IMAGE = "quay.io/biocontainers/sam2lca:1.1.4--pyhdfd78af_0"
# Upstream's process_medium is 6 CPUs and 36 GB; the test databases need far less.
DEFAULT_RESOURCES = flyte.Resources(cpu=2, memory="8Gi", disk="20Gi")

sam2lca_analyze_cmd = shell.create(
    name="sam2lca_analyze",
    image=SAM2LCA_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"bam": File, "bai": File, "database": File | None, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        BAM=({inputs.bam})
        BAI=({inputs.bai})
        shopt -s nullglob; DATABASE=({inputs.database}); shopt -u nullglob
        P={inputs.prefix}
        ARGS={inputs.args}
        OUT={outputs.results}
        W=$(mktemp -d); cd "$W"
        export HOME="$W"
        mkdir input
        NAME=$(basename "${BAM[0]}")
        ln -s "${BAM[0]}" "input/$NAME"
        ln -s "${BAI[0]}" "input/$NAME.bai"

        mkdir sam2lca_db
        DB=sam2lca_db
        if [ ${#DATABASE[@]} -gt 0 ]; then
            tar -xzf "${DATABASE[0]}" -C sam2lca_db
            ENTRIES=(sam2lca_db/*)
            if [ ${#ENTRIES[@]} -eq 1 ] && [ -d "${ENTRIES[0]}" ]; then DB="${ENTRIES[0]}"; fi
        fi

        sam2lca -d "$DB" analyze $ARGS -o "$P" "input/$NAME" || exit 1
        shopt -s nullglob
        mv ./*.csv ./*.json ./*.bam "$OUT"/
    """,
)

env = flyte.TaskEnvironment.from_task("sam2lca", sam2lca_analyze_cmd.as_task())


@dataclass
class Sam2lcaResult:
    csv: File  # <prefix>.sam2lca.csv
    json: File  # <prefix>.sam2lca.json
    bam: list[File]  # with -b: <prefix>.sam2lca.bam, and per-taxon BAMs with -r
    results: Dir  # everything sam2lca wrote


async def sam2lca_analyze(
    bam: File, bai: File, database: File | None = None, prefix: str = "sample", args: str = ""
) -> Sam2lcaResult:
    """sam2lca's LCA assignment of ``bam``'s reads, against ``database`` (a ``.tar.gz``) or built from ``args``."""
    results = await sam2lca_analyze_cmd(bam=bam, bai=bai, database=database, prefix=prefix, args=args)
    files = [f async for f in results.walk()]

    def named(suffix: str) -> list[File]:
        return sorted((f for f in files if f.path.endswith(suffix)), key=lambda f: f.path)

    csv, json = named(".csv"), named(".json")
    if len(csv) != 1 or len(json) != 1:
        raise FileNotFoundError(f"sam2lca wrote {len(csv)} CSVs and {len(json)} JSONs, expected one of each")
    return Sam2lcaResult(csv=csv[0], json=json[0], bam=named(".bam"), results=results)
