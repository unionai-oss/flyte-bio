"""pridepy — download files from the PRIDE proteomics archive.

Exposes:

- :func:`pridepy_downloadfile` — download one file of a PRIDE project by
  name, keeping its name.
- :func:`pridepy_fetchsdrf` — find a project's SDRF sample sheet
  (``sdrf.tsv`` or ``*.sdrf.tsv``) in its file metadata and download it as
  ``<prefix>.sdrf.tsv``.

Both download from PRIDE when they run. As upstream, an empty download (which
pridepy doesn't report as an error) fails the task.
"""

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

PRIDEPY_IMAGE = "quay.io/biocontainers/pridepy:0.0.15--pyhdfd78af_0"
# Upstream labels both process_single.
DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory="2Gi")

pridepy_downloadfile_cmd = shell.create(
    name="pridepy_downloadfile",
    image=PRIDEPY_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"file_name": str, "pride_accession": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        NAME={inputs.file_name}
        ACCESSION={inputs.pride_accession}
        ARGS={inputs.args}
        OUT={outputs.results}
        W=$(mktemp -d); cd "$W"
        export HOME="$W"
        pridepy download-file-by-name -a "$ACCESSION" -f "$NAME" -o . $ARGS || exit 1
        [ -s "$NAME" ] || { echo "ERROR: Downloaded file $NAME is empty" >&2; exit 1; }
        mv "$NAME" "$OUT"/
    """,
)

pridepy_fetchsdrf_cmd = shell.create(
    name="pridepy_fetchsdrf",
    image=PRIDEPY_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"pride_id": str, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        PRIDE_ID={inputs.pride_id}
        P={inputs.prefix}
        ARGS={inputs.args}
        OUT={outputs.results}
        W=$(mktemp -d); cd "$W"
        export HOME="$W"
        pridepy stream-files-metadata -a "$PRIDE_ID" -o files_metadata.json $ARGS || exit 1
        SDRF_NAME=$(awk -F'"' '/"fileName" *: *"(sdrf\.tsv|[^"]+\.sdrf\.tsv)"/ {print $4; exit}' files_metadata.json)
        [ -n "$SDRF_NAME" ] || { echo "ERROR: No SDRF file found for $PRIDE_ID" >&2; exit 1; }
        pridepy download-file-by-name -a "$PRIDE_ID" -f "$SDRF_NAME" -o . $ARGS || exit 1
        mv "$SDRF_NAME" "$OUT/$P.sdrf.tsv"
    """,
)

env = flyte.TaskEnvironment.from_task("pridepy", pridepy_downloadfile_cmd.as_task(), pridepy_fetchsdrf_cmd.as_task())


async def pridepy_pick(results: Dir, name: str) -> File:
    found = await results.get_file(name)
    if found is None:
        raise FileNotFoundError(f"pridepy wrote no {name}")
    return found


async def pridepy_downloadfile(file_name: str, pride_accession: str, args: str = "") -> File:
    """The file ``file_name`` of PRIDE project ``pride_accession``."""
    results = await pridepy_downloadfile_cmd(file_name=file_name, pride_accession=pride_accession, args=args)
    return await pridepy_pick(results, file_name)


async def pridepy_fetchsdrf(pride_id: str, prefix: str = "sample", args: str = "") -> File:
    """PRIDE project ``pride_id``'s SDRF sample sheet, as ``<prefix>.sdrf.tsv``."""
    results = await pridepy_fetchsdrf_cmd(pride_id=pride_id, prefix=prefix, args=args)
    return await pridepy_pick(results, f"{prefix}.sdrf.tsv")
