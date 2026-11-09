"""samshee — validation and reformatting of Illumina sample sheets.

Exposes :func:`samshee`, which validates an Illumina sample sheet (v1 or v2)
with samshee and writes the sheet as samshee formats it,
``<sample sheet name without extension>_formatted.csv``. samshee fails on a
sheet that doesn't validate.

Validation uses samshee's default checks plus, when given:

- ``schema``: one JSON schema, inline (e.g. ``{"required": ["Data"]}``) or a
  reference (``{"$ref": "urn:samshee:illuminav2/v1"}``);
- ``schema_file``: a JSON schema file, passed as ``{"$ref": "file:<path>"}``,
  as upstream.

``args`` takes samshee's other options, e.g. ``--output-format sectioned``
to write a v1 sheet.
"""

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

SAMSHEE_IMAGE = "community.wave.seqera.io/library/samshee_python:38088c103ef2751a"
# Upstream's process_single.
DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory="2Gi")

samshee_cmd = shell.create(
    name="samshee",
    image=SAMSHEE_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"samplesheet": File, "schema_file": File | None, "schema": str, "args": str},
    defaults={"schema": "", "args": ""},
    outputs={"results": Dir},
    script=r"""
        set -e
        SHEET=({inputs.samplesheet})
        shopt -s nullglob; SCHEMA_FILE=({inputs.schema_file}); shopt -u nullglob
        SCHEMA={inputs.schema}
        ARGS={inputs.args}
        SCHEMAS=()
        if [ -n "$SCHEMA" ]; then SCHEMAS+=(--schema "$SCHEMA"); fi
        if [ ${#SCHEMA_FILE[@]} -gt 0 ]; then SCHEMAS+=(--schema "{\"\$ref\": \"file:${SCHEMA_FILE[0]}\"}"); fi
        NAME=$(basename "${SHEET[0]}")
        python -m samshee "${SHEET[0]}" $ARGS "${SCHEMAS[@]}" > {outputs.results}/"${NAME%.*}_formatted.csv"
    """,
)

env = flyte.TaskEnvironment.from_task("samshee", samshee_cmd.as_task())


async def samshee(samplesheet: File, schema_file: File | None = None, schema: str = "", args: str = "") -> File:
    """``samplesheet`` validated and formatted by samshee, named ``<sheet name>_formatted.csv``."""
    results = await samshee_cmd(samplesheet=samplesheet, schema_file=schema_file, schema=schema, args=args)
    files = await results.list_files()
    if len(files) != 1:
        raise FileNotFoundError(f"samshee wrote {len(files)} files, expected one formatted sheet")
    return files[0]
