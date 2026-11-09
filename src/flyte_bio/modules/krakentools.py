"""krakentools — utilities for Kraken2 reports.

Exposes :func:`krakentools_kreport2krona`, which converts a Kraken2 report
into Krona's text format (``kreport2krona.py``): ``<prefix>.txt``, a count
then the lineage's names per line, for ``ktImportText``. Options such as
``--no-intermediate-ranks`` go in ``args``.
"""

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

KRAKENTOOLS_IMAGE = "quay.io/biocontainers/krakentools:1.2.1--pyh7e72e81_0"
# Upstream: process_single.
DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory="2Gi")

krakentools_kreport2krona_cmd = shell.create(
    name="krakentools_kreport2krona",
    image=KRAKENTOOLS_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"kreport": File, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        REPORT=({inputs.kreport})
        P={inputs.prefix}
        ARGS={inputs.args}
        OUT={outputs.results}
        W=$(mktemp -d); cd "$W"
        kreport2krona.py \
            -r "${REPORT[0]}" \
            -o "$P.txt" \
            $ARGS \
            || exit 1
        mv "$P.txt" "$OUT"/
    """,
)

env = flyte.TaskEnvironment.from_task("krakentools", krakentools_kreport2krona_cmd.as_task())


async def krakentools_kreport2krona(kreport: File, prefix: str = "sample", args: str = "") -> File:
    """The Kraken2 report ``kreport`` in Krona's text format: ``<prefix>.txt``."""
    results = await krakentools_kreport2krona_cmd(kreport=kreport, prefix=prefix, args=args)
    txt = await results.get_file(f"{prefix}.txt")
    if txt is None:
        raise FileNotFoundError(f"kreport2krona.py wrote no {prefix}.txt")
    return txt
