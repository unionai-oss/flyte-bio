"""krakentools — utilities for Kraken2 reports.

Exposes :func:`krakentools_kreport2krona`, which converts a Kraken2 report
into Krona's text format (``kreport2krona.py``): ``<prefix>.txt``, a count
then the lineage's names per line, for ``ktImportText``. Options such as
``--no-intermediate-ranks`` go in ``args``.

:func:`krakentools_combinekreports` combines several Kraken2 reports into
one (``combine_kreports.py``): ``<prefix>.txt``. Options such as
``--display-headers`` or ``--sample-names`` go in ``args``.
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

krakentools_combinekreports_cmd = shell.create(
    name="krakentools_combinekreports",
    image=KRAKENTOOLS_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"kreports": list[File], "kreports_order": str, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        REPORTS=({inputs.kreports})
        # Staged list inputs expand in name order; restore the order they were given in.
        read -ra REPORTS_NAMES <<< {inputs.kreports_order}; REPORTS_DIR=$(dirname "${REPORTS[0]}")
        REPORTS=(); for N in "${REPORTS_NAMES[@]}"; do REPORTS+=("$REPORTS_DIR/$N"); done
        P={inputs.prefix}
        ARGS={inputs.args}
        OUT={outputs.results}
        W=$(mktemp -d); cd "$W"
        NAMES=()
        for F in "${REPORTS[@]}"; do NAMES+=("$(basename "$F")"); ln -s "$F" "$(basename "$F")"; done
        combine_kreports.py \
            -r "${NAMES[@]}" \
            -o "$P.txt" \
            $ARGS \
            || exit 1
        mv "$P.txt" "$OUT"/
    """,
)

env = flyte.TaskEnvironment.from_task(
    "krakentools", krakentools_kreport2krona_cmd.as_task(), krakentools_combinekreports_cmd.as_task()
)


async def krakentools_kreport2krona(kreport: File, prefix: str = "sample", args: str = "") -> File:
    """The Kraken2 report ``kreport`` in Krona's text format: ``<prefix>.txt``."""
    results = await krakentools_kreport2krona_cmd(kreport=kreport, prefix=prefix, args=args)
    txt = await results.get_file(f"{prefix}.txt")
    if txt is None:
        raise FileNotFoundError(f"kreport2krona.py wrote no {prefix}.txt")
    return txt


async def krakentools_combinekreports(kreports: File | list[File], prefix: str = "sample", args: str = "") -> File:
    """Several Kraken2 reports combined into one: ``<prefix>.txt``."""
    reports = kreports if isinstance(kreports, list) else [kreports]
    results = await krakentools_combinekreports_cmd(
        kreports=reports, kreports_order=" ".join(f.name for f in reports), prefix=prefix, args=args
    )
    combined = await results.get_file(f"{prefix}.txt")
    if combined is None:
        raise FileNotFoundError(f"combine_kreports.py wrote no {prefix}.txt")
    return combined
