"""thermorawfileparser — convert Thermo RAW mass-spectrometry files to open formats.

Exposes :func:`thermorawfileparser`, which converts one Thermo ``.raw`` file
with ThermoRawFileParser and returns the spectra as ``<prefix>.<suffix>``. As
upstream, the suffix follows the output format chosen in ``args``
(``--format``/``-f``):

- ``0``: ``mgf``;
- ``1`` or ``2`` (indexed and plain mzML), or none: ``mzML``;
- ``3``: ``parquet``;

with ``.gz`` added when ``args`` has ``--gzip``. The spectra record the input
file's path (at the top) and a checksum (at the end).
"""

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

THERMORAWFILEPARSER_IMAGE = "quay.io/biocontainers/thermorawfileparser:1.4.5--h05cac1d_1"
# Upstream's process_low asks for 2 CPUs and 12 GB; a RAW file of tens of MB
# converts in well under that.
DEFAULT_RESOURCES = flyte.Resources(cpu=2, memory="4Gi")

thermorawfileparser_cmd = shell.create(
    name="thermorawfileparser",
    image=THERMORAWFILEPARSER_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"raw": File, "output_name": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        RAW=({inputs.raw})
        NAME={inputs.output_name}
        ARGS={inputs.args}
        OUT={outputs.results}
        W=$(mktemp -d); cd "$W"
        export HOME="$W"
        ThermoRawFileParser.sh \
            -i "${RAW[0]}" \
            -b "$OUT/$NAME" \
            $ARGS || exit 1
    """,
)

env = flyte.TaskEnvironment.from_task("thermorawfileparser", thermorawfileparser_cmd.as_task())


def spectra_suffix(args: str) -> str:
    """The spectra file's suffix for ``args``, as upstream picks it."""
    if "--format 0" in args or "-f 0" in args:
        suffix = "mgf"
    elif "--format 3" in args or "-f 3" in args:
        suffix = "parquet"
    else:
        suffix = "mzML"
    return f"{suffix}.gz" if "--gzip" in args else suffix


async def thermorawfileparser(raw: File, prefix: str = "sample", args: str = "") -> File:
    """``raw`` converted by ThermoRawFileParser: ``<prefix>.mzML`` (or the format ``args`` picks)."""
    name = f"{prefix}.{spectra_suffix(args)}"
    results = await thermorawfileparser_cmd(raw=raw, output_name=name, args=args)
    found = await results.get_file(name)
    if found is None:
        raise FileNotFoundError(f"ThermoRawFileParser wrote no {name}")
    return found
