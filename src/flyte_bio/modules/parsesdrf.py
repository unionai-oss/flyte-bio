"""parsesdrf — convert an SDRF proteomics sample sheet into a tool's input files.

Exposes :func:`parsesdrf_convert`, which runs sdrf-pipelines'
``parse_sdrf convert-<format>`` on an SDRF and returns the format's files,
in upstream's order:

- ``mhcquant``: ``<prefix>_samplesheet.tsv``, ``<prefix>_presets.tsv``;
- ``openms``: ``<prefix>_samplesheet.tsv`` and, when written,
  ``<prefix>_experimental_design.tsv``;
- ``maxquant`` (needs ``fasta`` and ``raw_folder``, the folder written into
  the parameter file): ``<prefix>.xml``, ``<prefix>_design.txt``;
- ``msstats``: ``<prefix>.csv``;
- ``normalyzerde``: ``<prefix>_design.csv``, ``<prefix>_comparisons.csv``;
- ``diann``: ``<prefix>.cfg``, ``<prefix>_design.tsv``.
"""

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

PARSESDRF_IMAGE = "quay.io/biocontainers/sdrf-pipelines:0.1.4--pyhdfd78af_0"
# Upstream labels it process_single.
DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory="2Gi")

# Each format's output names, with "{p}" for the prefix; the optional ones may be missing.
FORMAT_OUTPUTS = {
    "mhcquant": ["{p}_samplesheet.tsv", "{p}_presets.tsv"],
    "openms": ["{p}_samplesheet.tsv"],
    "maxquant": ["{p}.xml", "{p}_design.txt"],
    "msstats": ["{p}.csv"],
    "normalyzerde": ["{p}_design.csv", "{p}_comparisons.csv"],
    "diann": ["{p}.cfg", "{p}_design.tsv"],
}
OPTIONAL_OUTPUTS = {"openms": ["{p}_experimental_design.tsv"]}

parsesdrf_convert_cmd = shell.create(
    name="parsesdrf_convert",
    image=PARSESDRF_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"sdrf": File, "fasta": File | None, "raw_folder": str, "format": str, "prefix": str, "args": str},
    defaults={"raw_folder": "", "args": ""},
    outputs={"results": Dir},
    script=r"""
        SDRF=({inputs.sdrf})
        shopt -s nullglob; FA=({inputs.fasta}); shopt -u nullglob
        RAW_FOLDER={inputs.raw_folder}
        FORMAT={inputs.format}
        P={inputs.prefix}
        ARGS={inputs.args}
        OUT={outputs.results}
        W=$(mktemp -d); cd "$W"
        export HOME="$W"
        ln -s "${SDRF[0]}" .
        SDRF_NAME=$(basename "${SDRF[0]}")
        OUTPUT_ARGS=()
        case "$FORMAT" in
            mhcquant) OUTPUT_ARGS=(-os "${P}_samplesheet.tsv" -op "${P}_presets.tsv") ;;
            openms) ;;
            maxquant)
                ln -s "${FA[0]}" .
                OUTPUT_ARGS=(-f "$(basename "${FA[0]}")" -r "$RAW_FOLDER" -o1 "$P.xml" -o2 "${P}_design.txt") ;;
            msstats) OUTPUT_ARGS=(-o "$P.csv") ;;
            normalyzerde) OUTPUT_ARGS=(-o "${P}_design.csv" -oc "${P}_comparisons.csv") ;;
            diann) ;;
            *) echo "unsupported format $FORMAT" >&2; exit 1 ;;
        esac
        parse_sdrf convert-"$FORMAT" -s "$SDRF_NAME" "${OUTPUT_ARGS[@]}" $ARGS || exit 1
        if [ "$FORMAT" = openms ]; then
            mv openms.tsv "${P}_samplesheet.tsv"
            if [ -f experimental_design.tsv ]; then mv experimental_design.tsv "${P}_experimental_design.tsv"; fi
        elif [ "$FORMAT" = diann ]; then
            mv diann_config.cfg "$P.cfg" && mv diann_design.tsv "${P}_design.tsv"
        fi
        for f in "$P"*; do
            if [ -f "$f" ] && [ ! -L "$f" ]; then mv "$f" "$OUT"/; fi
        done
    """,
)

env = flyte.TaskEnvironment.from_task("parsesdrf", parsesdrf_convert_cmd.as_task())


async def parsesdrf_convert(
    sdrf: File,
    format: str,
    fasta: File | None = None,
    raw_folder: str = "",
    prefix: str = "sample",
    args: str = "",
) -> list[File]:
    """``sdrf`` converted to ``format``'s input files, in upstream's order (see the module docstring)."""
    if format not in FORMAT_OUTPUTS:
        raise ValueError(f"unsupported format {format!r}: expected one of {', '.join(FORMAT_OUTPUTS)}")
    if format == "maxquant" and (fasta is None or not raw_folder):
        raise ValueError("format 'maxquant' needs both fasta and raw_folder")
    results = await parsesdrf_convert_cmd(
        sdrf=sdrf, fasta=fasta, raw_folder=raw_folder, format=format, prefix=prefix, args=args
    )
    files = []
    for pattern in FORMAT_OUTPUTS[format]:
        name = pattern.format(p=prefix)
        found = await results.get_file(name)
        if found is None:
            raise FileNotFoundError(f"parse_sdrf wrote no {name}")
        files.append(found)
    for pattern in OPTIONAL_OUTPUTS.get(format, []):
        found = await results.get_file(pattern.format(p=prefix))
        if found is not None:
            files.append(found)
    return files
