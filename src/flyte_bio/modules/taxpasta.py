"""taxpasta — standardised taxonomic profiles across profilers.

- :func:`taxpasta_standardise` converts one profiler's output (``profiler``
  is taxpasta's name for it, e.g. ``kraken2``, ``bracken``, ``centrifuge``,
  ``ganon``, ``kmcp``, ``krakenuniq``, ``megan6``, ``metaphlan``, ``motus``)
  into taxpasta's standard table (``taxpasta standardise``):
  ``<prefix>.<format>``, with ``format`` one of ``tsv``, ``csv``, ``arrow``,
  ``parquet`` or ``biom``.
- :func:`taxpasta_merge` standardises several samples' profiles and merges
  them into one table (``taxpasta merge``): ``<prefix>.<format>``. Samples
  are named after the files or, given a ``samplesheet`` (a TSV of ``sample``
  and ``profile`` columns; the profiles by file name), after it.

Given a ``taxonomy`` Dir (NCBI-style ``names.dmp`` and ``nodes.dmp``),
taxpasta can add names and ranks (``--add-name`` etc. in ``args``).
"""

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

TAXPASTA_IMAGE = "quay.io/biocontainers/taxpasta:0.7.0--pyhdfd78af_1"
# Upstream labels both process_single.
DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory="4Gi")
FORMATS = ("tsv", "csv", "arrow", "parquet", "biom")

# Optional taxonomy Dir: an unset Dir input renders as a path that doesn't exist (or is empty).
TAXONOMY_OPTION = r"""
        TAXONOMY=()
        if [ -d "$TAX_DIR" ] && [ -n "$(ls -A "$TAX_DIR")" ]; then TAXONOMY=(--taxonomy "$TAX_DIR"); fi
""".strip()

taxpasta_standardise_cmd = shell.create(
    name="taxpasta_standardise",
    image=TAXPASTA_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"profile": File, "profiler": str, "format": str, "taxonomy": Dir | None, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        PROFILE=({{inputs.profile}})
        PROFILER={{inputs.profiler}}
        FORMAT={{inputs.format}}
        TAX_DIR={{inputs.taxonomy}}
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        {TAXONOMY_OPTION}
        taxpasta standardise \
            --profiler "$PROFILER" \
            --output "$P.$FORMAT" \
            $ARGS \
            "${{TAXONOMY[@]}}" \
            "${{PROFILE[0]}}" \
            || exit 1
        mv "$P.$FORMAT" "$OUT"/
    """,
)

taxpasta_merge_cmd = shell.create(
    name="taxpasta_merge",
    image=TAXPASTA_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "profiles": list[File],
        "profiles_order": str,
        "profiler": str,
        "format": str,
        "taxonomy": Dir | None,
        "samplesheet": File | None,
        "prefix": str,
        "args": str,
    },
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        PROFILES=({{inputs.profiles}})
        # Staged list inputs expand in name order; restore the order they were given in.
        read -ra PROFILES_NAMES <<< {{inputs.profiles_order}}; PROFILES_DIR=$(dirname "${{PROFILES[0]}}")
        PROFILES=(); for N in "${{PROFILES_NAMES[@]}}"; do PROFILES+=("$PROFILES_DIR/$N"); done
        shopt -s nullglob; SHEET=({{inputs.samplesheet}}); shopt -u nullglob
        PROFILER={{inputs.profiler}}
        FORMAT={{inputs.format}}
        TAX_DIR={{inputs.taxonomy}}
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        {TAXONOMY_OPTION}
        # Profiles are named after their files (and a samplesheet lists them by name), so link them here.
        NAMES=()
        for F in "${{PROFILES[@]}}"; do NAMES+=("$(basename "$F")"); ln -s "$F" "$(basename "$F")"; done
        SAMPLESHEET=()
        if [ ${{#SHEET[@]}} -gt 0 ]; then SAMPLESHEET=(-s "${{SHEET[0]}}"); fi
        taxpasta merge \
            --profiler "$PROFILER" \
            --output "$P.$FORMAT" \
            $ARGS \
            "${{TAXONOMY[@]}}" \
            "${{SAMPLESHEET[@]}}" \
            "${{NAMES[@]}}" \
            || exit 1
        mv "$P.$FORMAT" "$OUT"/
    """,
)

env = flyte.TaskEnvironment.from_task("taxpasta", taxpasta_standardise_cmd.as_task(), taxpasta_merge_cmd.as_task())


def check_format(format: str) -> None:
    if format not in FORMATS:
        raise ValueError(f"format must be one of {', '.join(FORMATS)}, got {format!r}")


async def taxpasta_standardise(
    profile: File,
    profiler: str,
    format: str = "tsv",
    taxonomy: Dir | None = None,
    prefix: str = "sample",
    args: str = "",
) -> File:
    """``profiler``'s ``profile`` as taxpasta's standard table: ``<prefix>.<format>``."""
    check_format(format)
    results = await taxpasta_standardise_cmd(
        profile=profile, profiler=profiler, format=format, taxonomy=taxonomy, prefix=prefix, args=args
    )
    table = await results.get_file(f"{prefix}.{format}")
    if table is None:
        raise FileNotFoundError(f"taxpasta standardise wrote no {prefix}.{format}")
    return table


async def taxpasta_merge(
    profiles: list[File],
    profiler: str,
    format: str = "tsv",
    taxonomy: Dir | None = None,
    samplesheet: File | None = None,
    prefix: str = "sample",
    args: str = "",
) -> File:
    """Several samples' ``profiler`` profiles standardised and merged: ``<prefix>.<format>``."""
    check_format(format)
    results = await taxpasta_merge_cmd(
        profiles=profiles,
        profiles_order=" ".join(f.name for f in profiles),
        profiler=profiler,
        format=format,
        taxonomy=taxonomy,
        samplesheet=samplesheet,
        prefix=prefix,
        args=args,
    )
    table = await results.get_file(f"{prefix}.{format}")
    if table is None:
        raise FileNotFoundError(f"taxpasta merge wrote no {prefix}.{format}")
    return table
