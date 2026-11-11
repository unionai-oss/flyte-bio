"""sylphtax — taxonomic profiles from sylph results.

Exposes :func:`sylphtax_taxprof`, which runs ``sylph-tax taxprof`` on a
``sylph profile`` table with taxonomy metadata files and writes upstream's
MetaPhlAn-style profile ``<prefix>.sylphmpa``. As upstream, sylph-tax's config
lives in ``/tmp/config.json``, so no taxonomy is downloaded unless ``-t``
names a built-in one.

:func:`sylphtax_merge` merges several samples' ``.sylphmpa`` profiles into one
table (``sylph-tax merge``), ``<prefix>.tsv``, of the column ``data_type``
(e.g. ``relative_abundance``).
"""

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

SYLPHTAX_IMAGE = "quay.io/biocontainers/sylph-tax:1.9.0--pyhdfd78af_0"
# Upstream labels it process_medium; it is single-threaded.
DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory="4Gi")

sylphtax_taxprof_cmd = shell.create(
    name="sylphtax_taxprof",
    image=SYLPHTAX_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"sylph_results": File, "taxonomy": list[File], "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        RESULTS=({inputs.sylph_results}); TAX=({inputs.taxonomy})
        P={inputs.prefix}
        ARGS={inputs.args}
        OUT={outputs.results}
        W=$(mktemp -d); cd "$W"
        ln -s "${RESULTS[0]}" "$(basename "${RESULTS[0]}")"
        TAXONOMY=()
        for F in "${TAX[@]}"; do TAXONOMY+=("$(basename "$F")"); ln -s "$F" "$(basename "$F")"; done
        export SYLPH_TAXONOMY_CONFIG="/tmp/config.json"
        sylph-tax \
            taxprof \
            "$(basename "${RESULTS[0]}")" \
            $ARGS \
            -t "${TAXONOMY[@]}" \
            || exit 1
        OUTPUTS=(./*.sylphmpa)
        if [ ${#OUTPUTS[@]} -ne 1 ] || [ ! -e "${OUTPUTS[0]}" ]; then
            echo "sylph-tax wrote ${OUTPUTS[*]}" >&2; exit 1
        fi
        mv "${OUTPUTS[0]}" "$OUT"/"$P".sylphmpa
    """,
)

sylphtax_merge_cmd = shell.create(
    name="sylphtax_merge",
    image=SYLPHTAX_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"reports": list[File], "reports_order": str, "data_type": str, "prefix": str},
    outputs={"results": Dir},
    script=r"""
        REPORTS=({inputs.reports})
        # Staged list inputs expand in name order; restore the order they were given in.
        read -ra REPORTS_NAMES <<< {inputs.reports_order}; REPORTS_DIR=$(dirname "${REPORTS[0]}")
        REPORTS=(); for N in "${REPORTS_NAMES[@]}"; do REPORTS+=("$REPORTS_DIR/$N"); done
        DATA_TYPE={inputs.data_type}
        P={inputs.prefix}
        OUT={outputs.results}
        W=$(mktemp -d); cd "$W"
        NAMES=()
        for F in "${REPORTS[@]}"; do NAMES+=("$(basename "$F")"); ln -s "$F" "$(basename "$F")"; done
        export SYLPH_TAXONOMY_CONFIG="/tmp/config.json"
        sylph-tax \
            merge \
            "${NAMES[@]}" \
            --column "$DATA_TYPE" \
            --output "$P.tsv" \
            || exit 1
        mv "$P.tsv" "$OUT"/
    """,
)

env = flyte.TaskEnvironment.from_task("sylphtax", sylphtax_taxprof_cmd.as_task(), sylphtax_merge_cmd.as_task())


async def sylphtax_taxprof(
    sylph_results: File, taxonomy: File | list[File], prefix: str = "sample", args: str = ""
) -> File:
    """The taxonomic profile of a ``sylph profile`` table, from ``taxonomy`` metadata: ``<prefix>.sylphmpa``."""
    taxonomies = taxonomy if isinstance(taxonomy, list) else [taxonomy]
    results = await sylphtax_taxprof_cmd(sylph_results=sylph_results, taxonomy=taxonomies, prefix=prefix, args=args)
    profile = await results.get_file(f"{prefix}.sylphmpa")
    if profile is None:
        raise FileNotFoundError(f"sylph-tax taxprof wrote no {prefix}.sylphmpa")
    return profile


async def sylphtax_merge(reports: File | list[File], data_type: str, prefix: str = "sample") -> File:
    """Several samples' ``.sylphmpa`` profiles merged into one table of ``data_type``: ``<prefix>.tsv``."""
    files = reports if isinstance(reports, list) else [reports]
    results = await sylphtax_merge_cmd(
        reports=files, reports_order=" ".join(f.name for f in files), data_type=data_type, prefix=prefix
    )
    merged = await results.get_file(f"{prefix}.tsv")
    if merged is None:
        raise FileNotFoundError(f"sylph-tax merge wrote no {prefix}.tsv")
    return merged
