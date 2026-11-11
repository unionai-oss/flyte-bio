"""shinyngs — validate, plot and package expression studies with shinyngs.

Exposes, each running one of shinyngs' command-line scripts:

- :func:`shinyngs_validatefomcomponents` (``validate_fom_components.R``) —
  check a study's sample sheet, assay matrices, feature metadata and
  contrasts (CSV, TSV or YAML) against each other, writing tidied copies:
  ``<sheet>.sample_metadata.tsv``, ``<matrix>.assay.tsv`` per matrix,
  ``<features>.feature_metadata.tsv`` and ``<contrasts>.contrasts_file.tsv``.
- :func:`shinyngs_staticexploratory` (``exploratory_plots.R``) — exploratory
  plots of the assays, coloured by a sample sheet variable: ``boxplot``,
  ``density``, ``pca2d``, ``pca3d``, ``mad_correlation`` and
  ``sample_dendrogram`` PNGs, and with ``--write_html`` in ``args``, HTML
  versions of all but the dendrogram.
- :func:`shinyngs_staticdifferential` (``differential_plots.R``) — a volcano
  plot of a differential results table: ``volcano.png`` (and
  ``volcano.html`` with ``--write_html``).
- :func:`shinyngs_app` (``make_app_from_files.R``) — a Shiny app of the
  study: its ``data.rds`` and ``app.R``, optionally with gene sets and
  enrichment tables (named in ``args``, e.g. ``--enrichment_gene_sets`` and
  ``--enrichment_filename_template``; they're staged under their own
  names). This port builds the app locally only: upstream can also deploy it
  to shinyapps.io (``--deploy_app``, with shinyapps.io credentials as
  secrets), which isn't wired here.

``args`` are upstream's ``ext.args``, parsed as the shell would (so quoted
values keep their spaces).
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

SHINYNGS_IMAGE = "community.wave.seqera.io/library/r-shinyngs:3.2.1--d43071e62bc500d3"
# Upstream: process_single.
DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory="8Gi")

# Stage each file of the named arrays into the working directory under its own name, and parse args
# as the shell would.
STAGE = r"""
        W=$(mktemp -d); cd "$W"
        stage() { for f in "$@"; do ln -s "$f" "$(basename "$f")"; done; }
        names() { local IFS=,; local out=(); for f in "$@"; do out+=("$(basename "$f")"); done; echo "${out[*]}"; }
        eval "ARGS=($ARGS_STRING)"
"""

shinyngs_validatefomcomponents_cmd = shell.create(
    name="shinyngs_validatefomcomponents",
    image=SHINYNGS_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "sample": File,
        "assay_files": list[File],
        "feature_meta": File | None,
        "contrasts": File,
        "prefix": str,
        "args": str,
    },
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        SAMPLE=({inputs.sample})
        ASSAYS=({inputs.assay_files})
        shopt -s nullglob; FEATURES=({inputs.feature_meta}); shopt -u nullglob
        CONTRASTS=({inputs.contrasts})
        P={inputs.prefix}
        ARGS_STRING={inputs.args}
        STAGE
        stage "${SAMPLE[@]}" "${ASSAYS[@]}" "${FEATURES[@]}" "${CONTRASTS[@]}"
        FEATURE=()
        if [ ${#FEATURES[@]} -gt 0 ]; then FEATURE=(--feature_metadata "$(basename "${FEATURES[0]}")"); fi
        validate_fom_components.R \
            --sample_metadata "$(basename "${SAMPLE[0]}")" \
            "${FEATURE[@]}" \
            --assay_files "$(names "${ASSAYS[@]}")" \
            --contrasts_file "$(basename "${CONTRASTS[0]}")" \
            --output_directory "$P" \
            "${ARGS[@]}" || exit 1
        mv "$P"/* {outputs.results}/
    """.replace("STAGE", STAGE.strip()),
)

shinyngs_staticexploratory_cmd = shell.create(
    name="shinyngs_staticexploratory",
    image=SHINYNGS_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "sample": File,
        "feature_meta": File,
        "assay_files": list[File],
        "variable": str,
        "prefix": str,
        "args": str,
    },
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        SAMPLE=({inputs.sample})
        FEATURES=({inputs.feature_meta})
        ASSAYS=({inputs.assay_files})
        P={inputs.prefix}
        ARGS_STRING={inputs.args}
        STAGE
        stage "${SAMPLE[@]}" "${FEATURES[@]}" "${ASSAYS[@]}"
        exploratory_plots.R \
            --sample_metadata "$(basename "${SAMPLE[0]}")" \
            --feature_metadata "$(basename "${FEATURES[0]}")" \
            --assay_files "$(names "${ASSAYS[@]}")" \
            --contrast_variable {inputs.variable} \
            --outdir "$P" \
            "${ARGS[@]}" || exit 1
        shopt -s nullglob
        mv "$P"/png/*.png "$P"/html/*.html {outputs.results}/
    """.replace("STAGE", STAGE.strip()),
)

shinyngs_staticdifferential_cmd = shell.create(
    name="shinyngs_staticdifferential",
    image=SHINYNGS_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"differential_result": File, "feature_meta": File, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        DIFFERENTIAL=({inputs.differential_result})
        FEATURES=({inputs.feature_meta})
        P={inputs.prefix}
        ARGS_STRING={inputs.args}
        STAGE
        stage "${DIFFERENTIAL[@]}" "${FEATURES[@]}"
        differential_plots.R \
            --differential_file "$(basename "${DIFFERENTIAL[0]}")" \
            --feature_metadata "$(basename "${FEATURES[0]}")" \
            --outdir "$P" \
            "${ARGS[@]}" || exit 1
        mv "$P"/png/volcano.png {outputs.results}/
        if [ -e "$P"/html/volcano.html ]; then mv "$P"/html/volcano.html {outputs.results}/; fi
    """.replace("STAGE", STAGE.strip()),
)

shinyngs_app_cmd = shell.create(
    name="shinyngs_app",
    image=SHINYNGS_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "sample": File,
        "feature_meta": File,
        "assay_files": list[File],
        "contrasts": File,
        "differential_results": list[File],
        "contrast_stats_assay": int,
        "gene_sets": list[File],
        "enrichment_results": list[File],
        "prefix": str,
        "args": str,
    },
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        SAMPLE=({inputs.sample})
        FEATURES=({inputs.feature_meta})
        ASSAYS=({inputs.assay_files})
        CONTRASTS=({inputs.contrasts})
        DIFFERENTIALS=({inputs.differential_results})
        shopt -s nullglob; GENE_SETS=({inputs.gene_sets}); ENRICHMENT=({inputs.enrichment_results}); shopt -u nullglob
        P={inputs.prefix}
        ARGS_STRING={inputs.args}
        STAGE
        stage "${SAMPLE[@]}" "${FEATURES[@]}" "${ASSAYS[@]}" "${CONTRASTS[@]}" "${DIFFERENTIALS[@]}" \
            "${GENE_SETS[@]}" "${ENRICHMENT[@]}"
        make_app_from_files.R \
            --sample_metadata "$(basename "${SAMPLE[0]}")" \
            --feature_metadata "$(basename "${FEATURES[0]}")" \
            --assay_files "$(names "${ASSAYS[@]}")" \
            --contrast_file "$(basename "${CONTRASTS[0]}")" \
            --contrast_stats_assay {inputs.contrast_stats_assay} \
            --differential_results "$(names "${DIFFERENTIALS[@]}")" \
            --output_dir "$P" \
            "${ARGS[@]}" || exit 1
        mv "$P"/data.rds "$P"/app.R {outputs.results}/
    """.replace("STAGE", STAGE.strip()),
)

env = flyte.TaskEnvironment.from_task(
    "shinyngs",
    shinyngs_validatefomcomponents_cmd.as_task(),
    shinyngs_staticexploratory_cmd.as_task(),
    shinyngs_staticdifferential_cmd.as_task(),
    shinyngs_app_cmd.as_task(),
)


def by_suffix(files: list[File], suffix: str) -> list[File]:
    return sorted((f for f in files if f.name.endswith(suffix)), key=lambda f: f.name)


def only(files: list[File], suffix: str, tool: str) -> File:
    found = by_suffix(files, suffix)
    if len(found) != 1:
        raise FileNotFoundError(f"{tool} wrote {len(found)} *{suffix} files, expected one")
    return found[0]


@dataclass
class ValidateFomComponentsResult:
    sample_meta: File  # <sheet>.sample_metadata.tsv
    feature_meta: File | None  # <features>.feature_metadata.tsv
    assays: list[File]  # <matrix>.assay.tsv, by name
    contrasts: File  # <contrasts>.contrasts_file.tsv


async def shinyngs_validatefomcomponents(
    sample: File,
    assay_files: list[File],
    contrasts: File,
    prefix: str,
    feature_meta: File | None = None,
    args: str = "",
) -> ValidateFomComponentsResult:
    """Validated, tidied copies of a study's sample sheet, assays, features and contrasts."""
    results = await shinyngs_validatefomcomponents_cmd(
        sample=sample, assay_files=assay_files, feature_meta=feature_meta, contrasts=contrasts, prefix=prefix, args=args
    )
    files = await results.list_files()
    tool = "validate_fom_components.R"
    features = by_suffix(files, ".feature_metadata.tsv")
    return ValidateFomComponentsResult(
        sample_meta=only(files, ".sample_metadata.tsv", tool),
        feature_meta=features[0] if features else None,
        assays=by_suffix(files, ".assay.tsv"),
        contrasts=only(files, ".contrasts_file.tsv", tool),
    )


@dataclass
class StaticExploratoryResult:
    boxplots_png: File  # boxplot.png
    densities_png: File  # density.png
    pca2d_png: File  # pca2d.png
    pca3d_png: File  # pca3d.png
    mad_png: File | None  # mad_correlation.png
    dendro: File  # sample_dendrogram.png
    boxplots_html: File | None  # boxplot.html, with --write_html
    densities_html: File | None  # density.html
    pca2d_html: File | None  # pca2d.html
    pca3d_html: File | None  # pca3d.html
    mad_html: File | None  # mad_correlation.html


async def shinyngs_staticexploratory(
    sample: File, feature_meta: File, assay_files: list[File], variable: str, prefix: str, args: str = ""
) -> StaticExploratoryResult:
    """Exploratory plots of ``assay_files``, coloured by the sample sheet's ``variable``."""
    results = await shinyngs_staticexploratory_cmd(
        sample=sample, feature_meta=feature_meta, assay_files=assay_files, variable=variable, prefix=prefix, args=args
    )
    by_name = {f.name: f for f in await results.list_files()}

    def required(name: str) -> File:
        if name not in by_name:
            raise FileNotFoundError(f"exploratory_plots.R wrote no {name}")
        return by_name[name]

    return StaticExploratoryResult(
        boxplots_png=required("boxplot.png"),
        densities_png=required("density.png"),
        pca2d_png=required("pca2d.png"),
        pca3d_png=required("pca3d.png"),
        mad_png=by_name.get("mad_correlation.png"),
        dendro=required("sample_dendrogram.png"),
        boxplots_html=by_name.get("boxplot.html"),
        densities_html=by_name.get("density.html"),
        pca2d_html=by_name.get("pca2d.html"),
        pca3d_html=by_name.get("pca3d.html"),
        mad_html=by_name.get("mad_correlation.html"),
    )


@dataclass
class StaticDifferentialResult:
    volcanos_png: File  # volcano.png
    volcanos_html: File | None  # volcano.html, with --write_html


async def shinyngs_staticdifferential(
    differential_result: File, feature_meta: File, prefix: str, args: str = ""
) -> StaticDifferentialResult:
    """A volcano plot of ``differential_result`` (the contrast's levels and columns come from ``args``)."""
    results = await shinyngs_staticdifferential_cmd(
        differential_result=differential_result, feature_meta=feature_meta, prefix=prefix, args=args
    )
    png = await results.get_file("volcano.png")
    if png is None:
        raise FileNotFoundError("differential_plots.R wrote no volcano.png")
    return StaticDifferentialResult(volcanos_png=png, volcanos_html=await results.get_file("volcano.html"))


@dataclass
class ShinyngsAppResult:
    data: File  # data.rds
    app: File  # app.R


async def shinyngs_app(
    sample: File,
    feature_meta: File,
    assay_files: list[File],
    contrasts: File,
    differential_results: list[File],
    contrast_stats_assay: int,
    prefix: str,
    gene_sets: list[File] | None = None,
    enrichment_results: list[File] | None = None,
    args: str = "",
) -> ShinyngsAppResult:
    """A shinyngs app (``data.rds`` and ``app.R``) of a study and its contrasts' differential results."""
    results = await shinyngs_app_cmd(
        sample=sample,
        feature_meta=feature_meta,
        assay_files=assay_files,
        contrasts=contrasts,
        differential_results=differential_results,
        contrast_stats_assay=contrast_stats_assay,
        gene_sets=gene_sets or [],
        enrichment_results=enrichment_results or [],
        prefix=prefix,
        args=args,
    )
    data = await results.get_file("data.rds")
    app = await results.get_file("app.R")
    if data is None or app is None:
        raise FileNotFoundError("make_app_from_files.R wrote no data.rds or app.R")
    return ShinyngsAppResult(data=data, app=app)
