"""custom — small table utilities for differential analysis.

Exposes:

- :func:`custom_filterdifferentialtable` — keep the rows of a differential
  results table (``.csv``, ``.tsv`` or ``.txt``) passing a fold-change and a
  statistic threshold: ``<prefix>_filtered.tsv``, and its up- and
  down-regulated rows, ``<prefix>_filtered_up.tsv`` and
  ``<prefix>_filtered_down.tsv``. The fold-change threshold is unlogged
  (its log2 is compared with the absolute log2 fold change).
- :func:`custom_matrixfilter` — keep an abundance matrix's features passing
  abundance and missing-value tests (``args`` as upstream's: e.g.
  ``--minimum_abundance 10 --minimum_samples 1``): ``<prefix>.filtered.tsv``,
  ``<prefix>.tests.tsv`` (each feature's test results),
  ``<prefix>.thresholds.tsv`` and ``R_sessionInfo.log``.
- :func:`custom_tabulartogseagct` — a tabular matrix as a GSEA ``.gct``.
- :func:`custom_tabulartogseacls` — a sample sheet's ``variable`` column as
  GSEA ``.cls`` phenotype labels.
- :func:`custom_tabulartogseachip` — a table's ID and symbol columns as a
  GSEA ``.chip`` mapping.

The filtering scripts are upstream's (the inline Python and the R template),
unchanged but for taking their templated values (inputs, ``args``, prefix) as
command-line arguments (and, so the shell template doesn't take them for
placeholders, ``{pd.__version__!s}`` for ``{pd.__version__}`` in the Python's
f-string); ``versions.yml`` names upstream's process, as upstream writes it.
The GSEA conversions are upstream's shell commands.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

PANDAS_IMAGE = "community.wave.seqera.io/library/pandas_python:3116f093aa26dc56"
R_BASE_IMAGE = "community.wave.seqera.io/library/r-base:4.2.1--d903e969fe33d78c"
GAWK_IMAGE = "community.wave.seqera.io/library/gawk:5.1.0--e2a485f195d5814e"
COREUTILS_IMAGE = "community.wave.seqera.io/library/coreutils:9.5--e94c8172b8bba3cf"
# Upstream: process_single.
DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory="6Gi")

FILTER_DIFFERENTIAL_TABLE_PY = r"""
#!/usr/bin/env python

import pandas as pd
import sys
from math import log2

if not sys.argv[1].lower().endswith(('.csv', '.tsv', '.txt')):
    sys.exit("Please provide a .csv, .tsv or .txt file!")

# Determine the separator based on file extension
sep = ',' if sys.argv[1].lower().endswith('.csv') else '\t'

# Read the input file
table = pd.read_csv(sys.argv[1], sep=sep)

# Calculate log2 fold change threshold
logfc_threshold = log2(float(sys.argv[2]))

# define evaluation
def evaluate_condition(x, threshold, cardinality):
    if cardinality == ">=":
        return x >= threshold
    elif cardinality == "<=":
        return x <= threshold
    elif cardinality == ">":
        return x > threshold
    elif cardinality == "<":
        return x < threshold
    else:
        raise ValueError(f"Invalid cardinality: {cardinality}")

# Apply filters
mask = (
    table[sys.argv[4]].notna() &
    table[sys.argv[7]].notna() &
    table[sys.argv[4]].abs().apply(lambda x: evaluate_condition(x, logfc_threshold, sys.argv[3])) &
    table[sys.argv[7]].apply(lambda x: evaluate_condition(x, float(sys.argv[5]), sys.argv[6]))
)
filtered_table = table[mask]
filtered_table_up = filtered_table[filtered_table[sys.argv[4]] > 0]
filtered_table_down = filtered_table[filtered_table[sys.argv[4]] < 0]

# Write the filtered table
filtered_table.to_csv(sys.argv[8] + "_filtered.tsv", sep='\t', index=False)
filtered_table_up.to_csv(sys.argv[8] + "_filtered_up.tsv", sep='\t', index=False)
filtered_table_down.to_csv(sys.argv[8] + "_filtered_down.tsv", sep='\t', index=False)

# Write versions
with open('versions.yml', 'w') as version_file:
    version_file.write('"CUSTOM_FILTERDIFFERENTIALTABLE":\n')
    version_file.write(f"    pandas: {pd.__version__!s}\n")
    version_file.write(f"    python: {sys.version.split()[0]}\n")
""".strip()

MATRIXFILTER_R = r"""
#!/usr/bin/env Rscript

cli_args <- commandArgs(trailingOnly = TRUE)

# Filter rows based on the number of columns passing the abundance threshold. By
# default this will be any row with a value of 1 or more, which would be a
# permissive threshold for RNA-seq data.
#
# In RNA-seq studies it's often not enough to just remove genes not expressed in
# any sample. We also want to remove anything likely to be part of noise, or
# which has sufficiently low expression that differential analysis would not be
# useful. For that reason we might require a higher threshold than 1, and
# require that more than one sample passes.
#
# Often we want to know that a gene is expressed in a substantial enough number
# of sample that differential analysis worthwhile, so we may pick a threshold
# sample number related to group size. Note that we do not filter with an
# awareness of the groups themselves, since this adds bias towards discovery
# between those groups.
#
# This script supports the following options:
#
# - `--minimum_abundance`: Minimum value threshold for feature filtering.
#   Set to a numeric value (e.g., 0.5) to enable, or to 'false'/'null' to disable.
# - `--minimum_samples`: Minimum number of samples passing the abundance threshold.
# - `--minimum_samples_not_na`: Minimum number of non-NA values per feature.
# - `--grouping_variable`: Optional column in sample sheet for group-specific filtering.
# - `--minimum_proportion`: Proportion-based filtering threshold.
# - `--minimum_proportion_not_na`: Minimum proportion of non-NA values required per feature.
# - `--most_variant_features`: Optional integer specifying the number of most-variant features to retain.

################################################
################################################
## Functions                                  ##
################################################
################################################

#' Parse out options from a string without recourse to optparse
#'
#' @param x Long-form argument list like --opt1 val1 --opt2 val2
#'
#' @return named list of options and values similar to optparse

parse_args <- function(x){
    args_list <- unlist(strsplit(x, ' ?--')[[1]])[-1]
    args_vals <- lapply(args_list, function(x) scan(text=x, what='character', quiet = TRUE))

    # Ensure the option vectors are length 2 (key/ value) to catch empty ones
    args_vals <- lapply(args_vals, function(z){ length(z) <- 2; z})

    parsed_args <- structure(lapply(args_vals, function(x) x[2]), names = lapply(args_vals, function(x) x[1]))
    parsed_args[! is.na(parsed_args)]
}

#' Flexibly read CSV or TSV files
#'
#' @param file Input file
#' @param header Passed to read.delim()
#' @param row.names Passed to read.delim()
#' @param nrows Passed to read.delim()
#'
#' @return output Data frame

read_delim_flexible <- function(file, header = TRUE, row.names = NULL, nrows = -1 ){

    ext <- tolower(tail(strsplit(basename(file), split = "\\.")[[1]], 1))

    if (ext == "tsv" || ext == "txt") {
        separator <- "\t"
    } else if (ext == "csv") {
        separator <- ","
    } else {
        stop(paste("Unknown separator for", ext))
    }

    read.delim(
        file,
        sep = separator,
        header = header,
        row.names = row.names,
        check.names = FALSE
    )
}

#' Convert NULL-like inputs to R NULL.
#'
#' @param v Any input to be checked and parsed.
#'
#' @return NULL object if input is NULL-like, otherwise the original input.
parse_null <- function(v) {
    if (is.null(v)) return(NULL)
    if (length(v) == 0) return(NULL)
    vv <- as.character(v)
    if (vv %in% c("null", "NULL", "", "false", "FALSE")) return(NULL)
    v
}

#' Identify rows that are among the top n most variant
#'
#' @param matrix_data Matrix object
#'
#' @return output Boolean vector

most_variant_test <- function(matrix_data) {

    # Determine the indices of the top variant rows based on variance
    top_indices <- order(-apply(matrix_data, 1, var, na.rm = TRUE))[1:opt$most_variant_features]

    # Return a boolean vector indicating if each row is among the top variant ones
    1:nrow(matrix_data) %in% top_indices
}

# Set up default options

opt <- list(
    abundance_matrix_file = cli_args[1],   # Path to input abundance matrix (TSV/CSV)
    sample_file = cli_args[2],           # Path to optional sample metadata file
    sample_id_col = NULL,                   # Column name in samplesheet matching matrix columns
    minimum_abundance = 1,                  # Abundance threshold (disabled if NULL/false)
    minimum_samples = 1,                    # Minimum number of samples passing abundance filter
    minimum_proportion = 0,                 # Proportion-based filtering threshold (optional)
    grouping_variable = NULL,               # Grouping variable for stratified filtering (optional)
    minimum_proportion_not_na = 0.5,        # Minimum proportion of non-NA values per feature
    minimum_samples_not_na = NULL,          # Minimum count of non-NA samples per feature (optional)
    most_variant_features = NULL            # Number of most variant rows to keep (optional)
)
opt_types <- lapply(opt, class)

args_opt <- parse_args(cli_args[3])
for ( ao in names(args_opt)){
    if (! ao %in% names(opt)){
        stop(paste("Invalid option:", ao))
    }else{

        # Allow explicit nulls from Nextflow CLI args
        args_opt[[ao]] <- parse_null(args_opt[[ao]])

        # Preserve classes from defaults where possible (only if not NULL now)
        if (!is.null(args_opt[[ao]]) && !is.null(opt[[ao]])) {
            args_opt[[ao]] <- as(args_opt[[ao]], opt_types[[ao]])
        }
        opt[[ao]] <- args_opt[[ao]]
    }
}

abundance_matrix <- read_delim_flexible(opt$abundance_matrix_file, row.names = 1)
feature_id_name <- colnames( read_delim_flexible(opt$abundance_matrix_file, nrows = 1)[1])

# If a sample sheet was specified, validate the matrix against it

if (opt$sample_file != ''){

    # Read the sample sheet and check against matrix

    samplesheet <- read_delim_flexible(
        opt$sample_file,
        row.names = ifelse(is.null(opt$sample_id_col), 1, opt$sample_id_col)
    )
    missing_samples <- setdiff(rownames(samplesheet), colnames(abundance_matrix))

    if (length(missing_samples) > 0){
        stop(
            paste(
                paste(missing_samples, collapse = ', '),
                'not represented in supplied abundance matrix'
            )
        )
    }else{
        abundance_matrix <- abundance_matrix[,rownames(samplesheet)]
    }
}else{

    # If we're not using a sample sheet to select columns, then at least make
    # sure the ones we have are numeric (some upstream things like the RNA-seq
    # workflow have annotation columns as well)

    numeric_columns <- unlist(lapply(1:ncol(abundance_matrix), function(x) is.numeric(abundance_matrix[,x])))
    abundance_matrix <- abundance_matrix[,numeric_columns]
}

# If we want to define n based on the levels of a grouping variable...

if ((opt$sample_file != '') && ( ! is.null(opt$grouping_variable))){

    # Pick a minimum number of samples to pass threshold based on group size

    if (! opt$grouping_variable %in% colnames(samplesheet)){
        stop(paste(opt$grouping_variable, "not in supplied sample sheet"))
    }else{
        opt$minimum_samples <- min(table(samplesheet[[opt$grouping_variable]]))
        if ( opt$minimum_proportion > 0){
            opt$minimum_samples <- opt$minimum_samples * opt$minimum_proportion
        }
    }
}else if (opt$minimum_proportion > 0){

    # Or if we want to define it based on a static proportion of the sample number

    opt$minimum_samples <- ncol(abundance_matrix) * opt$minimum_proportion
}

# Also set up filtering for NAs; use by default minimum_proportion_not_na; only
# use minimum_samples_not_na if it is provided (default NULL)
# -->NA test can always use minimum_samples_not_na as this will contain the correct
# value even if the proportion is to be used

if (is.null(opt$minimum_samples_not_na)) {
    opt$minimum_samples_not_na <- ncol(abundance_matrix) * opt$minimum_proportion_not_na
}

# Define the tests

tests <- list(
    'na' = function(x) !any(is.na(x)) || sum(!is.na(x)) >= opt$minimum_samples_not_na  # check if enough values in row are not NA
)

# Only apply abundance filter if a threshold was provided
if (!is.null(opt$minimum_abundance)) {
    tests[["abundance"]] = function(x) sum(x >= opt$minimum_abundance, na.rm = TRUE) >= opt$minimum_samples
}

# Apply the functions row-wise on the abundance_matrix and store the result in a boolean matrix

boolean_matrix <- do.call(
    rbind,
    lapply(seq_len(nrow(abundance_matrix)), function(i) {
        vapply(tests, function(f) f(abundance_matrix[i, ]), logical(1))
    })
)

rownames(boolean_matrix) <- rownames(abundance_matrix)
colnames(boolean_matrix) <- names(tests)
# Apply the 'most_variant_test' function to identify the most variant rows and add
# the result to the boolean matrix

if (! is.null(opt$most_variant_features)) {
    most_variant_vectors <- most_variant_test(abundance_matrix)
    boolean_matrix <- cbind(boolean_matrix, most_variant_vectors)
}

# We will retain features passing all tests

keep <- apply(boolean_matrix, 1, all)

# Write out the matrix retaining the specified rows and re-prepending the
# column with the feature identifiers

prefix = ifelse(cli_args[4] == 'null', cli_args[4], cli_args[4])

thresholds <- data.frame(
    rule = c(
        "minimum_samples_not_na",
        if (!is.null(opt$minimum_abundance)) c("minimum_abundance", "minimum_samples"),
        if (!is.null(opt$most_variant_features)) "most_variant_features"
    ),
    threshold = c(
        opt$minimum_samples_not_na,
        if (!is.null(opt$minimum_abundance)) c(opt$minimum_abundance, opt$minimum_samples),
        if (!is.null(opt$most_variant_features)) opt$most_variant_features
    ),
    check.names = FALSE
)

    write.table(
        thresholds,
        file = paste0(prefix, ".thresholds.tsv"),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )


write.table(
    data.frame(rownames(abundance_matrix)[keep], abundance_matrix[keep,,drop = FALSE]),
    file = paste0(
        prefix,
        '.filtered.tsv'
    ),
    col.names = c(feature_id_name, colnames(abundance_matrix)),
    row.names = FALSE,
    sep = '\t',
    quote = FALSE
)

# Write a boolean matrix returning specifying the status of each test

write.table(
    data.frame(rownames(abundance_matrix), boolean_matrix),
    file = paste0(prefix, '.tests.tsv'),
    col.names = c(feature_id_name, colnames(boolean_matrix)),
    row.names = FALSE,
    sep = '\t',
    quote = FALSE
)

################################################
################################################
## R SESSION INFO                             ##
################################################
################################################

sink("R_sessionInfo.log")
print(sessionInfo())
sink()

################################################
################################################
## VERSIONS FILE                              ##
################################################
################################################

r.version <- strsplit(version[['version.string']], ' ')[[1]][3]

writeLines(
    c(
        '"CUSTOM_MATRIXFILTER":',
        paste('    r-base:', r.version)
    ),
'versions.yml')

################################################
################################################
################################################
################################################
""".strip()  # noqa: E501

custom_filterdifferentialtable_cmd = shell.create(
    name="custom_filterdifferentialtable",
    image=PANDAS_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "table": File,
        "prefix": str,
        "logfc_column": str,
        "fc_threshold": str,
        "fc_cardinality": str,
        "stat_column": str,
        "stat_threshold": str,
        "stat_cardinality": str,
    },
    outputs={"results": Dir},
    script=r"""
        TABLE=({inputs.table})
        P={inputs.prefix}
        OUT={outputs.results}
        W=$(mktemp -d); cd "$W"
        NAME=$(basename "${TABLE[0]}"); ln -s "${TABLE[0]}" "$NAME"
        cat > filter_differential_table.py <<'END_OF_SCRIPT'
FILTER_DIFFERENTIAL_TABLE_PY
END_OF_SCRIPT
        python filter_differential_table.py "$NAME" {inputs.fc_threshold} {inputs.fc_cardinality} \
            {inputs.logfc_column} {inputs.stat_threshold} {inputs.stat_cardinality} {inputs.stat_column} "$P" \
            || exit 1
        mv "${P}_filtered.tsv" "${P}_filtered_up.tsv" "${P}_filtered_down.tsv" versions.yml "$OUT"/
    """.replace("FILTER_DIFFERENTIAL_TABLE_PY", FILTER_DIFFERENTIAL_TABLE_PY),
)

custom_matrixfilter_cmd = shell.create(
    name="custom_matrixfilter",
    image=R_BASE_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"abundance": File, "samplesheet": File | None, "prefix": str, "args": str},
    defaults={"args": "null"},
    outputs={"results": Dir},
    script=r"""
        ABUNDANCE=({inputs.abundance})
        shopt -s nullglob; SHEET=({inputs.samplesheet}); shopt -u nullglob
        P={inputs.prefix}
        OUT={outputs.results}
        W=$(mktemp -d); cd "$W"
        NAME=$(basename "${ABUNDANCE[0]}"); ln -s "${ABUNDANCE[0]}" "$NAME"
        SHEET_NAME=""
        if [ ${#SHEET[@]} -gt 0 ]; then SHEET_NAME=$(basename "${SHEET[0]}"); ln -s "${SHEET[0]}" "$SHEET_NAME"; fi
        cat > matrixfilter.R <<'END_OF_SCRIPT'
MATRIXFILTER_R
END_OF_SCRIPT
        Rscript matrixfilter.R "$NAME" "$SHEET_NAME" {inputs.args} "$P" || exit 1
        mv "$P.filtered.tsv" "$P.tests.tsv" "$P.thresholds.tsv" R_sessionInfo.log versions.yml "$OUT"/
    """.replace("MATRIXFILTER_R", MATRIXFILTER_R),
)

custom_tabulartogseagct_cmd = shell.create(
    name="custom_tabulartogseagct",
    image=COREUTILS_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"tabular": File, "prefix": str, "separator": str},
    outputs={"results": Dir},
    script=r"""
        TABULAR=({inputs.tabular})
        tabular="${TABULAR[0]}"
        separator={inputs.separator}
        W=$(mktemp -d); cd "$W"
        n_columns=$(head -n 1 ${tabular} | tr "${separator}" "\n" | wc -l)
        n_lines=$(wc -l < ${tabular})
        gct_file={inputs.prefix}.gct

        echo -e "#1.2$(printf '\t%.0s' {1..$n_columns})\n$((n_lines-1))\t$((n_columns-1))$(printf '\t%.0s' {1..$((n_columns-1))})" > $gct_file
        echo -e "NAME\tDESCRIPTION\t$(head -n 1 ${tabular} | cut -f1 -d"$(printf '%b' "${separator}")" --complement | awk -F"${separator}" 'BEGIN { OFS = "\t"}; {$1=$1}1' )" >> $gct_file
        cat ${tabular} | tail -n +2 | awk -F"${separator}" 'BEGIN { OFS = "\t"} {$1=$1"\tNA"}1' >> $gct_file
        mv $gct_file {outputs.results}/
    """,  # noqa: E501
)

custom_tabulartogseacls_cmd = shell.create(
    name="custom_tabulartogseacls",
    image=COREUTILS_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"samples": File, "prefix": str, "variable": str, "separator": str},
    outputs={"results": Dir},
    script=r"""
        SAMPLES=({inputs.samples})
        samples="${SAMPLES[0]}"
        separator={inputs.separator}
        variable={inputs.variable}
        W=$(mktemp -d); cd "$W"
        cls_file={inputs.prefix}.cls

        column_number=$(cat ${samples} | head -n 1 | tr "${separator}" "\n" | grep -En "^${variable}$" | awk -F':' '{print $1}')
        classes=$(tail -n +2 ${samples} | awk -F"${separator}" '{print $'$column_number'}' | sed 's/^$/empty/g')
        unique_classes=$(echo -e "$classes" | awk '!x[$0]++')

        echo -e "$(echo -e "$classes" | wc -l) $(echo -e "$unique_classes" | wc -l) 1" > $cls_file
        echo -e "#$(echo -e "$unique_classes" | tr '\n' ' ')" | sed "s/ $//" >> $cls_file
        echo -e "$classes" | tr '\n' ' ' | sed "s/ $//" >> $cls_file
        echo -e "\n" >> $cls_file
        mv $cls_file {outputs.results}/
    """,  # noqa: E501
)

custom_tabulartogseachip_cmd = shell.create(
    name="custom_tabulartogseachip",
    image=GAWK_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"tabular": File, "prefix": str, "id": str, "symbol": str},
    outputs={"results": Dir},
    script=r"""
        TABULAR=({inputs.tabular})
        tabular="${TABULAR[0]}"
        id={inputs.id}
        symbol={inputs.symbol}
        W=$(mktemp -d); cd "$W"
        function find_column_number {
            file=$1
            column=$2

            head -n 1 $file | tr '\t' '\n' | grep -n "^${column}$" | awk -F':' '{print $1}'
        }

        id_col=$(find_column_number $tabular $id)
        symbol_col=$(find_column_number $tabular $symbol)
        outfile={inputs.prefix}.chip

        echo -e "Probe Set ID\tGene Symbol\tGene Title" > ${outfile}.tmp
        tail -n +2 $tabular | awk -F'\t' -v id=$id_col -v symbol=$symbol_col '{print $id"\t"$symbol"\tNA"}' >> ${outfile}.tmp
        mv ${outfile}.tmp ${outfile}
        mv ${outfile} {outputs.results}/
    """,  # noqa: E501
)

# One image per environment: the four tools' images differ.
custom_pandas_env = flyte.TaskEnvironment.from_task("custom_pandas", custom_filterdifferentialtable_cmd.as_task())
custom_r_env = flyte.TaskEnvironment.from_task("custom_r", custom_matrixfilter_cmd.as_task())
custom_coreutils_env = flyte.TaskEnvironment.from_task(
    "custom_coreutils", custom_tabulartogseagct_cmd.as_task(), custom_tabulartogseacls_cmd.as_task()
)
custom_gawk_env = flyte.TaskEnvironment.from_task("custom_gawk", custom_tabulartogseachip_cmd.as_task())
env = flyte.TaskEnvironment(
    name="custom", depends_on=[custom_pandas_env, custom_r_env, custom_coreutils_env, custom_gawk_env]
)


async def pick(results: Dir, name: str, tool: str) -> File:
    found = await results.get_file(name)
    if found is None:
        raise FileNotFoundError(f"{tool} wrote no {name}")
    return found


@dataclass
class FilterDifferentialTableResult:
    filtered: File  # <prefix>_filtered.tsv
    filtered_up: File  # <prefix>_filtered_up.tsv
    filtered_down: File  # <prefix>_filtered_down.tsv
    versions: File  # versions.yml


async def custom_filterdifferentialtable(
    table: File,
    prefix: str,
    logfc_column: str,
    fc_threshold: float,
    fc_cardinality: str,
    stat_column: str,
    stat_threshold: float,
    stat_cardinality: str,
) -> FilterDifferentialTableResult:
    """Rows of ``table`` whose |``logfc_column``| and ``stat_column`` pass their thresholds (``>=``, ``<=``, …)."""
    results = await custom_filterdifferentialtable_cmd(
        table=table,
        prefix=prefix,
        logfc_column=logfc_column,
        fc_threshold=str(fc_threshold),
        fc_cardinality=fc_cardinality,
        stat_column=stat_column,
        stat_threshold=str(stat_threshold),
        stat_cardinality=stat_cardinality,
    )
    tool = "custom_filterdifferentialtable"
    return FilterDifferentialTableResult(
        filtered=await pick(results, f"{prefix}_filtered.tsv", tool),
        filtered_up=await pick(results, f"{prefix}_filtered_up.tsv", tool),
        filtered_down=await pick(results, f"{prefix}_filtered_down.tsv", tool),
        versions=await pick(results, "versions.yml", tool),
    )


@dataclass
class MatrixFilterResult:
    filtered: File  # <prefix>.filtered.tsv
    tests: File  # <prefix>.tests.tsv
    thresholds: File  # <prefix>.thresholds.tsv
    session_info: File  # R_sessionInfo.log
    versions: File  # versions.yml


async def custom_matrixfilter(
    abundance: File, prefix: str, samplesheet: File | None = None, args: str = ""
) -> MatrixFilterResult:
    """Features of ``abundance`` passing the tests ``args`` set; ``samplesheet`` selects and orders the samples."""
    # As upstream's template renders an unset ext.args.
    results = await custom_matrixfilter_cmd(
        abundance=abundance, samplesheet=samplesheet, prefix=prefix, args=args or "null"
    )
    tool = "custom_matrixfilter"
    return MatrixFilterResult(
        filtered=await pick(results, f"{prefix}.filtered.tsv", tool),
        tests=await pick(results, f"{prefix}.tests.tsv", tool),
        thresholds=await pick(results, f"{prefix}.thresholds.tsv", tool),
        session_info=await pick(results, "R_sessionInfo.log", tool),
        versions=await pick(results, "versions.yml", tool),
    )


def default_separator(name: str, comma_suffix: str) -> str:
    """As upstream: a comma for files ending ``comma_suffix``, else a tab (as ``\\t``, for tr and awk)."""
    return "," if name.endswith(comma_suffix) else "\\t"


async def custom_tabulartogseagct(tabular: File, prefix: str, separator: str | None = None) -> File:
    """``<prefix>.gct`` from a matrix whose first column holds the feature IDs (a ``.csv`` is comma-separated)."""
    sep = "\\t" if separator == "\t" else separator or default_separator(tabular.name or tabular.path, ".csv")
    results = await custom_tabulartogseagct_cmd(tabular=tabular, prefix=prefix, separator=sep)
    return await pick(results, f"{prefix}.gct", "custom_tabulartogseagct")


async def custom_tabulartogseacls(samples: File, prefix: str, variable: str, separator: str | None = None) -> File:
    """``<prefix>.cls`` with each sample's ``variable`` value (blank: ``empty``); a ``.tsv`` is tab-separated."""
    name = samples.name or samples.path
    sep = "\\t" if separator == "\t" else separator or ("\\t" if name.endswith(".tsv") else ",")
    results = await custom_tabulartogseacls_cmd(samples=samples, prefix=prefix, variable=variable, separator=sep)
    return await pick(results, f"{prefix}.cls", "custom_tabulartogseacls")


async def custom_tabulartogseachip(tabular: File, prefix: str, id: str, symbol: str) -> File:
    """``<prefix>.chip``: the tab-separated ``tabular``'s ``id`` and ``symbol`` columns."""
    results = await custom_tabulartogseachip_cmd(tabular=tabular, prefix=prefix, id=id, symbol=symbol)
    return await pick(results, f"{prefix}.chip", "custom_tabulartogseachip")
