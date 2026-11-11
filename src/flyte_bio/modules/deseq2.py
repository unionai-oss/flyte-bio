"""deseq2 — differential expression of a count matrix with DESeq2.

Exposes :func:`deseq2_differential`, which tests one contrast: either
``contrast_variable`` with its ``reference`` and ``target`` levels, or a
model ``formula`` with a ``comparison`` (a coefficient name, or a limma-style
contrast expression). Optional inputs: control genes (e.g. ERCC spikes;
stripped from the matrix, or used for size factors with
``--sizefactors_from_controls TRUE``) and transcript lengths (as DESeq2
normalisation factors). Outputs keep upstream's names:

- ``<prefix>.deseq2.results.tsv``, ``<prefix>.deseq2.model.txt``;
- ``<prefix>.deseq2.dispersion.{png,pdf}``, ``<prefix>.dds.rld.rds``;
- ``<prefix>.deseq2.sizefactors.tsv``, ``<prefix>.normalised_counts.tsv``;
- ``<prefix>.rlog.tsv`` and/or ``<prefix>.vst.tsv`` (per ``--vs_method``);
- ``<prefix>.R_sessionInfo.log`` and ``versions.yml``.

``args`` are upstream's ``ext.args``: ``--<option> <value>`` pairs (e.g.
``--sample_id_col sample --vs_method rlog --round_digits 5``). The R code is
upstream's template, unchanged but for taking the templated values as
command-line arguments; the prefix stands for both the templated
``ext.prefix`` and sample ID. As upstream renders unset values, a missing
contrast field is passed as ``null`` and a missing file as an empty string. On failure, the end
of R's error output becomes the task's error message.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

DESEQ2_IMAGE = (
    "community.wave.seqera.io/library/bioconductor-deseq2_bioconductor-limma_bioconductor-apeglm:7fccb58b6c45d084"
)
# Upstream: process_single.
DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory="8Gi")

DESEQ2_DIFFERENTIAL_R = r"""
#!/usr/bin/env Rscript

cli_args <- commandArgs(trailingOnly = TRUE)

################################################
################################################
## Functions                                  ##
################################################
################################################

#' Check for Non-Empty, Non-Whitespace String
#'
#' This function checks if the input is non-NULL and contains more than just whitespace.
#' It returns TRUE if the input is a non-empty, non-whitespace string, and FALSE otherwise.
#'
#' @param input A variable to check.
#' @return A logical value: TRUE if the input is a valid, non-empty, non-whitespace string; FALSE otherwise.
#' @examples
#' is_valid_string("Hello World") # Returns TRUE
#' is_valid_string("   ")         # Returns FALSE
#' is_valid_string(NULL)          # Returns FALSE

is_valid_string <- function(input) {
  !is.null(input) && nzchar(trimws(input))
}

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
#'
#' @return output Data frame

read_delim_flexible <- function(file, header = TRUE, row.names = NULL, check.names = TRUE){

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
    check.names = check.names
  )
}

#
#' Turn “null” or empty strings into actual NULL
#'
#' @param x Input option
#'
#' @return NULL or x
#'
nullify <- function(x) {
  if (is.character(x) && (tolower(x) == "null" || x == "")) NULL else x
}

#' Run DESeq2 differential expression results with shared settings
#'
#' A thin wrapper around DESeq2::results() that injects the common parameters
#' (lfcThreshold, altHypothesis, independentFiltering, alpha, pAdjustMethod, minmu)
#' from the global `opt` list and the `dds` object. Use the `...` to supply the
#' branch-specific argument (`name = …` or `contrast = …`).
#'
#' @param ... Additional arguments passed to `results()`, e.g. `name = opt$contrast_string`
#'   or `contrast = c(variable, target, reference)`.
#' @return A `DESeqResults` object containing the differential expression test results.
#'
run_results <- function(...) {
  results(
    dds,
    lfcThreshold         = opt$lfc_threshold,
    altHypothesis        = opt$alt_hypothesis,
    independentFiltering = opt$independent_filtering,
    alpha                = opt$alpha,
    pAdjustMethod        = opt$p_adjust_method,
    minmu                = opt$minmu,
    ...
  )
}
#' Shrink log2 fold changes with shared settings
#'
#' A thin wrapper around DESeq2::lfcShrink() that injects the common parameters
#' and the `dds` object. Use the `...` to supply the branch-specific argument
#' (`coef = …` or `contrast = …`).
#'
#' @param ... Additional arguments passed to `lfcShrink()`, e.g. `coef = opt$contrast_string`
#'   or `contrast = c(variable, target, reference)`.
#' @param shrink_type Type of shrinkage estimator ('ashr' or 'apeglm').
#' @param contrast_name Optional label of the contrast for error messages.
#' @return A `DESeqResults` object with shrunken log2 fold changes.
run_shrink <- function(..., shrink_type, contrast_name = NULL) {
  dots <- list(...)
  if (shrink_type == "apeglm" && "contrast" %in% names(dots)) {
    if (is.null(contrast_name)) {
      contrast_name <- paste(dots$contrast, collapse = ", ")
    }
    stop(
      "apeglm shrinkage requires a model coefficient (coef=), but contrast '",
      contrast_name,
      "' is not one. Use a formula with an intercept and a coefficient name from resultsNames() ",
      "as the comparison, or use --shrink_lfc_type ashr."
    )
  }
  lfcShrink(
    dds,
    type = shrink_type,
    ...
  )
}

################################################
################################################
## PARSE PARAMETERS FROM NEXTFLOW             ##
################################################
################################################

# I've defined these in a single array like this so that we could go back to an
# optparse-driven method in future with module bin/ directories, rather than
# the template

# Set defaults and classes

opt <- list(
  output_prefix = ifelse(cli_args[1] == 'null', cli_args[1], cli_args[1]),
  count_file = cli_args[2],
  sample_file = cli_args[3],
  contrast_variable = cli_args[4],
  reference_level = cli_args[5],
  target_level = cli_args[6],
  contrast_string    = cli_args[7],             # Optional full (complex) contrast expression comparison
  formula            = cli_args[8],              # User-specified formula
  blocking_variables = NULL,
  control_genes_file = cli_args[9],
  transcript_lengths_file = cli_args[10],
  sizefactors_from_controls = FALSE,
  gene_id_col = "gene_id",
  sample_id_col = "experiment_accession",
  subset_to_contrast_samples = FALSE,
  exclude_samples_col = NULL,
  exclude_samples_values = NULL,
  test = "Wald",
  fit_type = "parametric",
  sf_type = 'ratio',
  min_replicates_for_replace = 7,
  use_t = FALSE,
  lfc_threshold = 0,
  alt_hypothesis = 'greaterAbs',
  independent_filtering = TRUE,
  p_adjust_method = 'BH',
  alpha = 0.1,
  minmu = 0.5,
  vs_method = 'vst', # 'rlog', 'vst', or 'rlog,vst'
  shrink_lfc = TRUE,
  shrink_lfc_type = 'ashr',
  cores = 1,
  vs_blind = TRUE,
  vst_nsub = 1000,
  round_digits = NULL,
  seed = NULL
)
opt_types <- lapply(opt, class)

# Apply parameter overrides

args_opt <- parse_args(cli_args[11])
for ( ao in names(args_opt)){
  if (! ao %in% names(opt)){
    stop(paste("Invalid option:", ao))
  }else{

    # Preserve classes from defaults where possible
    if (! is.null(opt[[ao]])){
      args_opt[[ao]] <- as(args_opt[[ao]], opt_types[[ao]])
    }
    opt[[ao]] <- args_opt[[ao]]
  }
}


if ( ! is.null(opt$round_digits)){
  opt$round_digits <- as.numeric(opt$round_digits)
}

# If there is no option supplied, convert string "null" to NULL
keys <- c("formula", "contrast_string", "contrast_variable", "reference_level", "target_level", "seed", "blocking_variable", "transcript_lengths_file")
opt[keys] <- lapply(opt[keys], nullify)

if ( ! is.null(opt$seed)){
  set.seed(opt$seed)
}

# Check if required parameters have been provided
if (is_valid_string(opt$formula)) {
  required_opts <- c('output_prefix', 'contrast_string')
} else {
  required_opts <- c('contrast_variable', 'reference_level', 'target_level', 'output_prefix')
}
missing <- required_opts[!unlist(lapply(opt[required_opts], is_valid_string)) | !required_opts %in% names(opt)]

if (length(missing) > 0){
  stop(paste("Missing required options:", paste(missing, collapse=', ')))
}

if (opt$shrink_lfc && opt$shrink_lfc_type == 'apeglm' && is.null(opt$contrast_string)) {
  stop(
    "apeglm shrinkage requires a model coefficient, but variable/reference/target contrasts use a ",
    "'~ 0 + ...' design with no such coefficient. Use a formula with an intercept and a coefficient ",
    "name from resultsNames() as the comparison, or use --shrink_lfc_type ashr."
  )
}

# Check file inputs are valid

for (file_input in c('count_file', 'sample_file')){
  if (! is_valid_string(opt[[file_input]])) {
    stop(paste("Please provide", file_input), call. = FALSE)
  }

  if (! file.exists(opt[[file_input]])){
    stop(paste0('Value of ', file_input, ': ', opt[[file_input]], ' is not a valid file'))
  }
}

################################################
################################################
## Finish loading libraries                   ##
################################################
################################################

library(DESeq2)
library(BiocParallel)
library(limma)

################################################
################################################
## READ IN COUNTS FILE AND SAMPLE METADATA    ##
################################################
################################################

count.table <-
  read_delim_flexible(
    file = opt$count_file,
    header = TRUE,
    row.names = opt$gene_id_col,
    check.names = FALSE
  )
sample.sheet <- read_delim_flexible(file = opt$sample_file)

# Deal with spaces that may be in sample column
opt$sample_id_col <- make.names(opt$sample_id_col)

if (! opt$sample_id_col %in% colnames(sample.sheet)){
  stop(paste0("Specified sample ID column '", opt$sample_id_col, "' is not in the sample sheet"))
}

# Sample sheet can have duplicate rows for multiple sequencing runs, so uniqify
# before assigning row names

sample.sheet <- sample.sheet[! duplicated(sample.sheet[[opt$sample_id_col]]), ]
rownames(sample.sheet) <- sample.sheet[[opt$sample_id_col]]

# Check that all samples specified in the input sheet are present in the counts
# table. Assuming they are, subset and sort the count table to match the sample
# sheet

missing_samples <-
  sample.sheet[!rownames(sample.sheet) %in% colnames(count.table), opt$sample_id_col]

if (length(missing_samples) > 0) {
  stop(paste(
    length(missing_samples),
    'specified samples missing from count table:',
    paste(missing_samples, collapse = ',')
  ))
} else{
  # Save any non-count data, will gene metadata etc we might need later
  noncount.table <-
    count.table[, !colnames(count.table) %in% rownames(sample.sheet), drop = FALSE]
  count.table <- count.table[, rownames(sample.sheet)]
}

################################################
################################################
## CHECK CONTRAST SPECIFICATION               ##
################################################
################################################
if (! is_valid_string(opt$formula)) {
  contrast_variable <- make.names(opt$contrast_variable)
  blocking.vars <- c()

  if (!contrast_variable %in% colnames(sample.sheet)) {
    stop(
      paste0(
        'Chosen contrast variable \"',
        contrast_variable,
        '\" not in sample sheet'
      )
    )
  } else if (any(!c(opt$reference_level, opt$target_level) %in% sample.sheet[[contrast_variable]])) {
    stop(
      paste(
        'Please choose reference and treatment levels that are present in the',
        contrast_variable,
        'column of the sample sheet'
      )
    )
  } else if (is_valid_string(opt$blocking_variables)) {
    blocking.vars = make.names(unlist(strsplit(opt$blocking_variables, split = ';')))
    if (!all(blocking.vars %in% colnames(sample.sheet))) {
      missing_block <- paste(blocking.vars[! blocking.vars %in% colnames(sample.sheet)], collapse = ',')
      stop(
        paste(
          'Blocking variables', missing_block,
          'do not correspond to sample sheet columns.'
        )
      )
    }
  }

  # Optionally, subset to only the samples involved in the contrast
  if (opt$subset_to_contrast_samples){
    sample_selector <- sample.sheet[[contrast_variable]] %in% c(opt$target_level, opt$reference_level)
    selected_samples <- sample.sheet[sample_selector, opt$sample_id_col]
    count.table <- count.table[, selected_samples]
    sample.sheet <- sample.sheet[selected_samples, ]
  }

}
# Optionally, remove samples with specified values in a given field (probably
# don't use this as well as the above)

if ((is_valid_string(opt$exclude_samples_col)) && (is_valid_string(opt$exclude_samples_values))){
  exclude_values = unlist(strsplit(opt$exclude_samples_values, split = ';'))

  if (! opt$exclude_samples_col %in% colnames(sample.sheet)){
    stop(paste(opt$exclude_samples_col, ' specified to subset samples is not a valid sample sheet column'))
  }

  print(paste0('Excluding samples with values of ', opt$exclude_samples_values, ' in ', opt$exclude_samples_col))
  sample_selector <- ! sample.sheet[[opt$exclude_samples_col]] %in% exclude_values

  selected_samples <- sample.sheet[sample_selector, opt$sample_id_col]
  count.table <- count.table[, selected_samples]
  sample.sheet <- sample.sheet[selected_samples, ]
}

# Now specify the model. Use cell-means style so we can be explicit with the
# contrasts

if (is_valid_string(opt$formula)) {
  message("Using user-specified formula: ", opt$formula)
  model  <- opt$formula
}  else {
  model <- '~ 0'

  if (is_valid_string(opt$blocking_variables)) {
    model <- paste(model, paste(blocking.vars, collapse = ' + '), sep=' + ')
  }

  # Make sure all the appropriate variables are factors

  for (v in c(blocking.vars, contrast_variable)) {
    sample.sheet[[v]] <- as.factor(sample.sheet[[v]])
  }

  # Variable of interest goes last, see
  # https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#multi-factor-designs

  model <- paste(model, contrast_variable, sep = ' + ')
}
message("Final design formula: ", model)

################################################
################################################
## Run DESeq2 processes                       ##
################################################
################################################

if (opt$control_genes_file != ''){
  control_genes <- readLines(opt$control_genes_file)
  if (! opt$sizefactors_from_controls){
    count.table <- count.table[setdiff(rownames(count.table), control_genes),]
  }
}

dds <- DESeqDataSetFromMatrix(
  countData = round(count.table),
  colData = sample.sheet,
  design = as.formula(model)
)

# Build in transcript lengths. Copying what tximport does here:
# https://github.com/thelovelab/DESeq2/blob/6947d5bc629015fb8ffb2453a91b71665a164483/R/AllClasses.R#L409

if (!is.null(opt$transcript_lengths_file)){
  lengths <-
    read_delim_flexible(
      file = opt$transcript_lengths_file,
      header = TRUE,
      row.names = opt$gene_id_col,
      check.names = FALSE
    )
  lengths <- lengths[rownames(count.table), colnames(count.table)]
  dimnames(lengths) <- dimnames(dds)
  assays(dds)[["avgTxLength"]] <- lengths
}

if (opt$control_genes_file != '' && opt$sizefactors_from_controls){
  print(paste('Estimating size factors using', length(control_genes), 'control genes'))
  dds <- estimateSizeFactors(dds, controlGenes=rownames(count.table) %in% control_genes)
}

dds <- DESeq(
  dds,
  test = opt$test,
  fitType = opt$fit_type,
  minReplicatesForReplace = opt$min_replicates_for_replace,
  useT = opt$use_t,
  sfType = opt$sf_type,
  parallel=TRUE, BPPARAM=MulticoreParam(opt$cores)
)

if (!is.null(opt$contrast_string)) {
  coef_names <- resultsNames(dds)
  if (opt$contrast_string %in% coef_names) {
    # Direct coefficient name
    comp.results <- run_results(name = opt$contrast_string)
    if (opt$shrink_lfc) {
      comp.results <- run_shrink(coef = opt$contrast_string, shrink_type = opt$shrink_lfc_type)
    }
  } else {
    # Parse as limma-style contrast expression
    design_mat <- model.matrix(as.formula(model), data = as.data.frame(colData(dds)))
    colnames(design_mat) <- make.names(colnames(design_mat))
    numeric_contrast <- as.numeric(
      limma::makeContrasts(contrasts = opt$contrast_string, levels = colnames(design_mat))
    )

    # Run DESeq2 results with numeric contrast
    comp.results <- run_results(contrast = numeric_contrast)
    if (opt$shrink_lfc) {
      comp.results <- run_shrink(
        contrast      = numeric_contrast,
        shrink_type   = opt$shrink_lfc_type,
        contrast_name = opt$contrast_string
      )
    }
  }
} else {
  contrast_var_tg_ref <- c(contrast_variable,
                           opt$target_level,
                           opt$reference_level)

  comp.results <- run_results(contrast = contrast_var_tg_ref)
  if (opt$shrink_lfc) {
    comp.results <- run_shrink(contrast = contrast_var_tg_ref, shrink_type = opt$shrink_lfc_type)
  }
}

# See https://support.bioconductor.org/p/97676/

if (!is.null(opt$transcript_lengths_file)){
  size_factors = estimateSizeFactorsForMatrix(counts(dds) / assays(dds)[["avgTxLength"]])
} else {
  size_factors = sizeFactors(dds)
}

################################################
################################################
## Generate outputs                           ##
################################################
################################################

# if no explicit contrast_string was given, synthesize a two‐level default
if (! is_valid_string(opt$contrast_string)) {
  opt$contrast_string <- paste(opt$target_level, opt$reference_level, sep = "_vs_")
}

cat("Saving results for ", opt$contrast_string, " ...\n", sep = "")

# Differential expression table- note very limited rounding for consistency of
# results

if (! is.null(opt$round_digits)){
  comp.results <- apply(data.frame(comp.results), 2, function(x) round(x, opt$round_digits))
}
comp.results <- `colnames<-`(
  data.frame(
    gene_id = rownames(comp.results),
    comp.results,
    check.names = FALSE
  ),
  c(opt$gene_id_col, colnames(comp.results))  # Setting all column names
)

write.table(
  comp.results,
  file = paste(opt$output_prefix, 'deseq2.results.tsv', sep = '.'),
  col.names = TRUE,
  row.names = FALSE,
  sep = '\t',
  quote = FALSE
)

# Dispersion plot

png(
  file = paste(opt$output_prefix, 'deseq2.dispersion.png', sep = '.'),
  width = 600,
  height = 600
)
plotDispEsts(dds)
dev.off()

pdf(
  file = paste(opt$output_prefix, 'deseq2.dispersion.pdf', sep = '.')
)
plotDispEsts(dds)
dev.off()

# R object for other processes to use

saveRDS(dds, file = paste(opt$output_prefix, 'dds.rld.rds', sep = '.'))

# Size factors

sf_df = data.frame(
  sample = names(size_factors),
  data.frame(size_factors, check.names = FALSE),
  check.names = FALSE
)
colnames(sf_df) <- c('sample', 'sizeFactor')
write.table(
  sf_df,
  file = paste(opt$output_prefix, 'deseq2.sizefactors.tsv', sep = '.'),
  col.names = TRUE,
  row.names = FALSE,
  sep = '\t',
  quote = FALSE
)

# Write specified matrices

normalised_matrix <- counts(dds, normalized = TRUE)
if (! is.null(opt$round_digits)){
  normalised_matrix <- apply(normalised_matrix, 2, function(x) round(x, opt$round_digits))
}
normalised_matrix <- `colnames<-`(
  data.frame(
    gene_id = rownames(counts(dds)),  # First column with row names from counts(dds)
    normalised_matrix,                # Other columns
    check.names = FALSE
  ),
  c(opt$gene_id_col, colnames(normalised_matrix))  # Setting all column names
)

write.table(
  normalised_matrix,
  file = paste(opt$output_prefix, 'normalised_counts.tsv', sep = '.'),
  col.names = TRUE,
  row.names = FALSE,
  sep = '\t',
  quote = FALSE
)

# Note very limited rounding for consistency of results

for (vs_method_name in strsplit(opt$vs_method, ',')){
  if (vs_method_name == 'vst'){
    vs_mat <- assay(vst(dds, blind = opt$vs_blind, nsub = opt$vst_nsub))
  }else if (vs_method_name == 'rlog'){
    vs_mat <- assay(rlog(dds, blind = opt$vs_blind, fitType = opt$fit_type))
  }

  if (! is.null(opt$round_digits)){
    vs_mat <- apply(vs_mat, 2, function(x) round(x, opt$round_digits))
  }

  vs_mat <- `colnames<-`(
    data.frame(
      gene_id = rownames(counts(dds)),  # First column with row names from counts(dds)
      vs_mat,                           # Other columns from vs_mat
      check.names = FALSE
    ),
    c(opt$gene_id_col, colnames(vs_mat))  # Setting all column names
  )

  write.table(
    vs_mat,
    file = paste(opt$output_prefix, vs_method_name, 'tsv', sep = '.'),
    col.names = TRUE,
    row.names = FALSE,
    sep = '\t',
    quote = FALSE
  )
}

# Save model to file

write(model, file=paste(opt$output_prefix, 'deseq2.model.txt', sep = '.'))

################################################
################################################
## R SESSION INFO                             ##
################################################
################################################

sink(paste(opt$output_prefix, "R_sessionInfo.log", sep = '.'))
print(sessionInfo())
sink()

################################################
################################################
## VERSIONS FILE                              ##
################################################
################################################

r.version <- paste(R.version[['major']],R.version[['minor']], sep = ".")
deseq2.version <- as.character(packageVersion('DESeq2'))

writeLines(
  c(
    '"DESEQ2_DIFFERENTIAL":',
    paste('    r-base:', r.version),
    paste('    bioconductor-deseq2:', deseq2.version)
  ),
  'versions.yml')

################################################
################################################
################################################
################################################
""".strip()  # noqa: E501, RUF001

deseq2_differential_cmd = shell.create(
    name="deseq2_differential",
    image=DESEQ2_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "samplesheet": File,
        "counts": File,
        "control_genes": File | None,
        "transcript_lengths": File | None,
        "prefix": str,
        "contrast_variable": str,
        "reference": str,
        "target": str,
        "formula": str,
        "comparison": str,
        "args": str,
    },
    outputs={"results": Dir},
    script=r"""
        SHEET=({inputs.samplesheet})
        COUNTS=({inputs.counts})
        shopt -s nullglob; CONTROLS=({inputs.control_genes}); LENGTHS=({inputs.transcript_lengths}); shopt -u nullglob
        P={inputs.prefix}
        OUT={outputs.results}
        W=$(mktemp -d); cd "$W"
        stage() { if [ $# -gt 0 ]; then ln -s "$1" "$(basename "$1")"; basename "$1"; fi; }
        SHEET_NAME=$(stage "${SHEET[@]}")
        COUNTS_NAME=$(stage "${COUNTS[@]}")
        CONTROLS_NAME=$(stage "${CONTROLS[@]}")
        LENGTHS_NAME=$(stage "${LENGTHS[@]}")
        cat > deseq2_differential.R <<'END_OF_SCRIPT'
DESEQ2_DIFFERENTIAL_R
END_OF_SCRIPT
        Rscript deseq2_differential.R "$P" "$COUNTS_NAME" "$SHEET_NAME" {inputs.contrast_variable} \
            {inputs.reference} {inputs.target} {inputs.comparison} {inputs.formula} "$CONTROLS_NAME" \
            "$LENGTHS_NAME" {inputs.args} 2> >(tee r.err >&2) || { tail -n 20 r.err > /dev/termination-log; exit 1; }
        mv "$P".deseq2.* "$P".dds.rld.rds "$P".normalised_counts.tsv "$P".R_sessionInfo.log versions.yml "$OUT"/
        for f in "$P".rlog.tsv "$P".vst.tsv; do if [ -e "$f" ]; then mv "$f" "$OUT"/; fi; done
    """.replace("DESEQ2_DIFFERENTIAL_R", DESEQ2_DIFFERENTIAL_R),
)

env = flyte.TaskEnvironment.from_task("deseq2", deseq2_differential_cmd.as_task())


@dataclass
class Deseq2DifferentialResult:
    results: File  # <prefix>.deseq2.results.tsv
    dispersion_plot_png: File  # <prefix>.deseq2.dispersion.png
    dispersion_plot_pdf: File  # <prefix>.deseq2.dispersion.pdf
    rdata: File  # <prefix>.dds.rld.rds
    size_factors: File  # <prefix>.deseq2.sizefactors.tsv
    normalised_counts: File  # <prefix>.normalised_counts.tsv
    rlog_counts: File | None  # <prefix>.rlog.tsv
    vst_counts: File | None  # <prefix>.vst.tsv
    model: File  # <prefix>.deseq2.model.txt
    session_info: File  # <prefix>.R_sessionInfo.log
    versions: File  # versions.yml


async def deseq2_differential(
    samplesheet: File,
    counts: File,
    prefix: str,
    contrast_variable: str | None = None,
    reference: str | None = None,
    target: str | None = None,
    formula: str | None = None,
    comparison: str | None = None,
    control_genes: File | None = None,
    transcript_lengths: File | None = None,
    args: str = "",
) -> Deseq2DifferentialResult:
    """DESeq2 results for one contrast of ``counts`` (samples described by ``samplesheet``)."""
    results = await deseq2_differential_cmd(
        samplesheet=samplesheet,
        counts=counts,
        control_genes=control_genes,
        transcript_lengths=transcript_lengths,
        prefix=prefix,
        contrast_variable=contrast_variable or "null",
        reference=reference or "null",
        target=target or "null",
        formula=formula or "null",
        comparison=comparison or "null",
        args=args or "null",
    )

    async def pick(name: str) -> File:
        found = await results.get_file(name)
        if found is None:
            raise FileNotFoundError(f"DESeq2 wrote no {name}")
        return found

    return Deseq2DifferentialResult(
        results=await pick(f"{prefix}.deseq2.results.tsv"),
        dispersion_plot_png=await pick(f"{prefix}.deseq2.dispersion.png"),
        dispersion_plot_pdf=await pick(f"{prefix}.deseq2.dispersion.pdf"),
        rdata=await pick(f"{prefix}.dds.rld.rds"),
        size_factors=await pick(f"{prefix}.deseq2.sizefactors.tsv"),
        normalised_counts=await pick(f"{prefix}.normalised_counts.tsv"),
        rlog_counts=await results.get_file(f"{prefix}.rlog.tsv"),
        vst_counts=await results.get_file(f"{prefix}.vst.tsv"),
        model=await pick(f"{prefix}.deseq2.model.txt"),
        session_info=await pick(f"{prefix}.R_sessionInfo.log"),
        versions=await pick("versions.yml"),
    )
