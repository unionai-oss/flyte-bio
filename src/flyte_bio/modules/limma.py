"""limma — differential abundance of an intensity or count matrix with limma.

Exposes :func:`limma_differential`, which tests one contrast: either
``contrast_variable`` with its ``reference`` and ``target`` levels, or a
model ``formula`` with a ``comparison`` (a coefficient name or contrast
expression). For counts, ``--use_voom TRUE`` in ``args`` applies voom
(writing the normalised log-CPM matrix); ``--block`` fits a mixed model on
that variable. Outputs keep upstream's names:

- ``<prefix>.limma.results.tsv``, ``<prefix>.limma.model.txt``;
- ``<prefix>.limma.mean_difference.png``, ``<prefix>.MArrayLM.limma.rds``;
- ``<prefix>.normalised_counts.tsv`` (with voom);
- ``<prefix>.R_sessionInfo.log`` and ``versions.yml``.

``args`` are upstream's ``ext.args`` (``--<option> <value>`` pairs, e.g.
``--sample_id_col sample --probe_id_col gene_id``). The R code is upstream's
template, unchanged but for taking the templated values as command-line
arguments; the prefix stands for both the templated ``ext.prefix`` and
sample ID. As upstream renders unset values, a missing contrast field is
passed as ``null``.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

LIMMA_IMAGE = "community.wave.seqera.io/library/bioconductor-edger_bioconductor-limma:176c202c82450990"
# Upstream: process_single.
DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory="8Gi")

LIMMA_DE_R = r"""
#!/usr/bin/env Rscript

cli_args <- commandArgs(trailingOnly = TRUE)

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
    parsed_args[ ( ! parsed_args %in%  c('', 'null')) & ! is.na(parsed_args)]
}

#' Flexibly read CSV or TSV files
#'
#' @param file Input file
#' @param header Passed to read.delim()
#' @param row.names Passed to read.delim()
#'
#' @return output Data frame

read_delim_flexible <- function(file, header = TRUE, row.names = NULL, check.names = FALSE){

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
    formula                 = cli_args[7],
    contrast_string         = cli_args[8],
    blocking_variables = NULL,
    probe_id_col = "probe_id",
    sample_id_col = "experiment_accession",
    subset_to_contrast_samples = FALSE,
    exclude_samples_col = NULL,
    exclude_samples_values = NULL,
    use_voom = FALSE,
    number = Inf,
    ndups = NULL,                # lmFit
    spacing = NULL,              # lmFit
    block = NULL,                # lmFit
    correlation = NULL,          # lmFit
    method = 'ls',               # lmFit
    proportion = 0.01,           # eBayes
    stdev.coef.lim = '0.1,4',    # eBayes
    trend = FALSE,               # eBayes
    robust = FALSE,              # eBayes
    winsor.tail.p = '0.05,0.1',  # eBayes
    adjust.method = "BH",        # topTable
    p.value = 1,                 # topTable
    lfc = 0,                     # topTable
    confint = FALSE ,            # topTable
    round_digits = NULL,
    seed = NULL
)
opt_types <- lapply(opt, class)

# Apply parameter overrides

args_opt <- parse_args(cli_args[9])
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
if ( ! is.null(opt$seed)){
    opt$seed <- as.numeric(opt$seed)
    set.seed(opt$seed)
}

# If there is no option supplied, convert string "null" to NULL
keys <- c("formula", "contrast_string", "contrast_variable", "reference_level", "target_level")
opt[keys] <- lapply(opt[keys], nullify)

# Check if required parameters have been provided

required_opts <- c('output_prefix')
missing <- required_opts[unlist(lapply(opt[required_opts], is.null)) | ! required_opts %in% names(opt)]

if (length(missing) > 0){
    stop(paste("Missing required options:", paste(missing, collapse=', ')))
}

# Check file inputs are valid

for (file_input in c('count_file', 'sample_file')){
    if (is.null(opt[[file_input]])) {
        stop(paste("Please provide", file_input), call. = FALSE)
    }

    if (! file.exists(opt[[file_input]])){
        stop(paste0('Value of ', file_input, ': ', opt[[file_input]], ' is not a valid file'))
    }
}

# Convert things to vectors where needed
vector_opt <- c('stdev.coef.lim', 'winsor.tail.p')
opt[vector_opt] = lapply(strsplit(unlist(opt[vector_opt]), ','), as.numeric)

################################################
################################################
## Finish loading libraries                   ##
################################################
################################################

library(limma)
library(edgeR)

################################################
################################################
## READ IN COUNTS FILE AND SAMPLE METADATA    ##
################################################
################################################

intensities.table <-
    read_delim_flexible(
        file = opt$count_file,
        header = TRUE,
        row.names = opt$probe_id_col
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

# Check that all samples specified in the input sheet are present in the
# intensities table. Assuming they are, subset and sort the count table to
# match the sample sheet

missing_samples <-
    sample.sheet[!rownames(sample.sheet) %in% colnames(intensities.table), opt$sample_id_col]

if (length(missing_samples) > 0) {
    stop(paste(
        length(missing_samples),
        'specified samples missing from count table:',
        paste(missing_samples, collapse = ',')
    ))
} else{
    # Save any non-count data, will gene metadata etc we might need later
    nonintensities.table <-
        intensities.table[, !colnames(intensities.table) %in% rownames(sample.sheet), drop = FALSE]
    intensities.table <- intensities.table[, rownames(sample.sheet)]
}

################################################
################################################
## CHECK CONTRAST SPECIFICATION               ##
################################################
################################################

contrast_variable <- make.names(opt$contrast_variable)
blocking.vars <- c()

if (!is.null(opt$contrast_variable)) {
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
            'Please choose reference and target levels that are present in the',
            contrast_variable,
            'column of the sample sheet'
            )
        )
    } else if (!is.null(opt$blocking_variables)) {
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

    # Handle conflicts between blocking variables and block
    if (!is.null(opt$block) && !is.null(opt$blocking_variables)) {
        if (opt$block %in% blocking.vars) {
            warning(paste("Variable", opt$block, "is specified both as a fixed effect and a random effect. It will be treated as a random effect only."))
            blocking.vars <- setdiff(blocking.vars, opt$block)
            if (length(blocking.vars) == 0) {
                opt$blocking_variables <- NULL
            } else {
                opt$blocking_variables <- paste(blocking.vars, collapse = ';')
            }
        }
    }
}
# Optionally, subset to only the samples involved in the contrast

if (opt$subset_to_contrast_samples){
    sample_selector <- sample.sheet[[contrast_variable]] %in% c(opt$target_level, opt$reference_level)
    selected_samples <- sample.sheet[sample_selector, opt$sample_id_col]
    intensities.table <- intensities.table[, selected_samples]
    sample.sheet <- sample.sheet[selected_samples, ]
}

# Optionally, remove samples with specified values in a given field (probably
# don't use this as well as the above)

if ((! is.null(opt$exclude_samples_col)) && (! is.null(opt$exclude_samples_values))){
    exclude_values = unlist(strsplit(opt$exclude_samples_values, split = ';'))

    if (! opt$exclude_samples_col %in% colnames(sample.sheet)){
        stop(paste(opt$exclude_samples_col, ' specified to subset samples is not a valid sample sheet column'))
    }

    print(paste0('Excluding samples with values of ', opt$exclude_samples_values, ' in ', opt$exclude_samples_col))
    sample_selector <- ! sample.sheet[[opt$exclude_samples_col]] %in% exclude_values

    selected_samples <- sample.sheet[sample_selector, opt$sample_id_col]
    intensities.table <- intensities.table[, selected_samples]
    sample.sheet <- sample.sheet[selected_samples, ]
}

################################################
################################################
## Build the Model Formula and Run Limma      ##
################################################
################################################

if (!is.null(opt$formula)) {
    model <- opt$formula
    model_formula <- as.formula(model)
    cat("Using user-specified formula:\n   ", deparse(model_formula), "\n")
    design <- model.matrix(model_formula, data = sample.sheet)
    colnames(design) <- make.names(colnames(design))
    cat("Column names after make.names():\n   ", paste(colnames(design), collapse = ", "), "\n")

} else {
    # Put the contrast variable first in zero-intercept designs so each
    # contrast level is represented directly in the design matrix.
    model_vars <- contrast_variable

    if (!is.null(opt$blocking_variables)) {
        # Include blocking variables (including pairing variables if any)
        model_vars <- c(model_vars, blocking.vars)
    }

    # Construct the model formula
    model <- paste('~ 0 +', paste(model_vars, collapse = '+'))

    # Make sure all the appropriate variables are factors
    vars_to_factor <- model_vars  # All variables in the model need to be factors
    for (v in vars_to_factor) {
        sample.sheet[[v]] <- as.factor(sample.sheet[[v]])
    }

    # Generate the design matrix
    design <- model.matrix(
        as.formula(model),
        data=sample.sheet
    )

    # Adjust column names for the contrast variable
    colnames(design) <- sub(
        paste0('^', contrast_variable),
        paste0(contrast_variable, '.'),
        colnames(design)
    )

    # Adjust column names to be syntactically valid
    colnames(design) <- make.names(colnames(design))
    cat("Final column names after make.names():\n")
    print(colnames(design))
}

# Perform voom normalisation for RNA-seq data
if (!is.null(opt$use_voom) && opt$use_voom) {
    # Create a DGEList object for RNA-seq data
    dge <- DGEList(counts = intensities.table)

    # Normalize counts using TMM
    dge <- calcNormFactors(dge, method = "TMM")

    # Run voom to transform the data
    data_for_fit <- voom(dge, design)
} else {
    # Use as.matrix for regular microarray analysis
    data_for_fit <- as.matrix(intensities.table)
}

if (!is.null(opt$block)) {
    corfit = duplicateCorrelation(data_for_fit, design = design, block = sample.sheet[[opt$block]])
    if (!is.null(opt$use_voom) && opt$use_voom) {
        data_for_fit <- voom(counts = dge, design = design, plot = FALSE, correlation = corfit$consensus.correlation)
    }
}

# For Voom, write the normalized counts matrix to a TSV file
if (!is.null(opt$use_voom) && opt$use_voom) {
    normalized_counts <- data_for_fit$E
    if (! is.null(opt$round_digits)){
        normalized_counts <- apply(normalized_counts, 2, function(x) round(x, opt$round_digits))
    }
    normalized_counts_with_genes <- data.frame(Gene = rownames(normalized_counts), normalized_counts, check.names = FALSE, row.names = NULL)
    colnames(normalized_counts_with_genes)[1] <- opt$probe_id_col
    write.table(normalized_counts_with_genes,
        file = paste(opt$output_prefix, "normalised_counts.tsv", sep = '.'),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE)
}

# Prepare for and run lmFit()

lmfit_args <- list(
    object = data_for_fit,
    design = design
)

# Include optional parameters if provided
if (! is.null(opt$ndups)){
    lmfit_args[['ndups']] <- as.numeric(opt$ndups)
}
if (! is.null(opt$spacing)){
    lmfit_args[['spacing']] <- as.numeric(opt$spacing)
}
if (! is.null(opt$block)){
    lmfit_args[['block']] <- sample.sheet[[opt$block]]
}
if (! is.null(opt$correlation)){
    lmfit_args[['correlation']] <- as.numeric(opt$correlation)
} else if (! is.null(opt$block)){
    lmfit_args[['correlation']] <- corfit$consensus.correlation
}
if (! is.null(opt$method)){
    lmfit_args[['method']] <- opt$method
}

fit <- do.call(lmFit, lmfit_args)

# Contrasts bit

# Create the contrast string for the specified comparison
if (!is.null(opt$contrast_string)) {
    cat("Using contrast string:", opt$contrast_string, "\n")
    contrast_string <- as.character(opt$contrast_string)
    contrast.matrix <- makeContrasts(contrasts=contrast_string, levels=colnames(design))

} else {

    # Construct the expected column names for the target and reference levels in the design matrix
    treatment_target <- make.names(paste0(contrast_variable, ".", opt$target_level))
    treatment_reference <- make.names(paste0(contrast_variable, ".", opt$reference_level))

    # Determine how to construct the contrast string based on which levels are present in the design matrix
    if ((treatment_target %in% colnames(design)) && (treatment_reference %in% colnames(design))) {
        # Both target and reference levels are present in the design matrix
        # We can directly compare the two levels
        contrast_string <- paste0(treatment_target, "-", treatment_reference)
    } else if (treatment_target %in% colnames(design)) {
        # Only the target level is present in the design matrix
        # The reference level may have been omitted due to collinearity or being set as the baseline
        # We compare the target level to zero (implicit reference)
        contrast_string <- paste0(treatment_target, "- 0")
    } else if (treatment_reference %in% colnames(design)) {
        # Only the reference level is present in the design matrix
        # The target level may have been omitted from the design matrix
        # We compare zero (implicit target) to the reference level
        contrast_string <- paste0("0 - ", treatment_reference)
    } else {
        # Neither level is present in the design matrix
        # This indicates an error; the specified levels are not found
        stop(paste0(treatment_target, " and ", treatment_reference, " not found in design matrix"))
    }
    cat("Using contrast string:", contrast_string, "
")
    contrast.matrix <- makeContrasts(contrasts=contrast_string, levels=design)

}

# Create the contrast matrix
fit2 <- contrasts.fit(fit, contrast.matrix)

# Prepare for and run eBayes

ebayes_args <- list(
    fit = fit2
)

if (! is.null(opt$proportion)){
    ebayes_args[['proportion']] <- as.numeric(opt$proportion)
}
if (! is.null(opt$stdev.coef.lim)){
    ebayes_args[['stdev.coef.lim']] <- as.numeric(opt$stdev.coef.lim)
}
if (! is.null(opt$trend)){
    ebayes_args[['trend']] <- as.logical(opt$trend)
}
if (! is.null(opt$robust)){
    ebayes_args[['robust']] <- as.logical(opt$robust)
}
if (! is.null(opt$winsor.tail.p)){
    ebayes_args[['winsor.tail.p']] <- as.numeric(opt$winsor.tail.p)
}

fit2 <- do.call(eBayes, ebayes_args)

# Run topTable() to generate a results data frame
toptable_args <- list(
    fit = fit2,
    sort.by = 'none',
    number = nrow(intensities.table)
)

if (! is.null(opt$adjust.method)){
    toptable_args[['adjust.method']] <- opt$adjust.method
}
if (! is.null(opt$p.value)){
    toptable_args[['p.value']] <- as.numeric(opt$p.value)
}
if (! is.null(opt$lfc)){
    toptable_args[['lfc']] <- as.numeric(opt$lfc)
}
if (! is.null(opt$confint)){
    toptable_args[['confint']] <- as.logical(opt$confint)
}

comp.results <- do.call(topTable, toptable_args)[rownames(intensities.table),]

################################################
################################################
## Generate outputs                           ##
################################################
################################################

contrast.name <-
    paste(opt$target_level, opt$reference_level, sep = "_vs_")
cat("Saving results for ", contrast.name, " ...\n", sep = "")

# Differential expression table - note very limited rounding for consistency of
# results

if (! is.null(opt$round_digits)){
    comp.results <- apply(data.frame(comp.results), 2, function(x) round(x, opt$round_digits))
}
out_df <- cbind(
    setNames(data.frame(rownames(comp.results)), opt$probe_id_col),
    data.frame(comp.results[, !(colnames(comp.results) %in% opt$probe_id_col)], check.names = FALSE)
)
write.table(
    out_df,
    file = paste(opt$output_prefix, 'limma.results.tsv', sep = '.'),
    col.names = TRUE,
    row.names = FALSE,
    sep = '\t',
    quote = FALSE
)

# Dispersion plot

png(
    file = paste(opt$output_prefix, 'limma.mean_difference.png', sep = '.'),
    width = 600,
    height = 600
)
plotMD(fit2)
dev.off()

# R object for other processes to use

saveRDS(fit2, file = paste(opt$output_prefix, 'MArrayLM.limma.rds', sep = '.'))

# Save model to file

write(model, file=paste(opt$output_prefix, 'limma.model.txt', sep = '.'))

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

r.version <- strsplit(version[['version.string']], ' ')[[1]][3]
limma.version <- as.character(packageVersion('limma'))

writeLines(
    c(
        '"LIMMA_DIFFERENTIAL":',
        paste('    r-base:', r.version),
        paste('    bioconductor-limma:', limma.version)
    ),
'versions.yml')

################################################
################################################
################################################
################################################
""".strip()  # noqa: E501

limma_differential_cmd = shell.create(
    name="limma_differential",
    image=LIMMA_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "samplesheet": File,
        "intensities": File,
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
        INTENSITIES=({inputs.intensities})
        P={inputs.prefix}
        OUT={outputs.results}
        W=$(mktemp -d); cd "$W"
        SHEET_NAME=$(basename "${SHEET[0]}"); ln -s "${SHEET[0]}" "$SHEET_NAME"
        INTENSITIES_NAME=$(basename "${INTENSITIES[0]}"); ln -s "${INTENSITIES[0]}" "$INTENSITIES_NAME"
        cat > limma_de.R <<'END_OF_SCRIPT'
LIMMA_DE_R
END_OF_SCRIPT
        Rscript limma_de.R "$P" "$INTENSITIES_NAME" "$SHEET_NAME" {inputs.contrast_variable} {inputs.reference} \
            {inputs.target} {inputs.formula} {inputs.comparison} {inputs.args} || exit 1
        mv "$P".limma.* "$P".MArrayLM.limma.rds "$P".R_sessionInfo.log versions.yml "$OUT"/
        if [ -e "$P".normalised_counts.tsv ]; then mv "$P".normalised_counts.tsv "$OUT"/; fi
    """.replace("LIMMA_DE_R", LIMMA_DE_R),
)

env = flyte.TaskEnvironment.from_task("limma", limma_differential_cmd.as_task())


@dataclass
class LimmaDifferentialResult:
    results: File  # <prefix>.limma.results.tsv
    md_plot: File  # <prefix>.limma.mean_difference.png
    rdata: File  # <prefix>.MArrayLM.limma.rds
    model: File  # <prefix>.limma.model.txt
    normalised_counts: File | None  # <prefix>.normalised_counts.tsv, with voom
    session_info: File  # <prefix>.R_sessionInfo.log
    versions: File  # versions.yml


async def limma_differential(
    samplesheet: File,
    intensities: File,
    prefix: str,
    contrast_variable: str | None = None,
    reference: str | None = None,
    target: str | None = None,
    formula: str | None = None,
    comparison: str | None = None,
    args: str = "",
) -> LimmaDifferentialResult:
    """limma results for one contrast of ``intensities`` (samples described by ``samplesheet``)."""
    results = await limma_differential_cmd(
        samplesheet=samplesheet,
        intensities=intensities,
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
            raise FileNotFoundError(f"limma wrote no {name}")
        return found

    return LimmaDifferentialResult(
        results=await pick(f"{prefix}.limma.results.tsv"),
        md_plot=await pick(f"{prefix}.limma.mean_difference.png"),
        rdata=await pick(f"{prefix}.MArrayLM.limma.rds"),
        model=await pick(f"{prefix}.limma.model.txt"),
        normalised_counts=await results.get_file(f"{prefix}.normalised_counts.tsv"),
        session_info=await pick(f"{prefix}.R_sessionInfo.log"),
        versions=await pick("versions.yml"),
    )
