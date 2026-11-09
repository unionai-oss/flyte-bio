"""hasheddrops — demultiplex cell-hashing data with DropletUtils' hashedDrops.

Exposes :func:`hasheddrops`, which reads a 10x HTO count matrix directory and
calls each droplet's hashtag with ``hashedDrops``. With ``run_empty_drops``,
cells are first called from the RNA matrix with ``emptyDrops`` and only those
are demultiplexed, against the ambient HTO profile. Outputs keep upstream's
names:

- ``<prefix>_results_hasheddrops.csv`` and ``<prefix>_hasheddrops.rds``: the
  hashedDrops calls;
- ``<prefix>_id_to_hash.csv``: hashtag index to name;
- ``<prefix>_plot_hasheddrops.png``: log-fold changes of the calls;
- ``<prefix>_emptyDrops.{csv,rds,png}``: the emptyDrops results (empty when
  it isn't run);
- ``<prefix>_params_hasheddrops.csv``: the parameters used;
- ``versions.yml``.

The R code is upstream's template, unchanged but for taking the templated
values (matrix paths, ``runEmptyDrops``, prefix and ``args``) as command-line
arguments. ``args`` are upstream's ``ext.args``: ``--<option> <value>`` pairs
overriding the script's defaults (e.g. ``--niters 1000``). As upstream, the
matrices are linked into the working directory under ``hto_matrix_name`` and
``rna_matrix_name`` (upstream's staged names, ``hto`` and ``rna`` for the test
data); the parameters table records the HTO one.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

HASHEDDROPS_IMAGE = "community.wave.seqera.io/library/bioconductor-dropletutils_r-seurat:e1dff3a0fb7c5920"
# Upstream: process_low.
DEFAULT_RESOURCES = flyte.Resources(cpu=2, memory="8Gi")

HASHEDDROPS_R = r"""
#!/usr/bin/env Rscript

cli_args <- commandArgs(trailingOnly = TRUE)

################################################
################################################
## Functions                                  ##
################################################
################################################

# Helper function for NULL condition
string_to_null <- function(x, val = "NULL") if (x == val) NULL else x
null_to_string <- function(x, val = "NULL") if (is.null(x)) val else x

string_to_logical <- function(input) {
    if (input == "FALSE") {
        FALSE
    } else if (input == "TRUE") {
        TRUE
    } else {
        stop(paste0(input, " is not a valid logical. Use 'FALSE' or 'TRUE'."))
    }
}

#' Check for Non-Empty, Non-Whitespace String
#'
#' This function checks if the input is non-NULL and contains more than just whitespace.
#' It returns TRUE if the input is a non-empty, non-whitespace string, and FALSE otherwise.
#'
#' @param input A variable to check.
#' @return A logical value: TRUE if the input is a valid, non-empty, non-whitespace string; FALSE otherwise.

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

################################################
################################################
## PARSE PARAMETERS FROM NEXTFLOW             ##
################################################
################################################

# Set defaults and classes
opt <- list(
    # File inputs
    hto_matrix = cli_args[1],
    rna_matrix = cli_args[2],
    runEmptyDrops = string_to_logical(cli_args[3]),

    # emptyDrops Parameters
    lower = 100,        # A numeric scalar specifying the lower bound on the total UMI count, at or below which all barcodes are assumed to correspond to empty droplets.
    niters = 10000,     # An integer scalar specifying the number of iterations to use for the Monte Carlo p-value calculations.
    testAmbient = TRUE, # A logical scalar indicating whether results should be returned for barcodes with totals less than or equal to lower.
    ignore = NULL,      # A numeric scalar specifying the lower bound on the total UMI count, at or below which barcodes will be ignored.
    alpha = Inf,        # A numeric scalar specifying the scaling parameter for the Dirichlet-multinomial sampling scheme.
    round = TRUE,       # Logical scalar indicating whether to check for non-integer values in m and, if present, round them for ambient profile estimation.
    byRank = NULL,      # An integer scalar parametrizing an alternative method for identifying assumed empty droplets. If set, this is used to redefine lower and any specified value for lower is ignored.
    isCellFDR = 0.01,   # Threshold to filter the cells.

    # hashedDrops Parameters
    ambient = TRUE,          # Whether to use the relative abundance of each HTO in the ambient solution from emptyDrops, set TRUE only when test_ambient is TRUE.
    minProp = 0.05,          # Numeric scalar to be used to infer the ambient profile when ambient=NULL.
    pseudoCount = 5,         # A numeric scalar specifying the minimum pseudo-count when computing logfold changes.
    constantAmbient = FALSE, # Logical scalar indicating whether a constant level of ambient contamination should be used to estimate LogFC2 for all cells.
    doubletNmads = 3,        # A numeric scalar specifying the number of median absolute deviations (MADs) to use to identify doublets.
    doubletMin = 2,          # A numeric scalar specifying the minimum threshold on the log-fold change to use to identify doublets.
    doubletMixture = FALSE,  # Logical scalar indicating whether to use a 2-component mixture model to identify doublets.
    confidentNmads = 3,      # A numeric scalar specifying the number of MADs to use to identify confidently assigned singlets.
    confidentMin = 2,        # A numeric scalar specifying the minimum threshold on the log-fold change to use to identify singlets.
    combinations = NULL,     # An integer matrix specifying valid combinations of HTOs. Each row corresponds to a single sample and specifies the indices of rows in x corresponding to the HTOs used to label that sample.

    # others
    gene_col = 2,                                                                 # Specify which column of genes.tsv or features.tsv to use for gene names; default is 2.
    prefix = cli_args[4] # Prefix name for output files.
)
opt_types <- lapply(opt, class)

# Apply parameter overrides
args_string <- cli_args[5]
args_opt <- if (is_valid_string(args_string)) parse_args(args_string) else list()
for ( ao in names(args_opt)){
    if (! ao %in% names(opt)){
        stop(paste("Invalid option:", ao))
    }else{
        # Handle special cases for NULL values and logicals
        if (args_opt[[ao]] == "NULL") {
            opt[[ao]] <- NULL
        } else if (opt_types[[ao]] == "logical") {
            opt[[ao]] <- string_to_logical(args_opt[[ao]])
        } else if (! is.null(opt[[ao]])){
            # Preserve classes from defaults where possible
            opt[[ao]] <- as(args_opt[[ao]], opt_types[[ao]])
        } else {
            opt[[ao]] <- args_opt[[ao]]
        }
    }
}

# Set individual variables for backward compatibility and cleaner code
hto_matrix <- opt$hto_matrix
rna_matrix <- opt$rna_matrix
runEmptyDrops <- opt$runEmptyDrops
lower <- opt$lower
niters <- opt$niters
testAmbient <- opt$testAmbient
round <- opt$round
byRank <- opt$byRank
isCellFDR <- opt$isCellFDR
gene_col <- opt$gene_col
ignore <- opt$ignore
alpha <- opt$alpha
ambient <- opt$ambient
minProp <- opt$minProp
pseudoCount <- opt$pseudoCount
constantAmbient <- opt$constantAmbient
doubletNmads <- opt$doubletNmads
doubletMin <- opt$doubletMin
doubletMixture <- opt$doubletMixture
confidentNmads <- opt$confidentNmads
confidentMin <- opt$confidentMin
combinations <- opt$combinations
prefix <- opt$prefix

# check if the file exists
if (! file.exists(hto_matrix)){
    stop(paste0(hto_matrix, ' is not a valid file'))
}

################################################
################################################
## Finish loading libraries                   ##
################################################
################################################

library(Seurat)  # for Read10X()
library(DropletUtils) # for hashedDrops() and emptyDrops()

################################################
################################################
## Main Process                               ##
################################################
################################################

hto <- Read10X(data.dir = hto_matrix, gene.column = gene_col)

# determine hto_input and ambient_input
if (runEmptyDrops) {

    rna <- Read10X(data.dir = rna_matrix, gene.column = gene_col)

    emptyDrops_out <- emptyDrops(
        rna,
        lower = lower,
        niters = niters,
        test.ambient = testAmbient,
        ignore = NULL,
        alpha = alpha,
        round = round,
        by.rank = byRank
    )

    # which droplets are actual cells
    is.cell <- emptyDrops_out$FDR <= isCellFDR & !is.na(emptyDrops_out$FDR)
    hto_input <- hto[, which(is.cell)]

    if (ambient) {
        ambient_input <- metadata(emptyDrops_out)$ambient
    } else {
        ambient_input <- NULL
    }
} else {
    ambient_input <- NULL
    hto_input <- hto

    # only important for saving the results
    emptyDrops_out <- data.frame()
}

hashedDrops_out <- hashedDrops(
    hto_input,
    min.prop = minProp,
    ambient = ambient_input,
    pseudo.count = pseudoCount,
    constant.ambient = constantAmbient,
    doublet.nmads = doubletNmads,
    doublet.min = doubletMin,
    doublet.mixture = doubletMixture,
    confident.nmads = confidentNmads,
    confident.min = confidentMin,
    combinations = combinations
)

################################################
################################################
## SAVING RESULTS                             ##
################################################
################################################

#----- saving parameters in a dataframe  ------#
Argument <- c(
    "hto_matrix",
    "lower",
    "niters",
    "testAmbient",
    "ignore",
    "alpha",
    "round",
    "byRank",
    "isCellFDR",
    "gene_col",
    "ambient",
    "minProp",
    "pseudoCount",
    "constantAmbient",
    "doubletNmads",
    "doubletMin",
    "doubletMixture",
    "confidentNmads",
    "confidentMin",
    "combinations"
)

Value <- c(
    hto_matrix,
    lower,
    niters,
    testAmbient,
    null_to_string(ignore),
    null_to_string(alpha),
    round,
    null_to_string(byRank),
    isCellFDR,
    gene_col,
    ambient,
    null_to_string(minProp),
    pseudoCount,
    constantAmbient,
    doubletNmads,
    doubletMin,
    doubletMixture,
    confidentNmads,
    confidentMin,
    null_to_string(combinations)
)

params <- data.frame(Argument, Value)
write.csv(params, paste0(prefix ,"_params_hasheddrops.csv"))

#--------- save emptyDrops() results  ---------#

# create a plot with results from emptyDrops() or save an empty png
png(paste0(prefix, "_emptyDrops.png"))
if(runEmptyDrops){
    colors <- ifelse(is.cell, "red", "black")
    plot(emptyDrops_out$Total, -emptyDrops_out$LogProb, col = colors, xlab = "Total UMI count", ylab = "-Log Probability")
}else{
    plot.new()
}
dev.off()

write.csv(emptyDrops_out,paste0(prefix, "_emptyDrops.csv"))
saveRDS(emptyDrops_out,file = paste0(prefix, "_emptyDrops.rds"))

#--------- save hashedDrops() results ---------#


write.csv(params, paste0(prefix, "_params_hasheddrops.csv"))
write.csv(hashedDrops_out,paste0(prefix,"_results_hasheddrops.csv"))
saveRDS(hashedDrops_out,file = paste0(prefix,"_hasheddrops.rds"))

png(paste0(prefix, "_plot_hasheddrops.png"))
if (sum(is.na(hashedDrops_out$LogFC2)) != length(hashedDrops_out$LogFC2)) {

    colors <- ifelse(hashedDrops_out$Confident,
    "black",
    ifelse(hashedDrops_out$Doublet, "red", "grey")
    )

    plot(
    hashedDrops_out$LogFC,
    hashedDrops_out$LogFC2,
    col = colors,
    xlab = "Log-fold change between best and second HTO",
    ylab = "Log-fold change between second HTO and ambient"
    )
}else{

    plot.new()
}

# mapping from integers (e.g., in Best) to HTO names or combinations.
# If combinations are specified, the index will map to the joined HTO names separated by a "+".
# Otherwise, it will simply use the row name.
hto_names <- rownames(hto)
if (!is.null(combinations)){
    hto_names <- apply(combinations, 1, function(row) paste(row, collapse = "+"))
    # In some applications, samples are labelled with a combination of HTOs to enable achieve greater
    # multiplexing throughput. This is accommodated by passing combinations to specify the valid
    # HTO combinations that were used for sample labelling. Each row of combinations corresponds
    # to a sample and should contain non-duplicated row indices of x corresponding to the HTOs used in
    # that sample.
    # Source: https://bioconductor.statistik.tu-dortmund.de/packages/3.18/bioc/manuals/DropletUtils/man/DropletUtils.pdf

    # If combinations is specified, Best instead specifies the sample (i.e., row index of combinations).
    # Source: https://rdrr.io/github/MarioniLab/DropletUtils/man/hashedDrops.html
}

# Create a data frame mapping names to indices
hto_map <- data.frame(
Index = seq_along(hto_names),
HTO = hto_names
)

# Write to CSV
write.csv(hto_map, file = paste0(prefix,"_id_to_hash.csv"), row.names = FALSE)


dev.off()


################################################
################################################
## VERSIONS FILE                              ##
################################################
################################################

r.version <- paste(R.version[['major']],R.version[['minor']], sep = ".")
seurat.version <- as.character(packageVersion('Seurat'))
dropletutils.version <- as.character(packageVersion('DropletUtils'))

writeLines(
    c(
        '"HASHEDDROPS":',
        paste('    r-base:', r.version),
        paste('    r-seurat:', seurat.version),
        paste('    dropletutils:', dropletutils.version)
    ),
'versions.yml')

################################################
################################################
################################################
################################################
""".strip()  # noqa: E501

OUTPUT_SUFFIXES = (
    "_emptyDrops.png",
    "_emptyDrops.csv",
    "_emptyDrops.rds",
    "_results_hasheddrops.csv",
    "_id_to_hash.csv",
    "_hasheddrops.rds",
    "_plot_hasheddrops.png",
    "_params_hasheddrops.csv",
)

hasheddrops_cmd = shell.create(
    name="hasheddrops",
    image=HASHEDDROPS_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "hto_matrix": Dir,
        "hto_matrix_name": str,
        "rna_matrix": Dir | None,
        "rna_matrix_name": str,
        "run_empty_drops": str,
        "prefix": str,
        "args": str,
    },
    defaults={"hto_matrix_name": "hto", "rna_matrix_name": "rna", "args": "null"},
    outputs={"results": Dir},
    script=r"""
        HTO={inputs.hto_matrix}
        HTO_NAME={inputs.hto_matrix_name}
        RNA={inputs.rna_matrix}
        RNA_NAME={inputs.rna_matrix_name}
        P={inputs.prefix}
        OUT={outputs.results}
        W=$(mktemp -d); cd "$W"
        ln -s "$HTO" "$HTO_NAME"
        if [ -n "$RNA" ] && [ -d "$RNA" ]; then ln -s "$RNA" "$RNA_NAME"; else RNA_NAME=""; fi
        cat > hasheddrops.R <<'END_OF_SCRIPT'
HASHEDDROPS_R
END_OF_SCRIPT
        Rscript hasheddrops.R "$HTO_NAME" "$RNA_NAME" {inputs.run_empty_drops} "$P" {inputs.args} || exit 1
        for S in OUTPUT_SUFFIXES; do mv "$P$S" "$OUT"/; done
        mv versions.yml "$OUT"/
    """.replace("HASHEDDROPS_R", HASHEDDROPS_R).replace("OUTPUT_SUFFIXES", " ".join(OUTPUT_SUFFIXES)),
)

env = flyte.TaskEnvironment.from_task("hasheddrops", hasheddrops_cmd.as_task())


@dataclass
class HashedDropsResult:
    results: File  # <prefix>_results_hasheddrops.csv
    rds: File  # <prefix>_hasheddrops.rds
    id_to_hash: File  # <prefix>_id_to_hash.csv
    plot: File  # <prefix>_plot_hasheddrops.png
    empty_drops_csv: File  # <prefix>_emptyDrops.csv
    empty_drops_rds: File  # <prefix>_emptyDrops.rds
    empty_drops_plot: File  # <prefix>_emptyDrops.png
    params: File  # <prefix>_params_hasheddrops.csv
    versions: File  # versions.yml


async def hasheddrops(
    hto_matrix: Dir,
    prefix: str,
    rna_matrix: Dir | None = None,
    run_empty_drops: bool = False,
    hto_matrix_name: str = "hto",
    rna_matrix_name: str = "rna",
    args: str = "",
) -> HashedDropsResult:
    """hashedDrops calls for a 10x HTO matrix directory, optionally on the cells emptyDrops finds in ``rna_matrix``."""
    if run_empty_drops and rna_matrix is None:
        raise ValueError("hasheddrops: run_empty_drops needs rna_matrix")
    results = await hasheddrops_cmd(
        hto_matrix=hto_matrix,
        hto_matrix_name=hto_matrix_name,
        rna_matrix=rna_matrix,
        rna_matrix_name=rna_matrix_name,
        run_empty_drops="TRUE" if run_empty_drops else "FALSE",
        prefix=prefix,
        # As upstream's template renders an unset ext.args.
        args=args or "null",
    )

    async def pick(name: str) -> File:
        found = await results.get_file(name)
        if found is None:
            raise FileNotFoundError(f"hasheddrops wrote no {name}")
        return found

    return HashedDropsResult(
        results=await pick(f"{prefix}_results_hasheddrops.csv"),
        rds=await pick(f"{prefix}_hasheddrops.rds"),
        id_to_hash=await pick(f"{prefix}_id_to_hash.csv"),
        plot=await pick(f"{prefix}_plot_hasheddrops.png"),
        empty_drops_csv=await pick(f"{prefix}_emptyDrops.csv"),
        empty_drops_rds=await pick(f"{prefix}_emptyDrops.rds"),
        empty_drops_plot=await pick(f"{prefix}_emptyDrops.png"),
        params=await pick(f"{prefix}_params_hasheddrops.csv"),
        versions=await pick("versions.yml"),
    )
