"""bff — demultiplex cell-hashing data with BFF (cellhashR's bimodal flexible fitting).

Exposes :func:`bff`, which reads a 10x HTO count matrix directory and calls
each cell's hashtag with cellhashR's ``bff_raw`` and/or ``bff_cluster``
callers (``methods``: ``RAW``, ``CLUSTER`` or ``COMBINED``), optionally after
cellhashR's count-matrix preprocessing. Outputs keep upstream's names:

- ``<prefix>_assignment_bff.csv``: the per-cell calls and consensus;
- ``<prefix>_metrics_bff.csv``: cellhashR's metrics;
- ``<prefix>_params_bff.csv``: the parameters used;
- ``versions.yml``.

The R code is upstream's template, unchanged but for taking the templated
values (matrix path, methods, preprocessing, prefix and ``args``) as
command-line arguments. ``args`` are upstream's ``ext.args``: ``--<option>
<value>`` pairs overriding the script's defaults (e.g. ``--doHeatmap
false``). As upstream, the matrix is linked into the working directory under
``hto_matrix_name`` (upstream's staged directory name, ``hto`` for the test
data), which the parameters table records.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

BFF_IMAGE = "community.wave.seqera.io/library/bioconductor-cellhashr_r-seurat:25c4bc76749af5ac"
# Upstream: process_medium; the R code is single-threaded.
DEFAULT_RESOURCES = flyte.Resources(cpu=2, memory="8Gi")

BFF_R = r"""
#!/usr/bin/env Rscript

cli_args <- commandArgs(trailingOnly = TRUE)

################################################
################################################
## Fucntions                                  ##
################################################
################################################

# Helper function for NULL condition
null_to_string <- function(x, val = "NULL") if (is.null(x)) val else x

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
    # module input parameters
    hto_matrix                  = cli_args[1],    # Path to the HTO matrix file
    methods                     = cli_args[2],       # Methods to use: COMBINED, RAW, CLUSTER
    preprocessing               = cli_args[3], # Preprocessing method

    # ext.args for preprocessing
    barcodeWhitelist            = NULL,             # A vector of barcode names to retain, used for preprocessing step

    # ext.args for GenerateCellHashingCalls()
    cellbarcodeWhitelist        = NULL,             # A vector of expected cell barcodes. Allows reporting on the total set of expected barcodes, not just those in the filtered count matrix
    methodsForConsensus         = NULL,             # By default, a consensus call will be generated using all methods, NULL, RAW or CLUSTER
    metricsFile                 = NULL,             # Path to metrics file (output)
    doTSNE                      = FALSE,            # If true, tSNE will be run on the resulting hashing calls after each caller
    doHeatmap                   = TRUE,             # If true, Seurat::HTOHeatmap will be run on the results of each caller
    perCellSaturation           = NULL,             # An optional dataframe with the columns cellbarcode and saturation
    majorityConsensusThreshold  = NULL,             # This applies to calculating a consensus call when multiple algorithms are used
    chemistry                   = "10xV3",          # This string is passed to EstimateMultipletRate. Should be either 10xV2 or 10xV3
    callerDisagreementThreshold = NULL,             # If provided, the agreement rate will be calculated between each caller and the simple majority call, ignoring discordant and no-call cells

    prefix                      = cli_args[4] # Prefix name for output files
)
opt_types <- lapply(opt, class)

# Apply parameter overrides
args_string <- cli_args[5]
args_opt <- if (is_valid_string(args_string)) parse_args(args_string) else list()
for ( ao in names(args_opt)){
    if (! ao %in% names(opt)){
        stop(paste("Invalid option:", ao))
    }else{
        if (! is.null(opt[[ao]])){
            # Preserve classes from defaults where possible
            opt[[ao]] <- as(args_opt[[ao]], opt_types[[ao]])
        } else {
            opt[[ao]] <- args_opt[[ao]]
        }
    }
}

# Set individual variables for backward compatibility and cleaner code
hto_matrix                  <- opt$hto_matrix
methods                     <- opt$methods
preprocessing               <- opt$preprocessing
barcodeWhitelist            <- opt$barcodeWhitelist
cellbarcodeWhitelist        <- opt$cellbarcodeWhitelist
methodsForConsensus         <- opt$methodsForConsensus
metricsFile                 <- paste0(opt$prefix, "_metrics_bff.csv")
doTSNE                      <- opt$doTSNE
doHeatmap                   <- opt$doHeatmap
perCellSaturation           <- opt$perCellSaturation
majorityConsensusThreshold  <- opt$majorityConsensusThreshold
chemistry                   <- opt$chemistry
callerDisagreementThreshold <- opt$callerDisagreementThreshold
prefix                      <- opt$prefix

# check if the file exists
if (! file.exists(hto_matrix)){
    stop(paste0(hto_matrix, ' is not a valid file'))
}

################################################
################################################
## Finish loading libraries                   ##
################################################
################################################

library(Seurat) # for Read10X()
library(cellhashR) # for ProcessCountMatrix() and GenerateCellHashingCalls()

################################################
################################################
## Main Process                               ##
################################################
################################################

# perform a basic preprocessing or direcly read the data
if(preprocessing){

  if(is.null(barcodeWhitelist)){
    barcodes_list <- NULL
  }else{
    # separate the barcodes by comma and remove leading/trailing whitespace from each word
    barcodes_list <- trimws(strsplit(barcodeWhitelist, ",")[[1]])
  }

  # perform preprocessing
  counts <- ProcessCountMatrix(rawCountData = hto_matrix, barcodeWhitelist = barcodes_list)

  # Add '-1' suffix back to match original 10X barcodes (is removed by calling Read10X(strip.suffix = TRUE) in ProcessCountMatrix())
  barcodes <- readLines(gzfile(paste0(hto_matrix,"/barcodes.tsv.gz")))
  if(all(grepl("\\-1$", barcodes))) {
    colnames(counts) <- paste0(colnames(counts), "-1")
  }
}else{
    # perform preprocessing
  counts <- Read10X(hto_matrix)
}

# determine the methods that should be used
methods_input <- switch(
  methods,
  "COMBINED" = c("bff_raw", "bff_cluster"),
  "RAW" = c("bff_raw"),
  "CLUSTER" = c("bff_cluster"),
  stop("'methods' must be either 'RAW', 'CLUSTER', or 'COMBINED'")
)

# determine the final consensus call based on RAW, CLUSTER or BOTH (NULL)
methodsForConsensus_input = NULL
if (!is.null(methodsForConsensus)) {
    methodsForConsensus_input <- switch(methodsForConsensus,
        "RAW" = "bff_raw",
        "CLUSTER" = "bff_cluster",
        stop("'methodsForConsensus' must be either 'RAW', 'CLUSTER', or NULL")
    )
}

results <- GenerateCellHashingCalls(barcodeMatrix = counts, methods = methods_input, doTSNE = doTSNE, doHeatmap = doHeatmap, methodsForConsensus = methodsForConsensus_input,cellbarcodeWhitelist = cellbarcodeWhitelist,metricsFile = metricsFile, perCellSaturation = perCellSaturation, majorityConsensusThreshold = majorityConsensusThreshold, chemistry = chemistry, callerDisagreementThreshold = callerDisagreementThreshold)

################################################
################################################
## SAVING RESULTS                             ##
################################################
################################################

#saving parameters in a dataframe
Argument <- c("hto_matrix", "methods", "methodsForConsensus", "cellbarcodeWhitelist", "metricsFile", "perCellSaturation","majorityConsensusThreshold","callerDisagreementThreshold", "doTSNE","doHeatmap","chemistry")
Value <- c(hto_matrix, null_to_string(methods), null_to_string(methodsForConsensus), null_to_string(cellbarcodeWhitelist), metricsFile, null_to_string(perCellSaturation), null_to_string(majorityConsensusThreshold), null_to_string(callerDisagreementThreshold), doTSNE, doHeatmap, chemistry)
params <- data.frame(Argument, Value)
write.csv(params, paste0(prefix ,"_params_bff.csv"))

# save the results of GenerateCellHashingCalls or an empty dataframe if NULL
if(is.null(results)){
    results <- data.frame()
}
write.csv(results, paste0(prefix, "_assignment_bff.csv"), row.names=FALSE)

################################################
################################################
## VERSIONS FILE                              ##
################################################
################################################

r.version <- paste(R.version[['major']],R.version[['minor']], sep = ".")
seurat.version <- as.character(packageVersion('Seurat'))
cellhashR.version <- as.character(packageVersion('cellhashR'))

writeLines(
    c(
        '"BFF":',
        paste('    r-base:', r.version),
        paste('    r-seurat:', seurat.version),
        paste('    cellhashR:', cellhashR.version)
    ),
'versions.yml')

################################################
################################################
################################################
################################################
""".strip()  # noqa: E501

bff_cmd = shell.create(
    name="bff",
    image=BFF_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "hto_matrix": Dir,
        "hto_matrix_name": str,
        "methods": str,
        "preprocessing": str,
        "prefix": str,
        "args": str,
    },
    defaults={"hto_matrix_name": "hto", "args": "null"},
    outputs={"results": Dir},
    script=r"""
        MATRIX={inputs.hto_matrix}
        NAME={inputs.hto_matrix_name}
        P={inputs.prefix}
        OUT={outputs.results}
        W=$(mktemp -d); cd "$W"
        ln -s "$MATRIX" "$NAME"
        cat > bff.R <<'END_OF_SCRIPT'
BFF_R
END_OF_SCRIPT
        Rscript bff.R "$NAME" {inputs.methods} {inputs.preprocessing} "$P" {inputs.args} || exit 1
        mv "${P}_assignment_bff.csv" "${P}_metrics_bff.csv" "${P}_params_bff.csv" versions.yml "$OUT"/
    """.replace("BFF_R", BFF_R),
)

env = flyte.TaskEnvironment.from_task("bff", bff_cmd.as_task())


@dataclass
class BffResult:
    assignment: File  # <prefix>_assignment_bff.csv
    metrics: File  # <prefix>_metrics_bff.csv
    params: File  # <prefix>_params_bff.csv
    versions: File  # versions.yml


async def bff(
    hto_matrix: Dir,
    prefix: str,
    methods: str = "COMBINED",
    preprocessing: bool = True,
    hto_matrix_name: str = "hto",
    args: str = "",
) -> BffResult:
    """BFF hashtag calls for a 10x HTO matrix directory; ``methods`` is ``RAW``, ``CLUSTER`` or ``COMBINED``."""
    results = await bff_cmd(
        hto_matrix=hto_matrix,
        hto_matrix_name=hto_matrix_name,
        methods=methods,
        preprocessing="TRUE" if preprocessing else "FALSE",
        prefix=prefix,
        # As upstream's template renders an unset ext.args.
        args=args or "null",
    )

    async def pick(name: str) -> File:
        found = await results.get_file(name)
        if found is None:
            raise FileNotFoundError(f"bff wrote no {name}")
        return found

    return BffResult(
        assignment=await pick(f"{prefix}_assignment_bff.csv"),
        metrics=await pick(f"{prefix}_metrics_bff.csv"),
        params=await pick(f"{prefix}_params_bff.csv"),
        versions=await pick("versions.yml"),
    )
