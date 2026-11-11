"""affy — RMA-normalised expression from Affymetrix CEL files.

Exposes :func:`affy_justrma`, which runs affy's ``justRMA`` over the CEL
files a sample sheet lists (``.csv``, ``.tsv`` or ``.txt``; its ``file``
column names them, by default) in ``celfiles_dir``. Outputs keep upstream's
names:

- ``<prefix>_eset.rds``: the ExpressionSet;
- ``<prefix>_matrix.tsv``: the expression matrix (probe IDs, then one column
  per sample, rounded to 8 decimal places);
- ``<cdf>.annotation.tsv``: probe annotation, with
  ``--build_annotation TRUE`` in ``args``;
- ``R_sessionInfo.log`` and ``versions.yml``.

As upstream, the script installs the arrays' CDF package (and, for
annotation, its ``.db`` package) from Bioconductor when it runs, so the task
needs network access. ``args`` are upstream's ``ext.args`` (``--<option>
<value>`` pairs, e.g. ``--sample_name_col name``). The R code is upstream's
template, unchanged but for taking the templated values (sample sheet, CEL
directory, description file, ``args`` and prefix) as command-line arguments.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

AFFY_IMAGE = "community.wave.seqera.io/library/bioconductor-affy_r-base:dd8a5ecd6fc301b3"
# Upstream: process_single.
DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory="8Gi")

AFFY_JUSTRMA_R = r"""
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
    parsed_args[! is.na(parsed_args)]
}

#' Flexibly read CSV or TSV files
#'
#' @param file Input file
#' @param header Passed to read.delim()
#' @param row.names Passed to read.delim()
#'
#' @return output Data frame

read_delim_flexible <- function(file, header = TRUE, row.names = NULL){

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
        row.names = row.names
    )
}

#' Install the right CDF for a given cel file
#'
#' @param celfile A valid path to a CEL file
#' @param annotation Boolean indication whether to install the annotation
#'        package
#'
#' @return output The CDF environment or a list detailing the failed locations.

install_cdf_db <- function(celfile, annotation = FALSE){
    library(affyio)
    headdetails <- read.celfile.header(celfile)
    ref.cdfName <- headdetails[[1]]
    cleaned.cdfName <- cleancdfname(ref.cdfName, addcdf = FALSE)

    exts = 'cdf'
    if (annotation){
        exts <- c(exts, '.db')
    }
    options(timeout=600)
    for (package in paste0(cleaned.cdfName, exts)){
        install.packages(
            package,
            lib = 'libs',
            repos = BiocManager::repositories(),
            dependencies = c("Depends", "Imports")
        )
    }
    cleaned.cdfName
}

#' Round numeric dataframe columns to fixed decimal places by applying
#' formatting and converting back to numerics
#'
#' @param dataframe A data frame
#' @param columns Which columns to round (assumes all of them by default)
#' @param digits How many decimal places to round to?
#'
#' @return output Data frame

round_dataframe_columns <- function(df, columns = NULL, digits = 8){
    if (is.null(columns)){
        columns <- colnames(df)
    }

    df[,columns] <- format(data.frame(df[, columns]), nsmall = digits)

    # Convert columns back to numeric

    for (c in columns) {
        df[[c]][grep("^ *NA$", df[[c]])] <- NA
        df[[c]] <- as.numeric(df[[c]])
    }
    df
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
    sample_file = cli_args[1],
    celfiles_dir = cli_args[2],
    description = cli_args[3],
    file_name_col = 'file',
    sample_name_col = 'file',
    background = TRUE,
    normalize = TRUE,
    bgversion = 2,
    destructive=FALSE,
    cdfname = NULL,
    rm.mask = FALSE,
    rm.outliers = FALSE,
    rm.extra = FALSE,
    build_annotation = FALSE,
    keep.log2 = TRUE
)
if (opt$description == ''){
    opt$description = NULL
}
opt_types <- lapply(opt, class)

# Apply parameter overrides

args_opt <- parse_args(cli_args[4])
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

# Check file inputs are valid

for (file_input in c('sample_file')){
    if (is.null(opt[[file_input]])) {
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

library(affy)

################################################
################################################
## READ IN SAMPLE SHEET WITH CELL FILE NAMES  ##
################################################
################################################

sample.sheet <- read_delim_flexible(file = opt$sample_file)

if (! opt$file_name_col %in% colnames(sample.sheet)){
    stop(paste0("Specified file name column '", opt$file_name_col, "' is not in the sample sheet"))
}

################################################
################################################
## Run justRMA process                        ##
################################################
################################################

# Install the CDF in a path we can definitely write to (for some reason affy
# ignores .libPaths() left to its own devices)

dir.create('libs')
.libPaths('libs')
first_cel <- file.path(opt$celfiles_dir, sample.sheet[[opt$file_name_col]][1])
cdf_name <- install_cdf_db(first_cel, annotation = opt$build_annotation)

# Run the main function

rownames(sample.sheet) <- sample.sheet[[opt$file_name_col]]
eset <- justRMA(
    filenames = sample.sheet[[opt$file_name_col]],
    celfile.path = opt$celfiles_dir,
    phenoData = sample.sheet,
    description = opt$description,
    rm.mask = opt$rm.mask,
    rm.outliers = opt$rm.outliers,
    rm.extra = opt$rm.extra,
    sampleNames = sample.sheet[[opt$sample_name_col]],
    normalize = opt$normalize,
    background = opt$background,
    bgversion = opt$bgversion,
    destructive = opt$destructive,
    cdfname = opt$cdfname
)

# This should happen as part of the above, not sure why it doesn't
sampleNames(eset) <- sample.sheet[[opt$sample_name_col]]

# unlog2 the data, if keep.log2 is set to FALSE
if (!opt$keep.log2) exprs(eset) <- 2**exprs(eset)

################################################
################################################
## Generate outputs                           ##
################################################
################################################

if (opt$build_annotation){

    # Make some annotation

    dbname <- paste0(cdf_name, '.db')
    library(dbname, character.only = TRUE)
    anno <- select(
        get(dbname),
        keys=keys(get(dbname), keytype="PROBEID"),
        columns=c('ENSEMBL', 'ENTREZID', 'SYMBOL', 'GENENAME', 'GENETYPE'),
        keytype="PROBEID"
    )

    # Remove duplicates by probe
    anno <- do.call(
        rbind,
        lapply(
            split(
                anno,
                anno$PROBEID
            ),
            function(x) apply(x, 2, function(y) paste(unique(y), collapse=', '))
        )
    )

    write.table(
        anno,
        file = paste0(cdf_name, '.annotation.tsv'),
        col.names = TRUE,
        row.names = FALSE,
        sep = '\t',
        quote = FALSE
    )
}

# R object for other processes to use

output_prefix <- cli_args[5]
saveRDS(eset, file = paste0(output_prefix, '_eset.rds'))

# Write matrix

write.table(
    data.frame(
        probe_id = rownames(eset),
        round_dataframe_columns(as.data.frame(exprs(eset))),
        check.names = FALSE
    ),
    file = paste0(output_prefix, '_matrix.tsv'),
    col.names = TRUE,
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
affy.version <- as.character(packageVersion('affy'))

writeLines(
    c(
        '"AFFY_JUSTRMA":',
        paste('    r-base:', r.version),
        paste('    bioconductor-affy:', affy.version)
    ),
'versions.yml')

################################################
################################################
################################################
################################################
""".strip()

affy_justrma_cmd = shell.create(
    name="affy_justrma",
    image=AFFY_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"samplesheet": File, "celfiles_dir": Dir, "description": File | None, "prefix": str, "args": str},
    defaults={"args": "null"},
    outputs={"results": Dir},
    script=r"""
        SHEET=({inputs.samplesheet})
        shopt -s nullglob; DESCRIPTION=({inputs.description}); shopt -u nullglob
        P={inputs.prefix}
        OUT={outputs.results}
        W=$(mktemp -d); cd "$W"
        SHEET_NAME=$(basename "${SHEET[0]}"); ln -s "${SHEET[0]}" "$SHEET_NAME"
        DESCRIPTION_NAME=""
        if [ ${#DESCRIPTION[@]} -gt 0 ]; then
            DESCRIPTION_NAME=$(basename "${DESCRIPTION[0]}"); ln -s "${DESCRIPTION[0]}" "$DESCRIPTION_NAME"
        fi
        ln -s {inputs.celfiles_dir} celfiles
        cat > affy_justrma.R <<'END_OF_SCRIPT'
AFFY_JUSTRMA_R
END_OF_SCRIPT
        Rscript affy_justrma.R "$SHEET_NAME" celfiles "$DESCRIPTION_NAME" {inputs.args} "$P" || exit 1
        shopt -s nullglob
        mv "${P}_eset.rds" "${P}_matrix.tsv" ./*.annotation.tsv R_sessionInfo.log versions.yml "$OUT"/
    """.replace("AFFY_JUSTRMA_R", AFFY_JUSTRMA_R),
)

env = flyte.TaskEnvironment.from_task("affy", affy_justrma_cmd.as_task())


@dataclass
class AffyJustRmaResult:
    rds: File  # <prefix>_eset.rds
    expression: File  # <prefix>_matrix.tsv
    annotation: File | None  # <cdf>.annotation.tsv, with --build_annotation TRUE
    session_info: File  # R_sessionInfo.log
    versions: File  # versions.yml


async def affy_justrma(
    samplesheet: File, celfiles_dir: Dir, prefix: str, description: File | None = None, args: str = ""
) -> AffyJustRmaResult:
    """RMA expression of the CEL files ``samplesheet`` lists in ``celfiles_dir``."""
    # As upstream's template renders an unset ext.args.
    results = await affy_justrma_cmd(
        samplesheet=samplesheet, celfiles_dir=celfiles_dir, description=description, prefix=prefix, args=args or "null"
    )

    async def pick(name: str) -> File:
        found = await results.get_file(name)
        if found is None:
            raise FileNotFoundError(f"affy justRMA wrote no {name}")
        return found

    annotations = [f for f in await results.list_files() if f.name.endswith(".annotation.tsv")]
    return AffyJustRmaResult(
        rds=await pick(f"{prefix}_eset.rds"),
        expression=await pick(f"{prefix}_matrix.tsv"),
        annotation=annotations[0] if annotations else None,
        session_info=await pick("R_sessionInfo.log"),
        versions=await pick("versions.yml"),
    )
