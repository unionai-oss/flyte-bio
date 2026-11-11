"""geoquery — fetch a GEO series' expression data with GEOquery.

Exposes :func:`geoquery_getgeo`, which downloads a GEO series (``GSE...``)
with GEOquery's ``getGEO``, so the task needs network access to NCBI GEO.
Outputs keep upstream's names, after an optional prefix that's prepended
as is (upstream's tests use ``<id>.``):

- ``<prefix>eset.rds``: the series' first ExpressionSet (log2(x + 1)
  transformed if its maximum exceeds 20, as upstream);
- ``<prefix>matrix.tsv``: its expression matrix;
- ``<prefix>annotation.tsv``: its probe annotation (only the comma-separated
  columns ``--metacols`` lists in ``args``, if given);
- ``versions.yml``.

The R code is upstream's template, unchanged but for taking the templated
values (series, ``args``, prefix) as command-line arguments.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

GEOQUERY_IMAGE = "community.wave.seqera.io/library/bioconductor-geoquery_r-base:8547a00389cd568a"
# Upstream: process_single.
DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory="6Gi")

GETGEO_R = r"""
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

    df[,columns] <- format(
        data.frame(df[, columns], check.names = FALSE),
        nsmall = digits
    )

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

opt <- list(
    querygse = cli_args[1],
    metacols = NULL
)
args_opt <- parse_args(cli_args[2])
for ( ao in names(args_opt)){
    if (! ao %in% names(opt)){
        stop(paste("Invalid option:", ao))
    }else{
        opt[[ao]] <- args_opt[[ao]]
    }
}

################################################
################################################
## Finish loading libraries                   ##
################################################
################################################

library(GEOquery)

################################################
################################################
## Do the GEO query retrieval                 ##
################################################
################################################

# Fetch data for GSE number

eset <- getGEO(
    GEO = opt$querygse,
    destdir = getwd()
)[[1]]

# Write probeset annotation. If supplied, Parse metadata columns from nextflow
# parameters to subset on the feature metadata file

probeset_annotation = fData(eset)
if (! is.null(opt$metacols)){
    feature_cols = strsplit(opt$metacols,',')[[1]]
    probeset_annotation <- probeset_annotation[,feature_cols]
}

################################################
################################################
## Generate outputs                           ##
################################################
################################################

output_prefix <- ifelse(cli_args[3] == 'null', '', cli_args[3])

write.table(
    probeset_annotation,
    paste0(output_prefix,'annotation.tsv'),
    col.names=TRUE,
    row.names=FALSE,
    sep="\t",
    quote=FALSE
)

# If data is not log scale, transform it as needed for limma downstream

if(max(exprs(eset),na.rm=T) > 20) { # a bit dirty, needs proper solution later...
    exprs(eset) <- log2(exprs(eset) + 1)
}

saveRDS(eset, file = paste0(output_prefix, 'eset.rds'))

# Write intensity matrix (normalised)

write.table(
    data.frame(
        probe_id = rownames(eset),
        round_dataframe_columns(as.data.frame(exprs(eset))),
        check.names = FALSE
    ),
    file = paste0(output_prefix, 'matrix.tsv'),
    col.names = TRUE, row.names = FALSE,
    sep = '\t', quote = FALSE
)

################################################
################################################
## R SESSION INFO                             ##
################################################
################################################

sink(paste(output_prefix, "R_sessionInfo.log", sep = '.'))
print(sessionInfo())
sink()

################################################
################################################
## VERSIONS FILE                              ##
################################################
################################################

r.version <- strsplit(version[['version.string']], ' ')[[1]][3]
geoquery.version <- as.character(packageVersion("GEOquery"))

writeLines(
    c(
        '"GEOQUERY_GETGEO":',
        paste('    r-base:', r.version),
        paste('    bioconductor-geoquery:', geoquery.version)
    ),
    'versions.yml')

################################################
################################################
################################################
################################################
""".strip()

geoquery_getgeo_cmd = shell.create(
    name="geoquery_getgeo",
    image=GEOQUERY_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"querygse": str, "prefix": str, "args": str},
    outputs={"results": Dir},
    script=r"""
        P={inputs.prefix}
        OUT={outputs.results}
        W=$(mktemp -d); cd "$W"
        cat > getgeo.R <<'END_OF_SCRIPT'
GETGEO_R
END_OF_SCRIPT
        Rscript getgeo.R {inputs.querygse} {inputs.args} "$P" || exit 1
        if [ "$P" = null ]; then P=""; fi
        mv "${P}eset.rds" "${P}matrix.tsv" "${P}annotation.tsv" versions.yml "$OUT"/
    """.replace("GETGEO_R", GETGEO_R),
)

env = flyte.TaskEnvironment.from_task("geoquery", geoquery_getgeo_cmd.as_task())


@dataclass
class GetGeoResult:
    rds: File  # <prefix>eset.rds
    expression: File  # <prefix>matrix.tsv
    annotation: File  # <prefix>annotation.tsv
    versions: File  # versions.yml


async def geoquery_getgeo(querygse: str, prefix: str = "", args: str = "") -> GetGeoResult:
    """GEO series ``querygse``'s expression set, matrix and probe annotation."""
    # As upstream's template renders an unset ext.args and ext.prefix.
    results = await geoquery_getgeo_cmd(querygse=querygse, prefix=prefix or "null", args=args or "null")

    async def pick(name: str) -> File:
        found = await results.get_file(name)
        if found is None:
            raise FileNotFoundError(f"getGEO wrote no {name}")
        return found

    return GetGeoResult(
        rds=await pick(f"{prefix}eset.rds"),
        expression=await pick(f"{prefix}matrix.tsv"),
        annotation=await pick(f"{prefix}annotation.tsv"),
        versions=await pick("versions.yml"),
    )
