#!/usr/bin/env Rscript

# Vendored from upstream anndatar_convert. Core logic is upstream's; only the
# templated values are replaced with positional arguments: the input h5ad and
# the output prefix (written as <prefix>.seurat.rds and <prefix>.sce.rds,
# upstream's <id>_<input_type>_matrix). The versions.yml output and the empty
# per-sample directory are dropped.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
    stop("usage: anndatar_convert.r <input.h5ad> <output prefix>")
}
h5ad <- args[1]
prefix <- args[2]

# load libraries
library(anndataR)
library(SeuratObject)
library(SingleCellExperiment)

# read input
adata <- read_h5ad(h5ad)

# convert to Seurat
obj <- adata$as_Seurat()

# save files
saveRDS(obj, file = paste0(prefix, ".seurat.rds"))

# convert to SingleCellExperiment
obj <- adata$as_SingleCellExperiment()

# save files
saveRDS(obj, file = paste0(prefix, ".sce.rds"))
