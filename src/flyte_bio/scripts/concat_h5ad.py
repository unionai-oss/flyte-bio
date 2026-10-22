#!/usr/bin/env python
"""Concatenate per-sample h5ad files into one, annotated from the samplesheet (vendored from upstream concat_h5ad).

Core logic is upstream's, unchanged: it reads every ``*.h5ad`` under the
working directory (run it from a directory holding the per-sample files,
as upstream stages them) and keys each by its path minus
``_matrix.h5ad``. The templated values are replaced with an argparse CLI
(--samplesheet, --output) and the versions.yml output is dropped.
"""

# Set numba chache dir to current working directory (which is a writable mount also in containers)
import os

os.environ["NUMBA_CACHE_DIR"] = "."

import argparse
from pathlib import Path

import anndata as ad
import pandas as pd
import scanpy as sc


def read_samplesheet(samplesheet):
    df = pd.read_csv(samplesheet)
    df.set_index("sample")

    # samplesheet may contain replicates, when it has,
    # group information from replicates and collapse with commas
    # only keep unique values using set()
    df = df.groupby(["sample"]).agg(lambda column: ",".join(set(column.astype(str))))

    return df


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samplesheet", required=True, help="the pipeline's samplesheet CSV")
    parser.add_argument("--output", required=True, help="output h5ad path")
    args = parser.parse_args()

    # Open samplesheet as dataframe
    df_samplesheet = read_samplesheet(args.samplesheet)

    # find all h5ad and append to dict
    dict_of_h5ad = {str(path).replace("_matrix.h5ad", ""): sc.read_h5ad(path) for path in Path(".").rglob("*.h5ad")}

    # concat h5ad files
    adata = ad.concat(dict_of_h5ad, label="sample", merge="unique", index_unique="_")

    # merge with data.frame, on sample information
    adata.obs = adata.obs.join(df_samplesheet, on="sample", how="left").astype(str)
    adata.write_h5ad(args.output)

    print("Wrote h5ad file to {}".format(args.output))
