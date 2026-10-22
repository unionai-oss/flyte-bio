#!/usr/bin/env python
"""Convert a simpleaf quantification to an AnnData h5ad (vendored from upstream mtx_to_h5ad, simpleaf).

Core logic is upstream's, unchanged; the templated values are replaced with
an argparse CLI (--input: the af_quant directory, --sample, --output);
the versions.yml output and the empty per-sample directory are dropped.
"""

# Set numba chache dir to current working directory (which is a writable mount also in containers)
import os

os.environ["NUMBA_CACHE_DIR"] = "."

import argparse
import json

import scanpy as sc


def input_to_adata(
    input_data: str,
    output: str,
    sample: str,
):
    print(f"Reading in {input_data}")

    # open main data
    simpleaf_h5ad_path = f"{input_data}/alevin/quants.h5ad"

    # the simpleaf quant module exports an h5ad file.
    adata = sc.read_h5ad(simpleaf_h5ad_path)
    adata.obs_names = adata.obs["barcodes"].values
    adata.var_names = adata.var["gene_id"].values
    adata.obs["sample"] = sample

    # sort adata column- and row- wise to avoid positional differences
    adata = adata[adata.obs_names.sort_values(), adata.var_names.sort_values()].copy()

    # standard format
    # index are gene IDs and symbols are a column
    adata.var["gene_versions"] = adata.var["gene_id"]
    adata.var.index = adata.var["gene_versions"].str.split(".").str[0].values
    adata.var_names_make_unique()

    # sort adata column- and row- wise to avoid positional differences
    adata = adata[adata.obs_names.sort_values(), adata.var_names.sort_values()].copy()

    # Remove runtime to prevent hash changes
    simpleaf_map_info = json.loads(adata.uns["simpleaf_map_info"])
    simpleaf_map_info.pop("runtime_seconds")
    adata.uns["simpleaf_map_info"] = json.dumps(simpleaf_map_info, sort_keys=True)

    # write results
    adata.write_h5ad(f"{output}")
    print(f"Wrote h5ad file to {output}")


#
# Run main script
#

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="simpleaf af_quant directory")
    parser.add_argument("--sample", required=True, help="sample ID, recorded in obs['sample']")
    parser.add_argument("--output", required=True, help="output h5ad path")
    args = parser.parse_args()

    input_to_adata(input_data=args.input, output=args.output, sample=args.sample)
