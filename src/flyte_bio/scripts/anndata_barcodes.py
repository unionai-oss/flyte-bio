#!/usr/bin/env python3
"""Keep an h5ad's cells whose barcodes are listed in a CSV (vendored from upstream anndata/barcodes).

Core logic is upstream's, unchanged; the templated values are replaced with
an argparse CLI (--h5ad, --barcodes, --output) and the versions.yml output
is dropped.
"""

import argparse

import anndata as ad
import pandas as pd

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--h5ad", required=True, help="input h5ad")
    parser.add_argument("--barcodes", required=True, help="headerless CSV, one barcode per line")
    parser.add_argument("--output", required=True, help="output h5ad path")
    args = parser.parse_args()

    df = pd.read_csv(args.barcodes, header=None)
    adata = ad.read_h5ad(args.h5ad)

    adata = adata[df[0].values]

    adata.write_h5ad(args.output)
