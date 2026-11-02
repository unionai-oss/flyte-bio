"""Tests for flyte_bio.modules.anndata.

Upstream's test: three barcodes picked out of a filtered PBMC matrix, with
the output h5ad's md5 snapshotted. The h5ad's bytes don't match that md5
here (HDF5), so the test checks its content against the input instead.
"""

import tempfile
from pathlib import Path

import anndata
import numpy as np
from flyte.io import File

from flyte_bio.modules.anndata import anndata_barcodes
from tests.framework import env, fixture

BARCODES = ["AAACCTGAGATATGGT", "AAACCTGAGCAATATG", "AAACCTGAGCACCGCT"]


@env.task
async def test_anndata_barcodes() -> None:
    h5ad = await fixture("genomics/homo_sapiens/scrnaseq/h5ad/SRR28679759_filtered_matrix.h5ad")
    csv = Path(tempfile.mkdtemp()) / "barcodes.csv"
    csv.write_text("".join(f"{b}\n" for b in BARCODES))
    result = await anndata_barcodes(h5ad, await File.from_local(str(csv)), prefix="test")
    assert result.path.endswith("/test.h5ad"), result.path
    subset = anndata.read_h5ad(await result.download())
    full = anndata.read_h5ad(await h5ad.download())
    assert list(subset.obs_names) == BARCODES, list(subset.obs_names)
    assert list(subset.var_names) == list(full.var_names)
    assert np.array_equal(subset.X.toarray(), full[BARCODES].X.toarray())


tests = [test_anndata_barcodes]
