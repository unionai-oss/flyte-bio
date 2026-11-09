"""Tests for flyte_bio.modules.hasheddrops.

The case and expected values are upstream's hasheddrops test (modules commit
9a1cc9f): the HTO and RNA matrices, without emptyDrops. As upstream, the
matrices are staged as ``hto`` and ``rna``; every output is compared by md5,
the plots and RDS files included.
"""

from flyte_bio.modules.hasheddrops import hasheddrops
from tests.framework import assert_md5, env, fixture_dir

DATA = "genomics/homo_sapiens/10xgenomics/cellranger/hashing_demultiplexing/"


@env.task
async def test_hasheddrops_hto_rna() -> None:
    # upstream case: hasheddrops "hto_matrix - runEmptyDrops:false - rna_matrix"
    hto = await fixture_dir(DATA + "hto.tar.gz", "hto")
    rna = await fixture_dir(DATA + "rna.tar.gz", "rna")
    result = await hasheddrops(hto, prefix="test", rna_matrix=rna, run_empty_drops=False)
    expected = {
        "emptyDrops.csv": (result.empty_drops_csv, "a227d42afbcc590b4e949075cde4a5b6"),
        "emptyDrops.png": (result.empty_drops_plot, "2d7113a76a494b2b2c18b77c4a664063"),
        "emptyDrops.rds": (result.empty_drops_rds, "bda2e88b1e3e8c8ea5a22e19e3ab243d"),
        "id_to_hash.csv": (result.id_to_hash, "a015ca1178177038db21a97fa2070874"),
        "params_hasheddrops.csv": (result.params, "144fd97597ea5df4ad7bdddf2f6fabac"),
        "plot_hasheddrops.png": (result.plot, "2d7113a76a494b2b2c18b77c4a664063"),
        "hasheddrops.rds": (result.rds, "31ba0ecae1d94080bccf81cdbaa5afe8"),
        "results_hasheddrops.csv": (result.results, "d0f2af6d9bdb47a80bfbf08a85493d2e"),
        "versions.yml": (result.versions, "82ab6a6b5ac26c697551e18b1526daa5"),
    }
    for label, (file, md5) in expected.items():
        await assert_md5(file, md5, label=label)


tests = [test_hasheddrops_hto_rna]
