"""Tests for flyte_bio.modules.cnaqc.

The cases are upstream's cnaqc tests (modules commit 9a1cc9f, CNAqc 1.1.4),
with their test config's ``ext.args`` plus ``--qc_chr false`` / ``true``. As
upstream, the RDS and PDF outputs are checked by name and presence only (they
aren't byte-stable), and ``versions.yml`` by md5.
"""

from flyte_bio.modules.cnaqc import cnaqc
from tests.framework import assert_md5, assert_nonempty, env, fixture

DATA = "delete_me/cnaqc/HCC1395_WXS_HCC1395_HCC1395_HCC1395T_"
# Upstream's test config.
ARGS = (
    "--genome GRCh38 --matching_strategy rightmost --purity_error 0.05 --muts_per_karyotype 25 "
    "--min_absolute_karyotype_mutations 100 --min_VAF 0 --blacklist_indels TRUE"
)
VERSIONS_MD5 = "a847981274876c417351abe4568d0447"


async def run_cnaqc(qc_chr: str):
    snv = await fixture(DATA + "snv.rds")
    cna = await fixture(DATA + "cna.rds")
    result = await cnaqc(snv, cna, "HCC1395_HCC1395T", prefix="test", args=f"{ARGS} --qc_chr {qc_chr}")
    for file in (result.qc_rds, result.data_plot_rds, result.qc_plot_rds, result.data_pdf, result.qc_pdf):
        await assert_nonempty(file, label=file.name)
    await assert_md5(result.versions, VERSIONS_MD5, label="versions.yml")
    return result


@env.task
async def test_cnaqc_qc_chr_false() -> None:
    # upstream case: "test cnaqc qc_chr false"
    result = await run_cnaqc("false")
    assert result.qc_by_chr_rds is None and result.qc_by_chr_plot_rds is None and result.qc_by_chr_pdf is None


@env.task
async def test_cnaqc_qc_chr_true() -> None:
    # upstream case: "test cnaqc qc_chr true"
    result = await run_cnaqc("true")
    for file in (result.qc_by_chr_rds, result.qc_by_chr_plot_rds, result.qc_by_chr_pdf):
        assert file is not None, "cnaqc wrote no per-chromosome output"
        await assert_nonempty(file, label=file.name)


tests = [test_cnaqc_qc_chr_false, test_cnaqc_qc_chr_true]
