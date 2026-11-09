"""Tests for flyte_bio.modules.genomescope2.

The case and expected md5s are upstream's genomescope2 test (modules commit
9a1cc9f, GenomeScope 2.1.0) with its test config's ``-k 21 -p 1
--fitted_hist``, after the meryl count and histogram steps upstream sets up.
As upstream, the plots are only checked to exist.
"""

from flyte_bio.modules.genomescope2 import genomescope2
from flyte_bio.modules.meryl import meryl_count, meryl_histogram
from tests.framework import assert_md5, assert_nonempty, env, fixture

FRAGILIS = "genomics/prokaryotes/bacteroides_fragilis/illumina/fastq/test1_1.fastq.gz"


@env.task
async def test_genomescope2_fastq() -> None:
    dbs = await meryl_count([await fixture(FRAGILIS)], kvalue=21, prefix="test")
    hist = await meryl_histogram(dbs[0], kvalue=21, prefix="test")
    result = await genomescope2(hist, prefix="test", args="-k 21 -p 1 --fitted_hist")
    for plot in [
        result.linear_plot_png,
        result.transformed_linear_plot_png,
        result.log_plot_png,
        result.transformed_log_plot_png,
    ]:
        await assert_nonempty(plot, label=plot.name)
    assert result.fitted_histogram_png is not None, "genomescope2 wrote no fitted histogram"
    assert result.lookup_table is not None, "genomescope2 wrote no lookup table"
    await assert_md5(result.model, "b452ed84cc82ec5230d34ebe1c85b703", label="test_model.txt")
    await assert_md5(result.summary, "6d4b102b60643c9fdddb8b6931a5e5eb", label="test_summary.txt")
    await assert_md5(result.lookup_table, "0b2b660e343ba7139dceb9238429c4c9", label="test_lookup_table.txt")


tests = [test_genomescope2_fastq]
