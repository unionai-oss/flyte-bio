"""Tests for flyte_bio.modules.fmhfunprofiler.

The case is upstream's fmhfunprofiler test (modules commit 9a1cc9f,
fmh-funprofiler 1.1.1): single-end reads profiled at k=11, scale 1000 against
a sourmash sketch of the SARS-CoV-2 genome (sourmash sketch's defaults, as
upstream's setup step), with the table checked by md5.
"""

from flyte_bio.modules.fmhfunprofiler import fmhfunprofiler
from flyte_bio.modules.sourmash import sourmash_sketch
from tests.framework import assert_md5, env, fixture


@env.task
async def test_fmhfunprofiler_single_end_fastq() -> None:
    sketch = await sourmash_sketch([await fixture("genomics/sarscov2/genome/genome.fasta")], prefix="test")
    reads = await fixture("genomics/sarscov2/illumina/fastq/test_1.fastq.gz")
    table = await fmhfunprofiler(reads, sketch, 11, 1000, prefix="test")
    assert table.name == "test.csv", table.name
    await assert_md5(table, "16b7196a90c1d2fd3e3109c8d252cd08", label="csv")


tests = [test_fmhfunprofiler_single_end_fastq]
