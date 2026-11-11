"""Tests for flyte_bio.modules.krakentools.

The case and expected md5 are upstream's krakentools/kreport2krona test
(modules commit 9a1cc9f, KrakenTools 1.2.1): the Kraken2 report of
test_1.fastq.gz against upstream's test database, which the existing kraken2
module makes here as in upstream's setup.

The combine case is upstream's krakentools/combinekreports test (same
commit): the test data's Kraken2 report of test_1.fastq.gz, combined alone.
"""

from flyte_bio.modules.kraken2 import kraken2
from flyte_bio.modules.krakentools import krakentools_combinekreports, krakentools_kreport2krona
from tests.framework import assert_md5, env, fixture, fixture_dir


@env.task
async def test_krakentools_kreport2krona_single_end() -> None:
    # upstream case: krakentools/kreport2krona "sarscov2 illumina single end [fastq]"
    reads = await fixture("genomics/sarscov2/illumina/fastq/test_1.fastq.gz")
    db = await fixture_dir("genomics/sarscov2/genome/db/kraken2.tar.gz", "kraken2")
    report = (await kraken2(reads, db, prefix="test")).report
    txt = await krakentools_kreport2krona(report, prefix="test")
    assert txt.name == "test.txt", txt.name
    await assert_md5(txt, "c89a9db7acbdba9dea0fe246bcaa85c1", label="test.txt")


@env.task
async def test_krakentools_combinekreports() -> None:
    # upstream case: krakentools/combinekreports "sarscov2 - metagenome - kraken report"
    report = await fixture("genomics/sarscov2/metagenome/test_1.kraken2.report.txt")
    combined = await krakentools_combinekreports([report], prefix="test")
    assert combined.name == "test.txt", combined.name
    await assert_md5(combined, "c231bfcdde47a43336a418b21bb2f98d", label="test.txt")


tests = [test_krakentools_kreport2krona_single_end, test_krakentools_combinekreports]
