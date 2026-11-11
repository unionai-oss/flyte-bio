"""Tests for flyte_bio.modules.melon.

The case and expected md5s are upstream's melon test (modules commit
9a1cc9f, melon 0.2.5): Bacteroides fragilis nanopore reads profiled against
upstream's test database (``delete_me/melon_db.tar.gz``, extracted with
:func:`fixture_dir` as upstream untars it), without a Kraken2 database. As
upstream, the log is checked by name.
"""

from flyte_bio.modules.melon import melon
from tests.framework import assert_md5, env, fixture, fixture_dir


@env.task
async def test_melon_nanopore() -> None:
    # upstream case: melon "bacteroides_fragilis - nanopore fastq"
    reads = await fixture("genomics/prokaryotes/bacteroides_fragilis/nanopore/fastq/test.fastq.gz")
    database = await fixture_dir("delete_me/melon_db.tar.gz", "melon_test_db")
    result = await melon(reads, database, prefix="test")
    await assert_md5(result.tsv, "a040c7a1c87d39cebb17896f6f869328", label="test.tsv")
    await assert_md5(result.json, "5254e0dbe3aa85dbeda8afe934c08d4e", label="test.json")
    assert result.log.name == "test.log", result.log.name


tests = [test_melon_nanopore]
