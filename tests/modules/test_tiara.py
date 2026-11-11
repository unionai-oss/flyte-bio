"""Tests for flyte_bio.modules.tiara.

The case and md5s are upstream's tiara/tiara test (modules commit 9a1cc9f,
tiara 1.0.3), with its test config's ``--tf bac``. As upstream, the log's name
is checked but not its content.
"""

from flyte_bio.modules.tiara import file_name, tiara
from tests.framework import assert_md5, env, fixture


@env.task
async def test_tiara_contigs() -> None:
    # upstream case: tiara/tiara "illumina-contigs-fa"
    contigs = await fixture("genomics/prokaryotes/bacteroides_fragilis/illumina/fasta/test1.contigs.fa.gz")
    result = await tiara(contigs, prefix="test", args="--tf bac")
    await assert_md5(result.classifications, "1e4a4d8081c20b2540316fb9fafb5ad3", label="test.txt")
    assert file_name(result.log) == "log_test.txt", file_name(result.log)
    names = [file_name(f) for f in result.fasta]
    assert names == ["bacteria_test.fasta"], names
    await assert_md5(result.fasta[0], "d46834f0b781bc0277269f79226f5001", label="bacteria_test.fasta")


tests = [test_tiara_contigs]
