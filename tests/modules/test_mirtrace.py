"""Tests for flyte_bio.modules.mirtrace.

The cases and expected md5s are upstream's mirtrace/qc tests (modules commit
9a1cc9f, miRTrace 1.0.1). As upstream, the HTML report and the JSON results
are checked by content, not md5.
"""

import json
import tempfile
from pathlib import Path

from flyte.io import File

from flyte_bio.modules.mirtrace import MirtraceQcResult, mirtrace_qc
from tests.framework import assert_md5, env, fixture

READS = "genomics/homo_sapiens/illumina/fastq/test_rnaseq_1.fastq.gz"
STATS_MD5 = {
    "mirtrace-stats-contamination_basic.tsv": "ac69ca6d2a709854f1048b635d06e927",
    "mirtrace-stats-contamination_detailed.tsv": "ef80997ac12662c64cbcf5fe9851e786",
    "mirtrace-stats-length.tsv": "54ecc0698e8e83d5b1c979b2ee3b1512",
    "mirtrace-stats-mirna-complexity.tsv": "ab2a7600a2daa5c1797eea13d0abc2f0",
    "mirtrace-stats-phred.tsv": "44eaeae26ec629e71fb31e56bfb5a548",
    "mirtrace-stats-qcstatus.tsv": "623cf3a0c5e363488966844feb0dd978",
    "mirtrace-stats-rnatype.tsv": "ed8d3a76247a1432a365def87c3f4c67",
}


def file_name(file: File) -> str:
    return file.path.rstrip("/").rsplit("/", 1)[-1]


async def text(file: File) -> str:
    async with file.open("rb") as fh:
        return bytes(await fh.read()).decode()


async def check(result: MirtraceQcResult, fasta_md5: str, stats_md5: dict[str, str]) -> None:
    assert "This file is part of miRTrace." in await text(result.html), "html report lacks miRTrace's notice"
    stats = json.loads(await text(result.json))["results"][0]["stats"]
    assert stats["uniqueQCPassedSeqsCount"] == 912, f"uniqueQCPassedSeqsCount {stats['uniqueQCPassedSeqsCount']}"
    names = [file_name(f) for f in result.tsv]
    assert names == sorted(stats_md5), names
    for tsv in result.tsv:
        await assert_md5(tsv, stats_md5[file_name(tsv)], label=file_name(tsv))
    for label, fastas in (("all", result.all_fasta), ("rnatype_unknown", result.rnatype_unknown_fasta)):
        assert [file_name(f) for f in fastas] == ["test_rnaseq_1.fasta"], f"{label}: {[f.path for f in fastas]}"
        await assert_md5(fastas[0], fasta_md5, label=f"{label} test_rnaseq_1.fasta")


@env.task
async def test_mirtrace_qc_human_fastq() -> None:
    result = await mirtrace_qc([await fixture(READS)], "hsa")
    await check(result, "0181296141177654088bbc2a96b29560", STATS_MD5)


@env.task
async def test_mirtrace_qc_human_fastq_optional_config() -> None:
    path = Path(tempfile.mkdtemp()) / "mirtrace_config.txt"
    path.write_text("./test_rnaseq_1.fastq.gz,test_rnaseq_1,TGGAATTCTCGGGTGCCAAGG,33\n")
    config = await File.from_local(str(path))
    result = await mirtrace_qc([await fixture(READS)], "hsa", config=config)
    stats_md5 = STATS_MD5 | {"mirtrace-stats-length.tsv": "630b3a845953321ae4cf2fc4a4943ab5"}
    await check(result, "125963c9ee39a49b3a680903e28e7c9d", stats_md5)


tests = [test_mirtrace_qc_human_fastq, test_mirtrace_qc_human_fastq_optional_config]
