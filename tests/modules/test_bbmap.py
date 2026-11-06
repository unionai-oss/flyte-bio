"""Tests for flyte_bio.modules.bbmap.

Mirrors the upstream bbmap/bbsplit "split with prebuilt index" case: build
an index, then split reads with it. Expected md5s are the upstream snapshot
values (nf-test hashes the decompressed reads). The log check proves the
staged index was reused rather than rebuilt — the timestamp/path fix-up
upstream applies to a staged index is what makes that work.

The bbmap/bbduk cases are upstream's (modules commit 9a1cc9f, BBMap 39.18),
with their test config's ``trimq=10 qtrim=r`` and ``<id>.trim`` prefix. The
contaminant cases' snapshots record only the version, so they check the
output names and the log, as upstream asserts.
"""

from flyte.io import File

from flyte_bio.modules.bbmap import BBDukResult, bbduk, bbsplit, bbsplit_index
from tests.framework import assert_gunzipped_md5, assert_md5, assert_nonempty, env, fixture

FASTQ = "genomics/sarscov2/illumina/fastq/"
BBDUK_ARGS = "trimq=10 qtrim=r"


@env.task
async def test_bbsplit_prebuilt_index() -> None:
    # upstream case: bbmap/bbsplit "sarscov2_se_fastq_fasta_chr22_fasta - split with prebuilt index"
    primary = await fixture("genomics/sarscov2/genome/genome.fasta")
    human = await fixture("genomics/homo_sapiens/genome/chr22/sequence/chr22_23800000-23980000.fa")
    index = await bbsplit_index(primary, {"human": human})

    reads = await fixture("genomics/sarscov2/illumina/fastq/test_1.fastq.gz")
    result = await bbsplit(reads, None, index, prefix="test")
    await assert_gunzipped_md5(result.reads_1, "4161df271f9bfcd25d5845a1e220dbec", label="primary reads")
    await assert_md5(result.stats, "2cbf69b72e5f4f8508306b54e8fe2861", label="stats")
    async with result.log.open("rb") as fh:
        log = bytes(await fh.read()).decode()
    assert "If you wish to regenerate the index" in log, "bbsplit rebuilt the index instead of reusing it"


async def check_bbduk(result: BBDukResult, names: list[str], md5s: list[str] | None, mode: str) -> None:
    assert [f.name for f in result.reads] == names, [f.name for f in result.reads]
    for file, md5 in zip(result.reads, md5s or [None] * len(names), strict=True):
        if md5 is None:
            await assert_nonempty(file, label=file.name)
        else:
            await assert_gunzipped_md5(file, md5, label=file.name)
    async with result.log.open("rb") as fh:
        log = bytes(await fh.read()).decode()
    assert f"Input is being processed as {mode}" in log, f"bbduk log lacks 'processed as {mode}'"


async def bbduk_reads() -> tuple[File, File]:
    return await fixture(FASTQ + "test_1.fastq.gz"), await fixture(FASTQ + "test_2.fastq.gz")


@env.task
async def test_bbduk_single_end() -> None:
    # upstream case: bbmap/bbduk "sarscov2 - single end fastq - fastq"
    reads_1, _ = await bbduk_reads()
    result = await bbduk(reads_1, prefix="test.trim", args=BBDUK_ARGS)
    await check_bbduk(result, ["test.trim.fastq.gz"], ["4161df271f9bfcd25d5845a1e220dbec"], "unpaired")


@env.task
async def test_bbduk_paired_end() -> None:
    # upstream case: bbmap/bbduk "sarscov2 - paired end fastq - fastq"
    result = await bbduk(*await bbduk_reads(), prefix="test.trim", args=BBDUK_ARGS)
    await check_bbduk(
        result,
        ["test.trim_1.fastq.gz", "test.trim_2.fastq.gz"],
        ["4161df271f9bfcd25d5845a1e220dbec", "2ebae722295ea66d84075a3b042e2b42"],
        "paired",
    )


@env.task
async def test_bbduk_single_end_contaminants() -> None:
    # upstream case: bbmap/bbduk "sarscov2 - single end w/ contams [fastq,fasta] - fastq"
    reads_1, _ = await bbduk_reads()
    contaminants = await fixture("genomics/sarscov2/genome/transcriptome.fasta")
    result = await bbduk(reads_1, contaminants=contaminants, prefix="test.trim", args=BBDUK_ARGS)
    await check_bbduk(result, ["test.trim.fastq.gz"], None, "unpaired")


@env.task
async def test_bbduk_paired_end_contaminants() -> None:
    # upstream case: bbmap/bbduk "sarscov2 - paired end w/ contams [fastq,fasta] - fastq"
    contaminants = await fixture("genomics/sarscov2/genome/transcriptome.fasta")
    result = await bbduk(*await bbduk_reads(), contaminants=contaminants, prefix="test.trim", args=BBDUK_ARGS)
    await check_bbduk(result, ["test.trim_1.fastq.gz", "test.trim_2.fastq.gz"], None, "paired")


tests = [
    test_bbsplit_prebuilt_index,
    test_bbduk_single_end,
    test_bbduk_paired_end,
    test_bbduk_single_end_contaminants,
    test_bbduk_paired_end_contaminants,
]
