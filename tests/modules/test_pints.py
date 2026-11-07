"""Tests for flyte_bio.modules.pints.

The cases and expected md5s are upstream's pints/caller tests (modules commit
9a1cc9f, PINTS 1.1.13), with the ``--disable-small`` their config sets.
Upstream snapshots only the unidirectional peaks' names, and ignores the log
(its name has a timestamp); here the log's presence is checked too. Upstream's
failure case checks the process exits 1 when the reads carry no signal; here
the task's failure is caught.
"""

from flyte.io import File

from flyte_bio.modules.pints import PintsCallerResult, pints_caller
from tests.framework import assert_md5, env, fixture

GROHMM = "https://raw.githubusercontent.com/Kraus-Lab/groHMM/master/inst/extdata/"
ARGS = "--disable-small"


async def check(result: PintsCallerResult, divergent: list[str], bidirectional: list[str]) -> None:
    def names(files: list[File]) -> list[str]:
        return [f.name for f in files]

    n = len(divergent)
    expected = [f"test_{i}" for i in range(1, n + 1)]
    assert names(result.divergent) == [f"{e}_divergent_peaks.bed" for e in expected], names(result.divergent)
    assert names(result.bidirectional) == [f"{e}_bidirectional_peaks.bed" for e in expected]
    assert names(result.unidirectional) == [f"{e}_unidirectional_peaks.bed" for e in expected]
    for file, md5 in zip(result.divergent + result.bidirectional, divergent + bidirectional, strict=True):
        await assert_md5(file, md5, label=file.name)
    assert result.log.name.startswith("peakcalling_") and result.log.name.endswith(".log"), result.log.name


@env.task
async def test_pints_caller_two_bams() -> None:
    # upstream case: pints/caller "groHMM data - [bam, bam], fasta"
    bams = [await fixture(GROHMM + "S0mR1.bam"), await fixture(GROHMM + "S40mR1.bam")]
    result = await pints_caller(bams, "GROseq", prefix="test", args=ARGS)
    await check(
        result,
        ["d87cbbfa748e0a0ee01e26acb2922806", "a34c06a2fa7bd9490ef6a34b104651ce"],
        ["1c32c50def298e258b4ecd44f334cb9d", "18a1eac15b6ded3717523622948e9251"],
    )


@env.task
async def test_pints_caller_one_bam() -> None:
    # upstream case: pints/caller "groHMM data - [bam], fasta"
    result = await pints_caller([await fixture(GROHMM + "S0mR1.bam")], "GROseq", prefix="test", args=ARGS)
    await check(result, ["3aa8ae581240d4d26f63ea6f3ed1f233"], ["3c6ef0324e75082a6b172cb60549d2ad"])


@env.task
async def test_pints_caller_no_signal() -> None:
    # upstream case: pints/caller "Fails if no signal is found"
    bam = await fixture("genomics/sarscov2/illumina/bam/test.paired_end.sorted.bam")
    try:
        await pints_caller([bam], "GROseq", prefix="test", args=ARGS)
    except Exception:
        return
    raise AssertionError("PINTS called peaks on reads with no signal")


tests = [test_pints_caller_two_bams, test_pints_caller_one_bam, test_pints_caller_no_signal]
