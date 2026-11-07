"""Tests for flyte_bio.modules.seqtk.

The seqtk_sample cases and expected md5s are upstream's seqtk/sample tests (modules commit
9a1cc9f, seqtk 1.4), with their config's ``ext.args = '-s100'`` and
``ext.prefix = "<id>.sampled"``. Upstream's "paired-end" case samples the same
single file as the single-end case, so it is kept as it is. nf-test hashes
``.gz`` outputs decompressed.

The seqtk_trim cases are upstream's seqtk/trim tests (same commit and
container): single- and paired-end, with the default prefix (the sample id,
``test``) and no ``args``.
"""

from flyte.io import File

from flyte_bio.modules.seqtk import seqtk_sample, seqtk_trim
from tests.framework import assert_gunzipped_md5, env, fixture

READS = "genomics/sarscov2/illumina/fastq/test_1.fastq.gz"
READS_2 = "genomics/sarscov2/illumina/fastq/test_2.fastq.gz"


def name_of(file: File) -> str:
    return file.path.rstrip("/").rsplit("/", 1)[-1]


async def check(sample_size: int | float, md5: str) -> None:
    outs = await seqtk_sample([await fixture(READS)], sample_size, prefix="test.sampled", args="-s100")
    assert [name_of(o) for o in outs] == ["test.sampled_test_1.fastq.gz"], [name_of(o) for o in outs]
    await assert_gunzipped_md5(outs[0], md5, label="test.sampled_test_1.fastq.gz")


@env.task
async def test_seqtk_sample_single_end() -> None:
    await check(50, "e5f44fafd7351c5abb9925a075132941")


@env.task
async def test_seqtk_sample_paired_end() -> None:
    await check(50, "e5f44fafd7351c5abb9925a075132941")


@env.task
async def test_seqtk_sample_single_end_fraction() -> None:
    await check(0.5, "7552511dcd914f5183d80fae459bf917")


@env.task
async def test_seqtk_trim_single_end() -> None:
    outs = await seqtk_trim([await fixture(READS)], prefix="test")
    assert [name_of(o) for o in outs] == ["test_test_1.fastq.gz"], [name_of(o) for o in outs]
    await assert_gunzipped_md5(outs[0], "4161df271f9bfcd25d5845a1e220dbec", label="test_test_1.fastq.gz")


@env.task
async def test_seqtk_trim_paired_end() -> None:
    outs = await seqtk_trim([await fixture(READS), await fixture(READS_2)], prefix="test")
    names = [name_of(o) for o in outs]
    assert names == ["test_test_1.fastq.gz", "test_test_2.fastq.gz"], names
    await assert_gunzipped_md5(outs[0], "4161df271f9bfcd25d5845a1e220dbec", label=names[0])
    await assert_gunzipped_md5(outs[1], "2ebae722295ea66d84075a3b042e2b42", label=names[1])


tests = [
    test_seqtk_sample_single_end,
    test_seqtk_sample_paired_end,
    test_seqtk_sample_single_end_fraction,
    test_seqtk_trim_single_end,
    test_seqtk_trim_paired_end,
]
