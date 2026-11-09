"""Tests for flyte_bio.modules.find.

The cases and expected md5s are upstream's find/unpigz tests (modules commit
9a1cc9f, findutils 4.6.0, pigz 2.8). Upstream's failure case rejects files
without ``.gz`` before running anything; here the wrapper raises the same
message.
"""

from flyte_bio.modules.find import NOT_GZIPPED, find_unpigz
from tests.framework import assert_md5, env, fixture

GENOME = "genomics/sarscov2/genome/"
GFF3_MD5 = "357435a81a9981a0128e840ebe11051e"


@env.task
async def test_find_unpigz() -> None:
    # upstream case: find/unpigz "test_unpigz_success"
    files = await find_unpigz(
        [await fixture(GENOME + "genome.gff3.gz"), await fixture(GENOME + "alignment/last/contigs.genome.maf.gz")]
    )
    assert [f.name for f in files] == ["genome.gff3", "contigs.genome.maf"], [f.name for f in files]
    await assert_md5(files[0], GFF3_MD5, label="genome.gff3")
    await assert_md5(files[1], "caebfb2c47ffb50b1e64d5520dddd5d2", label="contigs.genome.maf")


@env.task
async def test_find_unpigz_single_file() -> None:
    # upstream case: find/unpigz "test_unpigz_success_single_file"
    files = await find_unpigz(await fixture(GENOME + "genome.gff3.gz"))
    assert [f.name for f in files] == ["genome.gff3"], [f.name for f in files]
    await assert_md5(files[0], GFF3_MD5, label="genome.gff3")


@env.task
async def test_find_unpigz_non_gzipped_files() -> None:
    # upstream case: find/unpigz "test_non_gzipped_files"
    files = [await fixture(GENOME + "genome.fasta"), await fixture(GENOME + "genome.sizes")]
    try:
        await find_unpigz(files)
    except ValueError as error:
        assert str(error) == NOT_GZIPPED, str(error)
    else:
        raise AssertionError("find_unpigz accepted files without .gz")


tests = [test_find_unpigz, test_find_unpigz_single_file, test_find_unpigz_non_gzipped_files]
