"""Tests for flyte_bio.modules.crisprcleanr.

The case is upstream's crisprcleanr/normalize test (modules commit 9a1cc9f,
CRISPRcleanR 3.0.0). Upstream snapshots only the table's name and the
versions, and checks that the table names its sgRNAs.
"""

from flyte.io import File

from flyte_bio.modules.crisprcleanr import crisprcleanr_normalize
from tests.framework import env, fixture

TSV = "genomics/homo_sapiens/genome/tsv/"


async def text(file: File) -> str:
    async with file.open("rb") as fh:
        return bytes(await fh.read()).decode()


@env.task
async def test_crisprcleanr_normalize() -> None:
    counts = await fixture(TSV + "functional_genomics_counts.tsv")
    library = await fixture(TSV + "library_functional_genomics.tsv")
    result = await crisprcleanr_normalize(counts, library, min_reads=30, min_targeted_genes=5, sample="test")
    assert result.norm_count_file.path.endswith("test_norm_table.tsv"), result.norm_count_file.path
    assert "sgRNA" in await text(result.norm_count_file), "normalised table names no sgRNA"
    versions = await text(result.versions)
    for line in ('"CRISPRCLEANR_NORMALIZE":', "    r-base: 4.2.2", "    crisprcleanr: 3.0.0"):
        assert line in versions.splitlines(), f"versions.yml lacks {line!r}: {versions!r}"


tests = [test_crisprcleanr_normalize]
