"""Tests for flyte_bio.modules.prokka.

The cases and expected md5s are upstream's prokka tests (modules commit
9a1cc9f, prokka 1.15.6). As upstream, the ``gbk``, ``log`` and ``sqn`` files,
which record the date and paths, are only checked to exist.
"""

from flyte.io import File

from flyte_bio.modules.prokka import ProkkaResult, prokka
from tests.framework import assert_md5, assert_nonempty, env, fixture

GENOME = "genomics/sarscov2/genome/"

EXPECTED_MD5S = {
    "err": "f01faf5b7ce06f777eb0b39124c69ca6",
    "faa": "a01ec93f915874c4e932b6b6a115ef22",
    "ffn": "b68515be1070f9e3c9de86866151e784",
    "fna": "787307f29a263e5657cc276ebbf7e2b3",
    "fsa": "71bbefcb7f12046bcd3263f58cfd5404",
    "gff": "df0955fc26abf54635720408300e9a3d",
    "tbl": "1cb043abd014f25129a43c45a72eb73d",
    "tsv": "ee0dc43ee778e6d6274cffe24fcbc95f",
    "txt": "b40e485ffc8eaf1feacf8d79d9751a33",
}


async def check_outputs(result: ProkkaResult) -> None:
    for ext, md5 in EXPECTED_MD5S.items():
        file: File = getattr(result, ext)
        await assert_md5(file, md5, label=f"test.{ext}")
    for ext in ("gbk", "log", "sqn"):
        await assert_nonempty(getattr(result, ext), label=f"test.{ext}")


@env.task
async def test_prokka_genome_fasta() -> None:
    result = await prokka(await fixture(GENOME + "genome.fasta"), prefix="test")
    await check_outputs(result)


@env.task
async def test_prokka_genome_fasta_gz() -> None:
    result = await prokka(await fixture(GENOME + "genome.fasta.gz"), prefix="test")
    await check_outputs(result)


tests = [
    test_prokka_genome_fasta,
    test_prokka_genome_fasta_gz,
]
