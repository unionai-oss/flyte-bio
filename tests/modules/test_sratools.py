"""Tests for flyte_bio.modules.sratools.

The cases and expected md5s are upstream's sratools/prefetch,
sratools/fasterqdump and custom/sratoolsncbisettings tests (modules commit
9a1cc9f, sra-tools 3.2.1). They download from NCBI, so they need network.
nf-test hashes the dumped reads decompressed.

Upstream tests the settings task with existing settings
(``NCBI_SETTINGS=/tmp/.ncbi/user-settings.mkfg``, a copy of the test-data
settings), so that case passes the same file as ``existing_settings``.
``test_sratools_ncbi_settings_generated`` also covers the branch upstream's
tests skip, where no settings exist and the task writes them.
"""

import tempfile
from pathlib import Path

from flyte.io import File

from flyte_bio.modules.sratools import sratools_fasterqdump, sratools_ncbi_settings, sratools_prefetch
from tests.framework import assert_gunzipped_md5, assert_md5, env, fixture, fixture_dir

SETTINGS = "generic/config/ncbi_user_settings.mkfg"
SRA = "genomics/sarscov2/illumina/sra/"


async def text(file: File) -> str:
    async with file.open("rb") as fh:
        return bytes(await fh.read()).decode()


def check_settings(content: str) -> None:
    assert any("/LIBS/GUID" in line for line in content.splitlines()), "settings lack /LIBS/GUID"
    assert any("/libs/cloud/report_instance_identity" in line for line in content.splitlines()), (
        "settings lack /libs/cloud/report_instance_identity"
    )


@env.task
async def test_sratools_ncbi_settings_existing() -> None:
    content = await text(await fixture(SETTINGS))
    local = Path(tempfile.mkdtemp()) / "user-settings.mkfg"
    local.write_text(content)
    settings = await sratools_ncbi_settings(existing_settings=await File.from_local(str(local)))
    assert settings.name == "user-settings.mkfg", f"settings named {settings.name!r}"
    check_settings(await text(settings))


@env.task
async def test_sratools_ncbi_settings_generated() -> None:
    settings = await sratools_ncbi_settings()
    assert settings.name == "user-settings.mkfg", f"settings named {settings.name!r}"
    check_settings(await text(settings))


@env.task
async def test_sratools_prefetch() -> None:
    sra = await sratools_prefetch("DRR000774", await fixture(SETTINGS))
    files = await sra.list_files()
    assert [f.name for f in files] == ["DRR000774.sra"], f"prefetch wrote {[f.name for f in files]}"
    await assert_md5(files[0], "7647dba20c89c0e3d7ad13842f060eb0", label="DRR000774.sra")


@env.task
async def test_sratools_prefetch_sralite() -> None:
    sra = await sratools_prefetch("SRR1170046", await fixture(SETTINGS))
    files = await sra.list_files()
    assert [f.name for f in files] == ["SRR1170046.sralite"], f"prefetch wrote {[f.name for f in files]}"
    await assert_md5(files[0], "7acfce556ca0951aff49d780899c105b", label="SRR1170046.sralite")


async def check_reads(reads: list[File], expected: dict[str, str]) -> None:
    names = [r.name for r in reads]
    assert names == list(expected), f"fasterq-dump wrote {names}, expected {list(expected)}"
    for read in reads:
        await assert_gunzipped_md5(read, expected[read.name], label=read.name)


@env.task
async def test_sratools_fasterqdump_single_end() -> None:
    sra = await fixture_dir(SRA + "SRR13255544.tar.gz", "SRR13255544")
    reads = await sratools_fasterqdump(sra, await fixture(SETTINGS), prefix="test_single_end", single_end=True)
    await check_reads(reads, {"test_single_end.fastq.gz": "674d78c1cc3c1308d6d39d6369a42887"})


@env.task
async def test_sratools_fasterqdump_paired_end() -> None:
    sra = await fixture_dir(SRA + "SRR11140744.tar.gz", "SRR11140744")
    reads = await sratools_fasterqdump(sra, await fixture(SETTINGS), prefix="test_paired_end", single_end=False)
    await check_reads(
        reads,
        {
            "test_paired_end_1.fastq.gz": "8573015c91d099b6e30789f8bab2f43c",
            "test_paired_end_2.fastq.gz": "37e6f719a022dc3c9994c80fbc20c311",
        },
    )


tests = [
    test_sratools_ncbi_settings_existing,
    test_sratools_ncbi_settings_generated,
    test_sratools_prefetch,
    test_sratools_prefetch_sralite,
    test_sratools_fasterqdump_single_end,
    test_sratools_fasterqdump_paired_end,
]
