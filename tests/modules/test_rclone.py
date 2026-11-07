"""Tests for flyte_bio.modules.rclone.

The cases and expected md5s are upstream's rclone/check, rclone/checksum and
rclone/copy tests (modules commit 9a1cc9f, rclone 1.74.3), with the same
``args``. Upstream's check and checksum cases read public S3 folders; here
rclone only ever reads local files, staged into the task as the ``staged:``
remote, and never a cloud remote or credentials. The reports list paths, not
contents, so each folder is rebuilt with the files upstream's folders hold
(same paths; same or differing contents where the case needs it), and the
reports' md5s are upstream's. The checksum match case needs the real files
its sum file lists, so those are test fixtures (seven small files of the
public snpEff cache, fetched over HTTPS like any fixture).

Upstream's copy cases are dry runs of copying a test-data file to ``/tmp/``;
their log is unstable, so upstream snapshots only its name. Here the source
is the same file, staged, and the log is also checked for the dry-run copy,
and for a filtered dry run of the folder holding the file. One more case,
not upstream's, copies for real.
"""

import asyncio
import tempfile
from pathlib import Path

from flyte.io import Dir, File

from flyte_bio.modules.rclone import RcloneCheckResult, rclone_check, rclone_checksum, rclone_copy
from tests.framework import assert_md5, env, fixture

EXIT_0 = "897316929176464ebc9ad085f31e7284"  # "0\n"
EXIT_1 = "b026324c6904b2a9cb4b88d6d61c81d1"  # "1\n"
SNPEFF_CACHE = "https://annotation-cache.s3.eu-west-1.amazonaws.com/snpeff_cache/WBcel235.105/WBcel235.105/"
SNPEFF_SEQUENCES = [f"sequence.{chromosome}.bin" for chromosome in ("I", "II", "III", "IV", "V", "X")] + [
    "sequence.bin"
]
GVCF = "genomics/homo_sapiens/illumina/gvcf/test.genome.vcf.gz"

# Upstream's S3 folders, by path. The two README.md files differ.
MALT = {"megan-nucl-Feb2022.db.zip": "megan", "test.rma6": "rma6"}
SNPEFF = {f"WBcel235.105/{name}": name for name in [*SNPEFF_SEQUENCES, "snpEffectPredictor.bin"]}
MIRTOP = {"README.md": "mirtop", "hairpin_mirtop.fa": ">hairpin", "hsa.gff3": "gff", "sim_isomir_sort.bam": "bam"}
METABULI = {
    "README.md": "metabuli readme",
    "acc2taxid": "acc",
    "metabuli_db.tar.gz": "db",
    "taxonomy/merged.dmp": "",
    "taxonomy/names.dmp": "names",
    "taxonomy/nodes.dmp": "nodes",
}


async def staged_dir(folders: dict[str, dict[str, str]], files: dict[str, File] | None = None) -> Dir:
    """A Dir of ``folders`` (folder → path → text) and of ``files`` (path → File)."""
    root = Path(tempfile.mkdtemp())
    for folder, contents in folders.items():
        for path, content in contents.items():
            (root / folder / path).parent.mkdir(parents=True, exist_ok=True)
            (root / folder / path).write_text(content)
    for path, file in (files or {}).items():
        (root / path).parent.mkdir(parents=True, exist_ok=True)
        await file.download(str(root / path))
    return await Dir.from_local(str(root))


async def text_file(name: str, content: str) -> File:
    path = Path(tempfile.mkdtemp()) / name
    path.write_text(content)
    return await File.from_local(str(path))


async def check_reports(result: RcloneCheckResult, expected: dict[str, str], label: str) -> None:
    """``expected`` maps each report that should exist to its md5; the others must be absent."""
    for name in ("combined", "differ", "missing_on_dst", "missing_on_src", "match", "error", "exit_code"):
        report = getattr(result, name)
        if name not in expected:
            assert report is None, f"{label}: unexpected {report and report.name}"
            continue
        assert report is not None and report.name == f"test.{name}.txt", f"{label}: no test.{name}.txt"
        await assert_md5(report, expected[name], label=f"{label} {report.name}")


async def text(file: File) -> str:
    async with file.open("rb") as fh:
        return bytes(await fh.read()).decode()


@env.task
async def test_rclone_check_match() -> None:
    # upstream case: rclone/check "match - malt" (a folder against itself)
    staged = await staged_dir({"malt": MALT})
    result = await rclone_check("staged:malt/", "staged:malt/", "test", staged=staged, args="--no-check-certificate")
    expected = {
        "combined": "6518ec1170e94cd1e9d2ef4a4ba5a9e8",
        "match": "66c7833b0fd04a059b29a116474c021c",
        "exit_code": EXIT_0,
    }
    await check_reports(result, expected, "match")


@env.task
async def test_rclone_check_different_files() -> None:
    # upstream case: rclone/check "different files" (no paths in common)
    staged = await staged_dir({"malt": MALT, "snpeff": SNPEFF})
    result = await rclone_check("staged:malt/", "staged:snpeff/", "test", staged=staged, args="--no-check-certificate")
    expected = {
        "combined": "750454d95092af361bdb1c3e8c113251",
        "missing_on_dst": "66c7833b0fd04a059b29a116474c021c",
        "missing_on_src": "1dd18b615ae0c9abf940ad0739fb8950",
        "exit_code": EXIT_1,
    }
    await check_reports(result, expected, "different files")


@env.task
async def test_rclone_check_differ_and_missing() -> None:
    # upstream case: rclone/check "differ and missing" (README.md differs, the rest is on one side only)
    staged = await staged_dir({"mirtop": MIRTOP, "metabuli": METABULI})
    result = await rclone_check(
        "staged:mirtop/", "staged:metabuli", "test", staged=staged, args="--no-check-certificate"
    )
    expected = {
        "combined": "74c774156fdfa62f1ba14c780adf89f4",
        "differ": "ac036608c5c402f92926505d0f0ab0e4",
        "missing_on_dst": "6435d60effdac7ca6fbdfad40ffe2535",
        "missing_on_src": "89fcd615fd22b904204dafd407daf015",
        "exit_code": EXIT_1,
    }
    await check_reports(result, expected, "differ and missing")


@env.task
async def test_rclone_checksum_match() -> None:
    # upstream case: rclone/checksum "test - md5 - match"
    sumfile, *sequences = await asyncio.gather(
        fixture("generic/txt/snpeff_cache_md5.txt"), *(fixture(SNPEFF_CACHE + name) for name in SNPEFF_SEQUENCES)
    )
    files = {f"snpeff/WBcel235.105/{name}": file for name, file in zip(SNPEFF_SEQUENCES, sequences, strict=True)}
    staged = await staged_dir({"snpeff": {"WBcel235.105/snpEffectPredictor.bin": "predictor"}}, files)
    result = await rclone_checksum(
        sumfile, "MD5", "staged:snpeff/", "test", staged=staged, args="--no-check-certificate --one-way"
    )
    expected = {
        "combined": "4b849f0b6831d6561da709ca6f626b7d",
        "match": "2458aa282a426dec085359d209d5077c",
        "exit_code": EXIT_0,
    }
    await check_reports(result, expected, "md5 match")


@env.task
async def test_rclone_checksum_differs() -> None:
    # upstream case: rclone/checksum "test - md5 - differs" (a sum for a file the destination lacks)
    sumfile = await text_file("test_bad.md5", "00000000000000000000000000000000  hello.txt")
    staged = await staged_dir({"snpeff": SNPEFF})
    result = await rclone_checksum(
        sumfile, "MD5", "staged:snpeff/", "test", staged=staged, args="--no-check-certificate --one-way"
    )
    expected = {
        "combined": "d79adc0e541e01c95735769c919ae993",
        "missing_on_dst": "a438474115db37daf9d8e1307c06eb4a",
        "exit_code": EXIT_1,
    }
    await check_reports(result, expected, "md5 differs")


async def staged_gvcf() -> Dir:
    return await staged_dir({}, {"test.genome.vcf.gz": await fixture(GVCF)})


@env.task
async def test_rclone_copy_dry_run() -> None:
    # upstream case: rclone/copy "homo_sapiens - gvcf - copy from https - dry-run"
    log = await rclone_copy(
        "staged:test.genome.vcf.gz",
        "/tmp/",
        "test",
        staged=await staged_gvcf(),
        args="--dry-run --no-check-certificate",
    )
    assert log.name == "test-rclone-copy.log", f"log is {log.name}"
    content = await text(log)
    assert "test.genome.vcf.gz: Skipped copy as --dry-run is set" in content, f"dry-run log: {content}"


@env.task
async def test_rclone_copy_filter_dry_run() -> None:
    # upstream case: rclone/copy "homo_sapiens - gvcf - filter - dry-run". Upstream's source is the file
    # itself, which rclone copies whatever the filters say (so upstream checks only the log's name); the
    # same filter on the folder holding it does exclude it.
    exclude, staged = await asyncio.gather(text_file("exclude.txt", "test.genome.vcf.gz"), staged_gvcf())
    args = "--dry-run --no-check-certificate --exclude-from exclude.txt"
    file_log, folder_log = await asyncio.gather(
        rclone_copy("staged:test.genome.vcf.gz", "/tmp/", "test", filter_file=exclude, staged=staged, args=args),
        rclone_copy("staged:", "/tmp/", "test", filter_file=exclude, staged=staged, args=args),
    )
    assert file_log.name == "test-rclone-copy.log", f"log is {file_log.name}"
    content = await text(folder_log)
    assert "Skipped copy" not in content, f"the excluded file was copied: {content}"


@env.task
async def test_rclone_copy_local() -> None:
    # Not an upstream case: a real copy, between two local folders.
    log = await rclone_copy("staged:test.genome.vcf.gz", "/tmp/copied/", "test", staged=await staged_gvcf(), args="-v")
    content = await text(log)
    assert "test.genome.vcf.gz: Copied (new)" in content, f"copy log: {content}"


tests = [
    test_rclone_check_match,
    test_rclone_check_different_files,
    test_rclone_check_differ_and_missing,
    test_rclone_checksum_match,
    test_rclone_checksum_differs,
    test_rclone_copy_dry_run,
    test_rclone_copy_filter_dry_run,
    test_rclone_copy_local,
]
