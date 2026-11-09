"""Tests for flyte_bio.modules.mgikit.

The case and expected md5s are upstream's mgikit/demultiplex test (modules
commit 9a1cc9f, mgikit 2.1.0), with its config's ``args``
(``--writing-buffer-size 131072 --threads 1``) and its default prefix
``out-test-1`` (id ``test``, lane 1). FASTQs are compared decompressed, as
nf-test hashes them.

Upstream's fixture run folder is a ``.zip``, which its test unzips into
``fc01/``; the module takes a ``.tar.gz``, so the test repacks the zip's
files as ``fc01.tar.gz``, which the module extracts into ``fc01/``.
"""

import asyncio
import gzip
import hashlib
import os
import tarfile
import tempfile
import zipfile
from pathlib import Path

from flyte.io import File

from flyte_bio.modules.mgikit import mgikit_demultiplex
from tests.framework import assert_gunzipped_md5, assert_md5, env, fixture

MGI = "genomics/homo_sapiens/mgi/"
FASTQS = {
    "Sample01_S1_L01_R1_001.fastq.gz": "19752026de713f514637fab29e4ffac4",
    "Sample01_S1_L01_R2_001.fastq.gz": "f9015fa4a73977dd093b8b2e4a67a316",
    "Sample02_S2_L01_R1_001.fastq.gz": "7f942bf3ef88ef04431ed0d40baae144",
    "Sample02_S2_L01_R2_001.fastq.gz": "7348f1bf066a69ee48ed3ee6afa2f91b",
    "Sample03_S3_L01_R1_001.fastq.gz": "56afcab8b820419bddec4982fc85b60f",
    "Sample03_S3_L01_R2_001.fastq.gz": "38ba38d0d71343e6125fc15b63483ad4",
    "Sample04_S4_L01_R1_001.fastq.gz": "dcf0ad68be6ff6c21ebceedaa4eaa4c4",
    "Sample04_S4_L01_R2_001.fastq.gz": "cd88081007614f5bed7975b04e25b7e3",
}
UNDETERMINED = {
    "Undetermined_L01_R1_001.fastq.gz": "63ce3d442f4f8c42e0d9addea9273315",
    "Undetermined_L01_R2_001.fastq.gz": "0af27de6fc08cce579e160c31fbe61b9",
}
REPORTS = {
    "FC01.L01.mgikit.general": "624b0d43c3995fde2c122c447e780191",
    "FC01.L01.mgikit.info": "69a72a1c43d47032828d8b97c6cf8807",
    "FC01.L01.mgikit.sample_stats": "20d5b9f4a430fd34579e0f70b42f73f0",
    "FC01.L01.mgikit.undetermined_barcode": "c1960b64d4cc5c141d742aa1b6f57a31",
    "FC01.L01.mgikit.undetermined_barcode.complete": "c1960b64d4cc5c141d742aa1b6f57a31",
}


def mgikit_local_archive_md5s(path: str) -> dict[str, str]:
    """Each regular file's md5 in the ``.tar.gz`` at ``path`` (``.gz`` members decompressed), by member path."""
    md5s = {}
    with tarfile.open(path) as tar:
        for member in tar:
            if not member.isfile():
                continue
            raw = tar.extractfile(member)
            assert raw is not None
            stream = gzip.GzipFile(fileobj=raw) if member.name.endswith(".gz") else raw
            h = hashlib.md5()
            for chunk in iter(lambda: stream.read(1 << 20), b""):
                h.update(chunk)
            md5s[member.name.removeprefix("./")] = h.hexdigest()
    return md5s


def mgikit_repack_zip(zip_path: str, name: str) -> str:
    """The files of the zip at ``zip_path`` repacked as ``<name>.tar.gz``, at its top level."""
    root = Path(tempfile.mkdtemp())
    with zipfile.ZipFile(zip_path) as zf:
        zf.extractall(root / "files")
    archive = root / f"{name}.tar.gz"
    with tarfile.open(archive, "w:gz") as tar:
        for path in sorted((root / "files").iterdir()):
            tar.add(path, arcname=path.name)
    return str(archive)


async def mgikit_run_folder() -> File:
    zipped = await fixture(MGI + "fc01.zip")
    archive = await asyncio.to_thread(mgikit_repack_zip, str(await zipped.download()), "fc01")
    return await File.from_local(archive)


@env.task
async def test_mgikit_demultiplex() -> None:
    # upstream case: "run mgikit demultiplex without errors"
    sheet, run_dir = await fixture(MGI + "fc01_sample_sheet.csv"), await mgikit_run_folder()
    result = await mgikit_demultiplex(
        sheet, run_dir, prefix="out-test-1", args="--writing-buffer-size 131072 --threads 1"
    )
    fastqs = await asyncio.to_thread(mgikit_local_archive_md5s, str(await result.fastq.download()))
    assert fastqs == FASTQS, f"sample FASTQs {fastqs}"
    assert result.ambiguous == [], result.ambiguous
    undetermined = {os.path.basename(f.path): f for f in result.undetermined}
    assert sorted(undetermined) == sorted(UNDETERMINED), sorted(undetermined)
    for name, md5 in UNDETERMINED.items():
        await assert_gunzipped_md5(undetermined[name], md5, label=name)
    for name, md5 in REPORTS.items():
        report = await result.results.get_file(name)
        assert report is not None, f"no {name}"
        await assert_md5(report, md5, label=name)
    top = [os.path.basename(f.path) for f in await result.results.list_files()]
    assert not any("ambiguous_barcode" in name for name in top), top


tests = [test_mgikit_demultiplex]
