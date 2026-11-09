"""Tests for flyte_bio.modules.bclconvert.

The cases and expected md5s are upstream's bclconvert tests (modules commit
9a1cc9f, BCL Convert 4.6.2), with their configs' ``args``. FASTQs are
compared decompressed, as nf-test hashes them. As upstream, the
undetermined reads are checked by name only, and the SARS-CoV-2 case leaves
out ``SampleZ``, whose barcode matches no reads.
"""

import asyncio
import gzip
import hashlib
import os
import tarfile

from flyte.io import Dir

from flyte_bio.modules.bclconvert import BclconvertResult, bclconvert
from tests.framework import env, fixture

UNDETERMINED = "Undetermined_S0_L001_R1_001.fastq.gz"
HUMAN_INTEROP = [
    "ControlMetricsOut.bin",
    "CorrectedIntMetricsOut.bin",
    "ErrorMetricsOut.bin",
    "ExtractionMetricsOut.bin",
    "IndexMetricsOut.bin",
    "QMetricsOut.bin",
    "TileMetricsOut.bin",
]
SARSCOV2_INTEROP = [
    "BasecallingMetricsOut.bin",
    "CorrectedIntMetricsOut.bin",
    "EmpiricalPhasingMetricsOut.bin",
    "EventMetricsOut.bin",
    "ExtendedTileMetricsOut.bin",
    "ExtractionMetricsOut.bin",
    "FWHMGridMetricsOut.bin",
    "ImageMetricsOut.bin",
    "IndexMetricsOut.bin",
    "OpticalModelMetricsOut.bin",
    "PFGridMetricsOut.bin",
    "QMetrics2030Out.bin",
    "QMetricsByLaneOut.bin",
    "QMetricsOut.bin",
    "RegistrationMetricsOut.bin",
    "TileMetricsOut.bin",
]


def bclconvert_local_archive_md5s(path: str) -> dict[str, str]:
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


async def bclconvert_names(results: Dir, folder: str) -> list[str]:
    files = [f async for f in results.walk()]
    prefix = f"/{folder}/"
    return sorted(os.path.basename(f.path) for f in files if prefix in f.path)


async def check_bclconvert(
    result: BclconvertResult, samples: dict[str, str], interop: list[str], skip: tuple[str, ...] = ()
) -> None:
    fastqs = await asyncio.to_thread(bclconvert_local_archive_md5s, str(await result.fastq.download()))
    assert UNDETERMINED in fastqs, f"no {UNDETERMINED}: {sorted(fastqs)}"
    found = {name: md5 for name, md5 in fastqs.items() if name != UNDETERMINED and name not in skip}
    assert found == samples, f"sample FASTQs {found}"
    assert not any("_I1_" in name or "_I2_" in name for name in fastqs), f"index FASTQs {sorted(fastqs)}"
    assert await bclconvert_names(result.results, "InterOp") == interop
    assert await bclconvert_names(result.results, "Reports"), "no Reports/"
    assert await bclconvert_names(result.results, "Logs"), "no Logs/"


@env.task
async def test_bclconvert_flowcell() -> None:
    # upstream case: "homo_sapiens illumina [bcl]"
    bcl = "genomics/homo_sapiens/illumina/bcl/"
    sheet, run_dir = await fixture(bcl + "flowcell_samplesheet.csv"), await fixture(bcl + "flowcell.tar.gz")
    result = await bclconvert(sheet, run_dir, args="--force --first-tile-only true")
    await check_bclconvert(
        result, {"Sample1_S1_L001_R1_001.fastq.gz": "0675fb6365322eaafb33c0f8e862b54b"}, HUMAN_INTEROP
    )


@env.task
async def test_bclconvert_sarscov2() -> None:
    # upstream case: "sars_cov2 illumina [bcl]"
    bcl = "genomics/sarscov2/illumina/bcl/"
    sheet = await fixture(bcl + "SampleSheet.csv")
    run_dir = await fixture(bcl + "200624_A00834_0183_BHMTFYDRXX.tar.gz")
    result = await bclconvert(sheet, run_dir, args="--force --tiles s_1_2101")
    samples = {
        "Sample1_S1_L001_R1_001.fastq.gz": "b5489d1964db8db5502eb742cc3ef3ec",
        "Sample23_S3_L001_R1_001.fastq.gz": "767a1091320320b140288066e29bccc5",
        "SampleA_S2_L001_R1_001.fastq.gz": "7de2ea88133409f34563f40a0d8c9e55",
        "sampletest_S4_L001_R1_001.fastq.gz": "c16c7de1b7bffb5e4503f4d94c40f881",
    }
    await check_bclconvert(result, samples, SARSCOV2_INTEROP, skip=("SampleZ_S5_L001_R1_001.fastq.gz",))


tests = [test_bclconvert_flowcell, test_bclconvert_sarscov2]
