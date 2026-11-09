"""Tests for flyte_bio.modules.bcl2fastq.

The case and expected md5s are upstream's bcl2fastq test (modules commit
9a1cc9f, bcl2fastq 2.20.0.422), with its config's ``--tiles s_1_1101``.
FASTQs are compared decompressed, as nf-test hashes them.

bcl2fastq's statistics depend on its thread count: with 4 or more processing
threads its unknown-barcode counts change, and with them ``Stats.json``,
``ConversionStats.xml``, ``DemuxSummaryF1L1.txt`` and three report pages.
Upstream's snapshot matches 1 to 3 threads, so the test pins 2 CPUs, salting
the cache (which ignores resources).

Upstream snapshots the ``Reports/`` pages as a nested list without their
folders; the test compares the same (name, md5) pairs. Upstream's ``xml``
output also picks up the run folder's ``RTAConfiguration.xml`` and
``config.xml`` (inputs); the module returns only its ``RunInfo.xml``.
"""

import asyncio
import gzip
import hashlib
import os
import tarfile

import flyte
from flyte.io import Dir, File

from flyte_bio.modules.bcl2fastq import bcl2fastq_cmd
from tests.framework import assert_md5, env, fixture

BCL = "genomics/homo_sapiens/illumina/bcl/"
TWO_CPUS = {"resources": flyte.Resources(cpu=2, memory="8Gi"), "cache": flyte.Cache(behavior="auto", salt="cpu=2")}

REPORTS = [
    ("lane.html", "794e48287f47a9f22dcb6b6d0c22c3eb"),
    ("laneBarcode.html", "2bbdae3ee57151eab520c966597d7438"),
    ("lane.html", "f741870307050dcea79a838f5971770f"),
    ("laneBarcode.html", "ffe2e863811c76cb9da27d5d124e0611"),
    ("lane.html", "cca97e771d3491dbc8cd3fe389595b09"),
    ("laneBarcode.html", "cca97e771d3491dbc8cd3fe389595b09"),
    ("lane.html", "c383b0768d9978733d3f5d3b91c15af0"),
    ("laneBarcode.html", "48842c23b9a2816aec540177df870968"),
    ("Report.css", "eb7d3eb68fc1539f411404987246b59b"),
    ("index.html", "5747c407854ae2c358d0ec201ce622d8"),
    ("tree.html", "a1b9bf592973ca829ec69ddf888b7e34"),
]
STATS = {
    "AdapterTrimming.txt": "48ed2b914b1246c0b5d8667525550946",
    "ConversionStats.xml": "8fe0f57f3f5d256a0762dba75ac62d05",
    "DemultiplexingStats.xml": "2047ff18f5b9107c084de06e9ff943ad",
    "DemuxSummaryF1L1.txt": "03e5fd0c1e3079c5f8c7b4d0501b37ff",
    "FastqSummaryF1L1.txt": "0c6f2d87ee183b84d1051cde9a5643d1",
    "Stats.json": "8e5f038b8aa9e465599d3575f930e604",
}
INTEROP = {
    "ControlMetricsOut.bin": "6d77b38d0793a6e1ce1e85706e488953",
    "CorrectedIntMetricsOut.bin": "2bbf84d3be72734addaa2fe794711434",
    "ErrorMetricsOut.bin": "38c88def138e9bb832539911affdb286",
    "ExtractionMetricsOut.bin": "7497c3178837eea8f09350b5cd252e99",
    "IndexMetricsOut.bin": "9e688c58a5487b8eaf69c9e1005ad0bf",
    "QMetricsOut.bin": "7e9f198d53ebdfbb699a5f94cf1ed51c",
    "TileMetricsOut.bin": "83891751ec1c91a425a524b476b6ca3c",
}


def bcl2fastq_local_archive_md5s(path: str) -> dict[str, str]:
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


async def bcl2fastq_archive_md5s(archive: File) -> dict[str, str]:
    return await asyncio.to_thread(bcl2fastq_local_archive_md5s, str(await archive.download()))


async def bcl2fastq_check_files(results: Dir, folder: str, expected: dict[str, str]) -> None:
    for name, md5 in expected.items():
        found = await results.get_file(f"{folder}/{name}")
        assert found is not None, f"no {folder}/{name}"
        await assert_md5(found, md5, label=f"{folder}/{name}")


@env.task
async def test_bcl2fastq_flowcell() -> None:
    # upstream case: "homo sapiens illumina [bcl]"
    sheet, run_dir = await fixture(BCL + "flowcell_samplesheet.csv"), await fixture(BCL + "flowcell.tar.gz")
    task = bcl2fastq_cmd.as_task().override(**TWO_CPUS)
    results = await task(samplesheet=sheet, run_dir=run_dir, args="--tiles s_1_1101", args2="", args3="")
    fastq, reports = await asyncio.gather(results.get_file("fastq.tar.gz"), results.get_file("Reports.tar.gz"))
    assert fastq is not None and reports is not None, "bcl2fastq wrote no fastq.tar.gz or Reports.tar.gz"
    fastqs = await bcl2fastq_archive_md5s(fastq)
    assert fastqs == {"Sample1_S1_L001_R1_001.fastq.gz": "0675fb6365322eaafb33c0f8e862b54b"}, fastqs
    pages = sorted((os.path.basename(path), md5) for path, md5 in (await bcl2fastq_archive_md5s(reports)).items())
    assert pages == sorted(REPORTS), f"Reports pages {pages}"
    await bcl2fastq_check_files(results, "Stats", STATS)
    await bcl2fastq_check_files(results, "InterOp", INTEROP)
    run_info = await results.get_file("RunInfo.xml")
    assert run_info is not None, "no RunInfo.xml"
    await assert_md5(run_info, "03038959f4dd181c86bc97ae71fe270a", label="RunInfo.xml")


tests = [test_bcl2fastq_flowcell]
