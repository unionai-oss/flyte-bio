"""Tests for flyte_bio.modules.bases2fastq.

The case and expected md5s are upstream's bases2fastq test (modules commit
9a1cc9f, bases2fastq 2.3.0) on its simulated run folder. FASTQs are compared
decompressed, as nf-test hashes them. As upstream, the QC and MultiQC
reports, ``RunStats.json`` and the sample's stats are checked by name only.
"""

import asyncio
import gzip
import hashlib
import os
import tarfile

from flyte_bio.modules.bases2fastq import bases2fastq
from tests.framework import assert_md5, env, fixture

SIM_DATA = "https://github.com/nf-core/test-datasets/raw/demultiplex/testdata/sim-data/"
MANIFEST = "https://raw.githubusercontent.com/nf-core/test-datasets/demultiplex/testdata/sim-data/RunManifest.csv"


def bases2fastq_local_archive_md5s(path: str) -> dict[str, str]:
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


@env.task
async def test_bases2fastq_sim_data() -> None:
    # upstream case: "sim-data.tar.gz - bcl"
    run_dir, manifest = await fixture(SIM_DATA + "sim-data.tar.gz"), await fixture(MANIFEST)
    result = await bases2fastq(run_dir, manifest, prefix="sim-data")
    await assert_md5(result.metrics, "1db1769b9b5917b8353e05bf8fd699c7", label="Metrics.csv")
    await assert_md5(result.run_manifest, "d584cdf00f577937a6cb691bb5331fe8", label="RunManifest.json")
    await assert_md5(result.unassigned, "11c1693830ce941b8cfb8d2431a59097", label="UnassignedSequences.csv")
    top = sorted(os.path.basename(f.path) for f in await result.results.list_files())
    for name in ("20230404-Bases2Fastq-Sim_QC.html", "multiqc_report.html", "RunStats.json"):
        assert name in top, f"no {name} in {top}"
    members = await asyncio.to_thread(bases2fastq_local_archive_md5s, str(await result.samples.download()))
    fastqs = {os.path.basename(p): md5 for p, md5 in members.items() if p.endswith(".fastq.gz")}
    assert fastqs == {
        "DefaultSample_R1.fastq.gz": "44ff0899ad4568a5cab7ccb7df9928e8",
        "DefaultSample_R2.fastq.gz": "f57acd0136bac5780d1592a5de7005ef",
    }, fastqs
    assert any(p.endswith("/DefaultSample_stats.json") for p in members), sorted(members)


tests = [test_bases2fastq_sim_data]
