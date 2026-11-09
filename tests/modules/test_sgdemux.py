"""Tests for flyte_bio.modules.sgdemux.

The case and expected md5s are upstream's sgdemux test (modules commit
9a1cc9f, sgdemux 1.1.1) on simulated dual-indexed paired-end reads. As
upstream, the FASTQs (whose read order depends on the thread count) and
``most_frequent_unmatched.tsv`` are checked by name, and the empty
``sample_barcode_hop_metrics.tsv`` by name.

Upstream's test writes into the folder holding its input FASTQs, so its
snapshot also lists the two inputs ``out_L001_R1_001.fastq.gz`` and
``out_L001_R2_001.fastq.gz``; the module writes into its own folder, so they
aren't among its outputs.
"""

import os
import tarfile

from flyte_bio.modules.sgdemux import sgdemux
from tests.framework import assert_md5, env, fixture

SIM_DATA = "https://raw.githubusercontent.com/nf-core/test-datasets/demultiplex/testdata/sim-data/"
SAMPLES = ["Undetermined_S25"] + [f"s{n}_S{n}" for n in range(1, 25)]
FASTQS = sorted(f"{sample}_L001_R{read}_001.fastq.gz" for sample in SAMPLES for read in (1, 2))


def sgdemux_local_archive_names(path: str) -> list[str]:
    with tarfile.open(path) as tar:
        return sorted(m.name.removeprefix("./") for m in tar if m.isfile())


@env.task
async def test_sgdemux_sim_data() -> None:
    # upstream case: "test_sgdemux"
    sheet, fastqs = await fixture(SIM_DATA + "out.sample_meta.csv"), await fixture(SIM_DATA + "fastq.tar.gz")
    result = await sgdemux(sheet, fastqs, prefix="sim-data")
    names = sgdemux_local_archive_names(str(await result.fastq.download()))
    assert names == FASTQS, f"FASTQs {names}"
    await assert_md5(result.metrics, "7627e4bd2a56ea551fd74b4a2b5cb4b2", label="metrics.tsv")
    await assert_md5(result.per_project_metrics, "7393bc580c8ed8431107ad3dd00a7d2a", label="per_project_metrics.tsv")
    await assert_md5(result.per_sample_metrics, "fb8ce024ad2508e7b22a4af65fe44358", label="per_sample_metrics.tsv")
    for file, name in (
        (result.most_frequent_unmatched, "most_frequent_unmatched.tsv"),
        (result.sample_barcode_hop_metrics, "sample_barcode_hop_metrics.tsv"),
    ):
        assert os.path.basename(file.path) == name, file.path


tests = [test_sgdemux_sim_data]
