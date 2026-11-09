"""Tests for flyte_bio.modules.fqtk.

The case and expected md5s are upstream's fqtk test (modules commit 9a1cc9f,
fqtk 0.4.0) on simulated dual-indexed paired-end reads. FASTQs are compared
decompressed, as nf-test hashes them.
"""

import asyncio
import gzip
import hashlib
import tarfile

from flyte_bio.modules.fqtk import fqtk
from tests.framework import assert_md5, env, fixture

SIM_DATA = "https://github.com/nf-core/test-datasets/raw/demultiplex/testdata/sim-data/"
READ_STRUCTURES = [
    ("out_L001_I1_001.fastq.gz", "8B"),
    ("out_L001_I2_001.fastq.gz", "8B"),
    ("out_L001_R1_001.fastq.gz", "150T"),
    ("out_L001_R2_001.fastq.gz", "150T"),
]
FASTQS = {
    "s1.R1.fq.gz": "1eadc3622432e260b2d1c43d669b621a",
    "s1.R2.fq.gz": "fe11bc516e32e2ee6c3064dd7d157bc8",
    "s10.R1.fq.gz": "68d72892789df70b17b8886d6a616999",
    "s10.R2.fq.gz": "abe44a80eebd1a4a7ec432ee3142f599",
    "s11.R1.fq.gz": "1e63515691b7694b3237f69e95dacdd7",
    "s11.R2.fq.gz": "3431b7a6d5840afb4460529a6f9611e3",
    "s12.R1.fq.gz": "a27e66093a904a269bfcc3c3f98a1238",
    "s12.R2.fq.gz": "25b5b8ced9353f202880077aabd2a721",
    "s13.R1.fq.gz": "353b2cbd1394d0e1c03cd3a4a7f4f623",
    "s13.R2.fq.gz": "3e3c6ed5933d864759cceee197f8fa64",
    "s14.R1.fq.gz": "83744d94d2f5d3722c960566b7a5267e",
    "s14.R2.fq.gz": "1aecf34deb4bc27069f4ff83a0fb9f5d",
    "s15.R1.fq.gz": "ea1f11b52fc604e4aa168a78bf086892",
    "s15.R2.fq.gz": "f059ac2b1125d364570be8a5010e94d8",
    "s16.R1.fq.gz": "150db1d3462b6cdd8022011b462d190d",
    "s16.R2.fq.gz": "a939a31d634df5a6c5741a79ad84fea5",
    "s17.R1.fq.gz": "eb2a43e737d5bc941866957d878f52cb",
    "s17.R2.fq.gz": "affeeeb95f4ab3eb73ff3fd49d62a460",
    "s18.R1.fq.gz": "413c245496557c3973bcb2c8fd738fed",
    "s18.R2.fq.gz": "50d854cd559fb641967e8eda19006f8f",
    "s19.R1.fq.gz": "9581c36425dd36b1bf4a844d4c3fd434",
    "s19.R2.fq.gz": "13e966cdfcf3fd272a1d739dbe82e05c",
    "s2.R1.fq.gz": "e5a9598bef0d53087aa25530667ccf79",
    "s2.R2.fq.gz": "7cbede4e3e641207da50ebb96653d0df",
    "s20.R1.fq.gz": "32144c46d24d5c264d98be49b8d8fb1d",
    "s20.R2.fq.gz": "585a1ff0e75a693c7464f6d0966c5fb9",
    "s21.R1.fq.gz": "1e92094589206300ef884273ee9ee265",
    "s21.R2.fq.gz": "9f9278e60ba2dd2c5bb14d021268c613",
    "s22.R1.fq.gz": "5abfd96f6a1eb0391208fe7f4dcecade",
    "s22.R2.fq.gz": "868ef4eb6ca676ed93d2389188ed3d81",
    "s23.R1.fq.gz": "0d0ec11e03749ac7a48d00a46781a74f",
    "s23.R2.fq.gz": "3ca8a6141f9e4f9595a90075ac390e7f",
    "s24.R1.fq.gz": "9af8ba800a69fc8637cfe864870bade7",
    "s24.R2.fq.gz": "1e988f6a217b67271801f2041be5c66b",
    "s3.R1.fq.gz": "aa1fe4107eb39313f97c1c045a3a15ed",
    "s3.R2.fq.gz": "22332ab06eed7ff915ec1ffba0100610",
    "s4.R1.fq.gz": "f6818e27d3118f66168b8af2f155f4b0",
    "s4.R2.fq.gz": "c489acd3fd4775cb9a69053d4b26250f",
    "s5.R1.fq.gz": "2c111ddbdc64a9f4a9fc7a27f2e2aaea",
    "s5.R2.fq.gz": "6ae8e46cf5fa66f6998a72f15caa9c58",
    "s6.R1.fq.gz": "cc6c0fb00166827307c846c8e9661c47",
    "s6.R2.fq.gz": "6d5556ca178aab7202cf860fff94f574",
    "s7.R1.fq.gz": "36ab856e02e7d1c17c3a122ec9179cbb",
    "s7.R2.fq.gz": "0764aae83716030f4ea05cca5d4511a4",
    "s8.R1.fq.gz": "f37f97582eafb96aae79a1d635e4585e",
    "s8.R2.fq.gz": "b5e48a29663e3e92354794290e2ba419",
    "s9.R1.fq.gz": "42e79f8a4e8db043b76782538982d1ea",
    "s9.R2.fq.gz": "d3af5fce6432934a479770f16c8f50b1",
    "unmatched.R1.fq.gz": "8f77e8d5adc4b16dd499fe3bd092f66d",
    "unmatched.R2.fq.gz": "a2c623dcaa11fe3faaecc563c376f24b",
}


def fqtk_local_archive_md5s(path: str) -> dict[str, str]:
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
async def test_fqtk_sim_data() -> None:
    # upstream case: "tsv - fastq"
    sheet = await fixture(SIM_DATA + "fqtk_sample_metadata_subset.tsv")
    fastqs = await fixture(SIM_DATA + "fastq.tar.gz")
    result = await fqtk(sheet, fastqs, READ_STRUCTURES)
    await assert_md5(result.metrics, "1d587fa959f9129155314cf531103347", label="demux-metrics.txt")
    found = await asyncio.to_thread(fqtk_local_archive_md5s, str(await result.fastq.download()))
    assert found == FASTQS, f"FASTQs differ: {sorted(set(found.items()) ^ set(FASTQS.items()))}"


tests = [test_fqtk_sim_data]
