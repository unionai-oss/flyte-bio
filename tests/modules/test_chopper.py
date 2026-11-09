"""Tests for flyte_bio.modules.chopper.

The cases are upstream's chopper tests (modules commit 9a1cc9f, chopper
0.14.1, pigz 2.8): Nanopore reads filtered with and without the SARS-CoV-2
genome as the contaminant reference. chopper's output order depends on its
threads, so, as upstream, the tests check the output's line count and whether
one read's header is present (a read that maps to the genome is gone).
Decompressed line by line, to keep the check memory-bounded.
"""

import gzip

from flyte.io import File

from flyte_bio.modules.chopper import chopper
from tests.framework import env, fixture

NANOPORE = "genomics/sarscov2/nanopore/fastq/"
CONTAMINANT_READ = (
    "@a52a642e-88d0-4584-babd-414ea84db484 runid=71c83ae0021f873e29b130c6562a4c27185f93b8 read=2768 ch=489 "
    "start_time=2021-08-11T12:07:39Z flow_cell_id=FAQ57606 protocol_group_id=210811_47CoV_SA sample_id=CS5 "
    "barcode=barcode04 barcode_alias=barcode04"
)
KEPT_READ = (
    "@2109d790-67ec-4fd1-8931-6c7e61908ff3 runid=97ca62ca093ff43533aa34c38a10b1d6325e7e7b read=52274 ch=243 "
    "start_time=2021-02-05T23:27:30Z flow_cell_id=FAP51364 protocol_group_id=data sample_id=RN20097 "
    "barcode=barcode01 barcode_alias=barcode01"
)


async def lines_of(fastq: File) -> tuple[int, set[str]]:
    """The gzipped FASTQ's line count and its header lines."""
    local = await fastq.download()
    count, headers = 0, set()
    with gzip.open(local, "rt") as fh:
        for i, line in enumerate(fh):
            count += 1
            if i % 4 == 0:
                headers.add(line.rstrip("\n"))
    return count, headers


@env.task
async def test_chopper_with_contaminants() -> None:
    out = await chopper(
        await fixture(NANOPORE + "test_2.fastq.gz"),
        await fixture("genomics/sarscov2/genome/genome.fasta.gz"),
        prefix="test_out.chopper",
    )
    assert out.name == "test_out.chopper.fastq.gz", out.name
    count, headers = await lines_of(out)
    assert count == 15984, f"{count} lines"
    assert CONTAMINANT_READ not in headers, "a read mapping to the contaminant reference was kept"


@env.task
async def test_chopper_without_contaminants() -> None:
    out = await chopper(await fixture(NANOPORE + "test.fastq.gz"), prefix="test_out.chopper")
    assert out.name == "test_out.chopper.fastq.gz", out.name
    count, headers = await lines_of(out)
    assert count == 400, f"{count} lines"
    assert KEPT_READ in headers, "the expected read is missing"


tests = [test_chopper_with_contaminants, test_chopper_without_contaminants]
