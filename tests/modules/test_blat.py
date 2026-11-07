"""Tests for flyte_bio.modules.blat.

The case and expected md5 are upstream's blat test (modules commit 9a1cc9f,
BLAT 39x1). Upstream's setup step, ``seqtk seq -A`` (FASTQ to gzipped FASTA),
is done here in Python.
"""

import gzip
import tempfile
from pathlib import Path

from flyte.io import File

from flyte_bio.modules.blat import blat
from tests.framework import assert_md5, env, fixture


async def fastq_to_fasta_gz(fastq: File, name: str) -> File:
    """``seqtk seq -A``: each record's header line and sequence, as gzipped FASTA."""
    local = await fastq.download()
    out = Path(tempfile.mkdtemp()) / name
    with gzip.open(local, "rt") as fin, gzip.open(out, "wt") as fout:
        for i, line in enumerate(fin):
            if i % 4 == 0:
                fout.write(">" + line[1:])
            elif i % 4 == 1:
                fout.write(line)
    return await File.from_local(str(out))


@env.task
async def test_blat_sarscov2() -> None:
    reads = await fixture("genomics/sarscov2/illumina/fastq/test_1.fastq.gz")
    query = await fastq_to_fasta_gz(reads, "test.seqtk-seq.fasta.gz")
    psl = await blat(query, await fixture("genomics/sarscov2/genome/genome.fasta"), prefix="test")
    await assert_md5(psl, "6e2e5b3be48c84877f3c54b32bb9ec33", label="test.psl")


tests = [test_blat_sarscov2]
