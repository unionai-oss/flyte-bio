"""Tests for flyte_bio.modules.gtfgenefilter.

Upstream has no module tests for this local script, so the test checks the
filtering itself: lines on the genome's sequences are kept as they are
(including gene records without a transcript_id), others dropped.
"""

import tempfile
from pathlib import Path

from flyte.io import File

from flyte_bio.modules.gtfgenefilter import gtf_gene_filter
from tests.framework import env, fixture


@env.task
async def test_gtf_gene_filter() -> None:
    fasta = await fixture("genomics/sarscov2/genome/genome.fasta")
    lines = [
        'MT192765.1\tsrc\tgene\t1\t100\t.\t+\t.\tgene_id "g1";\n',
        'MT192765.1\tsrc\texon\t1\t100\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n',
        'chrUn\tsrc\texon\t1\t100\t.\t+\t.\tgene_id "g2"; transcript_id "t2";\n',
    ]
    path = Path(tempfile.mkdtemp()) / "genes.gtf"
    path.write_text("".join(lines))
    filtered = await gtf_gene_filter(fasta, await File.from_local(str(path)))
    assert filtered.path.endswith("/genome_genes.gtf"), filtered.path
    async with filtered.open("rb") as fh:
        kept = bytes(await fh.read()).decode().splitlines(keepends=True)
    assert kept == lines[:2], kept


tests = [test_gtf_gene_filter]
