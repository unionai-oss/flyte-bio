"""Tests for flyte_bio.modules.svim.

The case and expected lines are upstream's svim/alignment test (modules
commit 9a1cc9f, SVIM 2.0.0), with its params' ``--min_mapq 5``. SVIM's VCF
records the date on line 2, so as upstream the test compares line 1 and
lines 3-21.
"""

import asyncio

from flyte_bio.modules.svim import svim_alignment
from tests.framework import env, fixture

HUMAN = "genomics/homo_sapiens/"
HEADER = [
    "##source=SVIM-v2.0.0",
    "##contig=<ID=chr22,length=40001>",
    '##ALT=<ID=DEL,Description="Deletion">',
    '##ALT=<ID=INV,Description="Inversion">',
    '##ALT=<ID=DUP,Description="Duplication">',
    '##ALT=<ID=DUP:TANDEM,Description="Tandem Duplication">',
    '##ALT=<ID=DUP:INT,Description="Interspersed Duplication">',
    '##ALT=<ID=INS,Description="Insertion">',
    '##ALT=<ID=BND,Description="Breakend">',
    '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">',
    (
        "##INFO=<ID=CUTPASTE,Number=0,Type=Flag,"
        'Description="Genomic origin of interspersed duplication seems to be deleted">'
    ),
    '##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">',
    '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">',
    '##INFO=<ID=SUPPORT,Number=1,Type=Integer,Description="Number of reads supporting this variant">',
    '##INFO=<ID=STD_SPAN,Number=1,Type=Float,Description="Standard deviation in span of merged SV signatures">',
    '##INFO=<ID=STD_POS,Number=1,Type=Float,Description="Standard deviation in position of merged SV signatures">',
    '##INFO=<ID=STD_POS1,Number=1,Type=Float,Description="Standard deviation of breakend 1 position">',
    '##INFO=<ID=STD_POS2,Number=1,Type=Float,Description="Standard deviation of breakend 2 position">',
    '##FILTER=<ID=hom_ref,Description="Genotype is homozygous reference">',
]


@env.task
async def test_svim_alignment() -> None:
    # upstream case: svim/alignment "homo_sapiens - [bam, bai], fasta - vcf"
    bam, bai, fasta = await asyncio.gather(
        fixture(HUMAN + "nanopore/bam/test.sorted.bam"),
        fixture(HUMAN + "nanopore/bam/test.sorted.bam.bai"),
        fixture(HUMAN + "genome/genome.fasta"),
    )
    vcf = await svim_alignment(bam, bai, fasta, prefix="test", args="--min_mapq 5")
    assert vcf.path.endswith("/test.vcf"), vcf.path
    async with vcf.open("rb") as fh:
        lines = bytes(await fh.read()).decode().splitlines()
    assert lines[0] == "##fileformat=VCFv4.2", lines[0]
    assert lines[2:21] == HEADER, lines[2:21]


tests = [test_svim_alignment]
