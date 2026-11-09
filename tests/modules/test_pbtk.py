"""Tests for flyte_bio.modules.pbtk.

The case and expected values are upstream's pbtk/pbmerge test (modules
commit 9a1cc9f, pbtk 3.1.1): merging the Iso-Seq clustered reads and their
singletons. As upstream, it checks the merged BAM's header lines (read with
``samtools view -H --no-PG``, so samtools adds none), its reads' sequences
(nft-bam's ``getReadsMD5``) and the ``.pbi``'s name. The BAM's ``@HD`` line
says ``VN:1.5``; upstream's snapshot shows ``VN:1.6`` because htsjdk, which
nft-bam reads the header with, writes its own SAM version there, so that one
field is compared as htsjdk reports it.
"""

import re

from flyte_bio.modules.pbtk import pbtk_pbmerge
from flyte_bio.modules.samtools import samtools_view
from tests.framework import assert_reads_md5, env, fixture

BAMS = "genomics/homo_sapiens/pacbio/bam/alz.ccs.fl.NEB_5p--NEB_Clontech_3p.flnc.clustered"
HEADER = [
    "@HD\tVN:1.6\tSO:unknown\tpb:3.0.1",
    "@RG\tID:e4927d21\tPL:PACBIO\tDS:READTYPE=TRANSCRIPT\tPU:transcript\tPM:SEQUEL",
    (
        "@PG\tID:ccs-6.0.0\tPN:ccs\tVN:6.0.0\tDS:Generate circular consensus sequences (ccs) from subreads.\t"
        "CL:ccs ccs alz.1perc.subreads.100000.bam alz.1perc.subreads.100000.1.ccs.bam "
        "--report-file alz.1perc.subreads.100000.1.ccs_report.txt "
        "--report-json alz.1perc.subreads.100000.1.ccs_report.json "
        "--metrics-json alz.1perc.subreads.100000.1.zmw_metrics.json.gz --chunk 1/1 -j 2 --min-rq 0.9"
    ),
    (
        "@PG\tID:cluster\tVN:3.4.0 (commit v3.4.0)\tCL:cluster "
        "alz.1perc.subreads.100000.1.ccs.fl.NEB_5p--NEB_Clontech_3p.flnc.bam "
        "alz.1perc.subreads.100000.1.ccs.fl.NEB_5p--NEB_Clontech_3p.flnc.clustered.bam "
        "--singletons --singletons --use-qvs --verbose"
    ),
    (
        "@PG\tID:lima\tVN:2.2.0 (commit v2.2.0)\tCL:lima alz.1perc.subreads.100000.1.ccs.bam primers.fasta "
        "alz.1perc.subreads.100000.1.ccs.fl.bam -j 2 --isoseq --peek-guess"
    ),
    (
        "@PG\tID:refine\tVN:3.4.0 (commit v3.4.0)\tCL:refine -j 2 "
        "alz.1perc.subreads.100000.1.ccs.fl.NEB_5p--NEB_Clontech_3p.bam primers.fasta "
        "alz.1perc.subreads.100000.1.ccs.fl.NEB_5p--NEB_Clontech_3p.flnc.bam"
    ),
    "@PG\tID:pbmerge-3.1.1\tPN:pbmerge\tVN:3.1.1",
]


@env.task
async def test_pbtk_pbmerge_bam() -> None:
    bams = [await fixture(BAMS + ".bam"), await fixture(BAMS + ".singletons.bam")]
    result = await pbtk_pbmerge(bams, prefix="test")
    assert result.bam.name == "test.bam", result.bam.name
    assert result.pbi is not None and result.pbi.name == "test.bam.pbi", result.pbi
    header = await samtools_view(bam=result.bam, args="-H --no-PG")
    async with header.open("rb") as fh:
        lines = bytes(await fh.read()).decode().splitlines()
    # htsjdk reports its own SAM version in @HD.
    lines = [re.sub(r"^(@HD\tVN:)[0-9.]+", r"\g<1>1.6", line) for line in lines]
    assert lines == HEADER, lines
    await assert_reads_md5(result.bam, "1ec5c2c66d16d0236c20d8e26d0109ea", label="test.bam")


tests = [test_pbtk_pbmerge_bam]
