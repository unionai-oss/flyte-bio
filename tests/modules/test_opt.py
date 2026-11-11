"""Tests for flyte_bio.modules.opt.

The cases and expected md5s are upstream's opt/flip, opt/track and opt/stat
tests (modules commit 9a1cc9f, opt 0.0.1): a probe panel checked against the
human minigenome and its GTF. As upstream, the track test runs flip first, and
the stat test runs flip and track first, with the gene synonym table.
"""

from flyte.io import File

from flyte_bio.modules.opt import opt_flip, opt_stat, opt_track
from tests.framework import assert_md5, env, fixture

PROBES = "spatial_omics/probe_detection/probes.fa"
SYNONYMS = "spatial_omics/probe_detection/gene_synonyms.csv"
GTF = "genomics/homo_sapiens/genome/minigenome.gtf"
FASTA = "genomics/homo_sapiens/genome/minigenome.fa"


async def opt_flipped(prefix: str) -> tuple[File, File, File]:
    gtf, fasta = await fixture(GTF), await fixture(FASTA)
    return await opt_flip(await fixture(PROBES), gtf, fasta, prefix=prefix), gtf, fasta


@env.task
async def test_opt_flip() -> None:
    # upstream case: "testrun opt/flip panel probe sequences"
    flipped, _, _ = await opt_flipped("testrun_opt_flip")
    await assert_md5(flipped, "6e880a8dbc3a6cc998caaeacd35003fa", label="fwd_oriented.fa")


@env.task
async def test_opt_track() -> None:
    # upstream case: "testrun opt/track panel probe sequences"
    flipped, gtf, fasta = await opt_flipped("testrun_opt_track")
    targets = await opt_track(flipped, gtf, fasta, prefix="testrun_opt_track")
    await assert_md5(targets, "e15465df3845d7a6acf64dd3be04391b", label="probe2targets.tsv")


@env.task
async def test_opt_stat() -> None:
    # upstream case: "testrun opt/stat panel probe sequences"
    flipped, gtf, fasta = await opt_flipped("testrun_opt_stat")
    targets = await opt_track(flipped, gtf, fasta, prefix="testrun_opt_stat")
    summary = await opt_stat(targets, flipped, await fixture(SYNONYMS), prefix="testrun_opt_stat")
    await assert_md5(summary, "b2884d9c8c89124d3cbbbc1223d81c99", label="collapsed_summary.tsv")


tests = [test_opt_flip, test_opt_track, test_opt_stat]
