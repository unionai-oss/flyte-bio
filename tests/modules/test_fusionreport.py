"""Tests for flyte_bio.modules.fusionreport.

The cases are upstream's fusionreport/download and fusionreport/detect tests
(modules commit 9a1cc9f, fusion-report 4.1.2). Both run without COSMIC
(``--no-cosmic``), as upstream's do: COSMIC needs an account.

- download: upstream snapshots nothing of the databases (``unstableKeys``),
  as they're downloaded live. FusionGDB2's server
  (``compbio.uth.edu/FusionGDB/combined_tables``) answers every scripted
  request with HTTP 429 (Too Many Requests), from the cluster and elsewhere
  (checked 2026-11-09), which fails ``fusion_report download``. So this test
  also leaves out FusionGDB2 (``--no-fusiongdb2``) and checks that the
  Mitelman database came.
- detect: the merged Arriba, STAR-Fusion and FusionCatcher test calls with a
  tool cutoff of 1 and ``--export csv``. Upstream snapshots the md5s of both
  fusion lists, the CSV and MultiQC's JSON, and the names of the HTML pages
  (here the members of the pages' archive). ``fusion_report run`` needs the
  FusionGDB2 table (it fails with "no such table: fusiongdb2" without it,
  even with ``--no-fusiongdb2``), so while FusionGDB2 can't be downloaded this
  test can't run: it is kept, as upstream's case, in ``blocked_tests``, not
  ``tests``. The databases are downloaded live, so its md5s hold while the
  upstream databases don't change.
"""

import tarfile

from flyte_bio.modules.fusionreport import fusionreport_detect, fusionreport_download
from tests.framework import assert_md5, env, fixture

HUMAN = "genomics/homo_sapiens/genome/"
HTML_PAGES = [
    "AKAP9_BRAF.html",
    "BRD4-1_NUTM1.html",
    "BRD4_NUTM1.html",
    "CD74_AL132671.2.html",
    "CD74_ROS1.html",
    "CIC_DUX4.html",
    "DUX4_IGH@.html",
    "EML4_ALK.html",
    "ETV6_NTRK3.html",
    "EWSR1_ATF1.html",
    "EWSR1_FLI1.html",
    "FGFR3_TACC3.html",
    "FIP1L1_PDGFRA.html",
    "GOPC_ROS1.html",
    "HOOK3_RET.html",
    "IGH@_CRLF2.html",
    "MALT1_IGH@.html",
    "NPM1_ALK.html",
    "TMPRSS2_ETV1.html",
    "test_sample_fusionreport_index.html",
]


@env.task
async def test_fusionreport_download() -> None:
    # upstream case: fusionreport/download "Download fusionreport databases", without FusionGDB2 (see above)
    db = await fusionreport_download(args="--no-cosmic --no-fusiongdb2")
    names = sorted(f.name for f in await db.list_files())
    assert "mitelman.db" in names, names
    assert "cosmic.db" not in names and "fusiongdb2.db" not in names, names


@env.task
async def test_fusionreport_detect() -> None:
    # upstream case: fusionreport/detect "FUSIONREPORT - arriba - starfusion - fusioncatcher"
    db = await fusionreport_download(args="--no-cosmic")
    result = await fusionreport_detect(
        db,
        prefix="test_sample",
        tools_cutoff=1,
        arriba_fusions=await fixture(HUMAN + "arriba.tsv"),
        starfusion_fusions=await fixture(HUMAN + "starfusion.tsv"),
        fusioncatcher_fusions=await fixture(HUMAN + "fusioncatcher.txt"),
        args="--no-cosmic --export csv",
    )
    await assert_md5(result.fusion_list, "3593b7021f26cc5427fdc96f0d1c72f0", label="test_sample.fusionreport.tsv")
    await assert_md5(
        result.fusion_list_filtered, "3593b7021f26cc5427fdc96f0d1c72f0", label="test_sample.fusionreport_filtered.tsv"
    )
    assert result.csv is not None and result.csv.name == "test_sample.fusions.csv", result.csv
    await assert_md5(result.csv, "74b4672cfa6f619df964cb4c7ec07735", label="test_sample.fusions.csv")
    assert [f.name for f in result.json] == ["fusion_genes_mqc.json"], [f.name for f in result.json]
    await assert_md5(result.json[0], "cc715e5007c345b831e9f60ec6c55018", label="fusion_genes_mqc.json")
    assert result.report.name == "test_sample_fusionreport_index.html", result.report.name
    with tarfile.open(await result.html.download()) as tar:
        pages = sorted(name.removeprefix("./") for name in tar.getnames())
    assert pages == HTML_PAGES, pages


tests = [test_fusionreport_download]
# Needs FusionGDB2, whose server refuses scripted downloads for now (see above): run it when it serves again.
blocked_tests = [test_fusionreport_detect]
