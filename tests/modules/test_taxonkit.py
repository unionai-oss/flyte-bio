"""Tests for flyte_bio.modules.taxonkit.

The cases and expected md5 are upstream's taxonkit/lineage tests (modules
commit 9a1cc9f, TaxonKit 0.18.0): SARS-CoV-2's lineage (taxon 2697049), given
as an ID and as a file of IDs, against NCBI's current taxdump. Upstream
fetches the taxdump over FTP; here it's the same file over HTTPS, passed as
the archive. As upstream's, these tests follow NCBI's live taxonomy, so the md5
would change if NCBI renamed a taxon in the lineage.
"""

import tempfile
from pathlib import Path

from flyte.io import File

from flyte_bio.modules.taxonkit import taxonkit_lineage
from tests.framework import assert_md5, env, fixture

TAXDUMP = "https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz"
LINEAGE_MD5 = "da0998a741219814ae2cfad153f0ecc2"


@env.task
async def test_taxonkit_lineage_taxid() -> None:
    out = await taxonkit_lineage(await fixture(TAXDUMP), taxid="2697049", prefix="test")
    assert out.name == "test.tsv", f"output is {out.name}"
    await assert_md5(out, LINEAGE_MD5, label="test.tsv")


@env.task
async def test_taxonkit_lineage_taxid_list() -> None:
    # As upstream's collectFile: the ID and a newline, in taxids.txt.
    ids = Path(tempfile.mkdtemp()) / "taxids.txt"
    ids.write_text("2697049\n")
    out = await taxonkit_lineage(await fixture(TAXDUMP), taxidfile=await File.from_local(str(ids)), prefix="test")
    assert out.name == "test.tsv", f"output is {out.name}"
    await assert_md5(out, LINEAGE_MD5, label="test.tsv")


tests = [
    test_taxonkit_lineage_taxid,
    test_taxonkit_lineage_taxid_list,
]
