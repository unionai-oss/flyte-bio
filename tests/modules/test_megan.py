"""Tests for flyte_bio.modules.megan.

The case is upstream's megan/rma2info test (modules commit 9a1cc9f, MEGAN
Community Edition 6.25.9): upstream's test RMA6 file
(``delete_me/malt/test.rma6``) reported with ``-c2c Taxonomy`` and a MEGAN
summary. As upstream (both outputs are unstable), it checks the outputs by
name.
"""

from flyte_bio.modules.megan import megan_rma2info
from tests.framework import env, fixture


@env.task
async def test_megan_rma2info() -> None:
    # upstream case: megan/rma2info "malt - rma"
    rma6 = await fixture("delete_me/malt/test.rma6")
    result = await megan_rma2info(rma6, megan_summary=True, prefix="test", args="-c2c Taxonomy")
    assert result.txt.name == "test.txt.gz", result.txt.name
    assert result.megan_summary is not None, "rma2info wrote no test.megan"
    assert result.megan_summary.name == "test.megan", result.megan_summary.name


tests = [test_megan_rma2info]
