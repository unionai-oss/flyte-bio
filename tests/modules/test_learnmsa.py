"""Tests for flyte_bio.modules.learnmsa.

The case is upstream's learnmsa/align test (modules commit 9a1cc9f, learnMSA
2.0.16) on the seatoxin sequences. Upstream treats the alignment as unstable
and snapshots only its name; this checks the name and that every sequence is
in the alignment.
"""

from flyte_bio.modules.learnmsa import learnmsa_align
from tests.framework import env, fixture

SEQUENCES = "https://raw.githubusercontent.com/nf-core/test-datasets/multiplesequencealign/testdata/setoxin-ref.fa"
SEATOXIN_IDS = ("1apf", "1ahl", "1atx", "1sh1", "1bds")


@env.task
async def test_learnmsa_seatoxin_fa() -> None:
    out = await learnmsa_align(await fixture(SEQUENCES), prefix="test")
    assert out.name == "test.aln", f"output named {out.name}"
    async with out.open("rb") as fh:
        text = bytes(await fh.read()).decode()
    missing = [i for i in SEATOXIN_IDS if f">{i}" not in text]
    assert not missing, f"test.aln lacks {missing}"


tests = [test_learnmsa_seatoxin_fa]
