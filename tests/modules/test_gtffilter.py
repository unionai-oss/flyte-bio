"""Tests for flyte_bio.modules.gtffilter.

The filtered GTF is a deterministic subset of the input lines, so it's a
clean md5 anchor (and the first exercise of the vendored-script mechanism).
md5 from the upstream custom/gtffilter snapshot (pinned version).
"""



from flyte_bio.modules.gtffilter import gtf_filter
from tests.framework import assert_md5, env, fixture


@env.task
async def test_gtf_filter() -> None:
    # upstream case: custom/gtffilter "test_custom_gtffilter"
    gtf = await fixture("genomics/sarscov2/genome/genome.gtf")
    fasta = await fixture("genomics/sarscov2/genome/genome.fasta")
    out = await gtf_filter(gtf=gtf, fasta=fasta)
    await assert_md5(out, "aa8b2aa1e0b5fbbba3b04d471e1b0535", label="gtf filter")


tests = [test_gtf_filter]
