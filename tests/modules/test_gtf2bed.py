"""Tests for flyte_bio.modules.gtf2bed.

The BED output is deterministic (sorted by chr/pos/transcript), so it's an
md5 anchor. md5 from the upstream ea-utils/gtf2bed snapshot (pinned
version); the vendored perl was md5-verified locally before embedding.
"""



from flyte_bio.modules.gtf2bed import gtf2bed
from tests.framework import assert_md5, env, fixture


@env.task
async def test_gtf2bed() -> None:
    # upstream case: ea-utils/gtf2bed "homo_sapiens - gtf" (default, no -x)
    gtf = await fixture("genomics/homo_sapiens/genome/genome.gtf")
    out = await gtf2bed(gtf=gtf)
    await assert_md5(out, "ef93285cfee828e5f7b3301b2f49608f", label="gtf2bed")


tests = [test_gtf2bed]
