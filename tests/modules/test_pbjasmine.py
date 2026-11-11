"""Tests for flyte_bio.modules.pbjasmine.

The case and expected md5 are upstream's pbjasmine test (modules commit
9a1cc9f, jasmine 26.1.3). nf-test hashes the output BAM as it is, and its
header's ``@PG`` command line records ``--num-threads``: upstream's snapshot
is of a 2-thread run (reproduced with 2 threads, not 1 or 4), so the test
runs on 2 CPUs.
"""

import flyte

from flyte_bio.modules.pbjasmine import pbjasmine_cmd
from tests.framework import assert_md5, env, fixture

PACBIO_BAM = "genomics/homo_sapiens/pacbio/bam/NA03697B2_downsampled.pbmm2.repeats.bam"


@env.task
async def test_pbjasmine() -> None:
    # upstream case: pbjasmine "homo_sapiens - bam"
    # The cache key ignores resources, so salt it to keep the 2-CPU run apart from a default one.
    task = pbjasmine_cmd.override(
        resources=flyte.Resources(cpu=2, memory="8Gi"), cache=flyte.Cache(behavior="auto", salt="cpu=2")
    )
    results = await task(bam=await fixture(PACBIO_BAM), prefix="test_jasmine", args="")
    bam = await results.get_file("test_jasmine.bam")
    assert bam is not None, "jasmine wrote no test_jasmine.bam"
    await assert_md5(bam, "45b4823896043a3f606cda869054f915", label="test_jasmine.bam")


tests = [test_pbjasmine]
