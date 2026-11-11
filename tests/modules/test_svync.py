"""Tests for flyte_bio.modules.svync.

The case is upstream's svync test (modules commit 9a1cc9f, svync 0.3.0),
with its inline ``config.yaml``: nft-vcf's ``variantsMD5`` of the output (the
md5 of its record lines joined without newlines) and the index's name.
"""

import gzip
import hashlib
import tempfile
from pathlib import Path

from flyte.io import File

from flyte_bio.modules.svync import svync
from tests.framework import env, fixture

SIMULATED_SV = "genomics/homo_sapiens/illumina/vcf/chr21/simulated_sv.vcf.gz"
CONFIG = """info:
  SVMETHOD:
    value: svync
    description: Type of approach used to detect SV
    number: 1
    type: string
format:
  GT:
    value: 1/1
    description: Genotype
    number: 1
    type: string
"""


@env.task
async def test_svync_vcf_config() -> None:
    config = Path(tempfile.mkdtemp()) / "config.yaml"
    config.write_text(CONFIG)
    vcf = await fixture(SIMULATED_SV)
    result = await svync(
        vcf, await File.from_local(str(config)), prefix="test", tbi=await fixture(SIMULATED_SV + ".tbi")
    )
    assert (result.vcf.name, result.tbi.name) == ("test.vcf.gz", "test.vcf.gz.tbi"), (result.vcf, result.tbi)
    async with result.vcf.open("rb") as fh:
        lines = gzip.decompress(bytes(await fh.read())).decode().splitlines()
    records = [x for x in lines if not x.startswith("#")]
    md5 = format(int(hashlib.md5("".join(records).encode()).hexdigest(), 16), "x")
    assert md5 == "2bce0e43fb6a93283f76f1464130bfc3", f"variantsMD5 {md5}"


tests = [test_svync_vcf_config]
