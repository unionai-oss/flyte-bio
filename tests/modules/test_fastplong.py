"""Tests for flyte_bio.modules.fastplong.

The case and expected md5s are upstream's fastplong test (modules commit
9a1cc9f, fastplong 0.4.1). nf-test computes the reads' md5 over the
decompressed FASTQ.

The JSON report records the command line, ``--thread`` included; upstream's
matches 2 threads, so the test runs on 2 CPUs.
"""

import flyte
from flyte.io import File

from flyte_bio.modules.fastplong import fastplong_cmd, read_fastplong_results
from tests.framework import assert_gunzipped_md5, assert_md5, env, fixture


async def text(file: File) -> str:
    async with file.open("rb") as fh:
        return bytes(await fh.read()).decode()


@env.task
async def test_fastplong_pacbio() -> None:
    reads = await fixture("genomics/homo_sapiens/pacbio/fastq/alz.ccs.fastq.gz")
    # The cache key ignores resources, so salt it to keep the 2-CPU run apart from a default one.
    task = fastplong_cmd.override(
        resources=flyte.Resources(cpu=2, memory="8Gi"), cache=flyte.Cache(behavior="auto", salt="cpu=2")
    )
    results = await task(
        reads=reads,
        adapter_fasta=None,
        prefix="alz.ccs",
        discard_trimmed_pass=False,
        save_trimmed_fail=False,
        args="",
    )
    result = await read_fastplong_results(results, "alz.ccs")
    total = "<tr><td class='col1'>total bases:</td><td class='col2'>272.128000 K</td></tr>"
    assert total in await text(result.html), "fastplong html lacks the total bases"
    assert "reads passed filter: 100" in await text(result.log), "fastplong log lacks 'reads passed filter: 100'"
    await assert_md5(result.json, "5d8fb4ca4a7e1282d843040baae98f43", label="alz.ccs.fastplong.json")
    assert result.reads is not None and result.reads_fail is None
    await assert_gunzipped_md5(result.reads, "23a6d8f301d0fd4d2da21f86ff0afac6", label="alz.ccs.fastplong.fastq.gz")


tests = [test_fastplong_pacbio]
