"""Tests for flyte_bio.modules.unzip.

The expected md5 is the upstream unzip snapshot value.
"""

import os

from flyte_bio.modules.unzip import unzip
from tests.framework import assert_md5, env, fixture


@env.task
async def test_unzip_tar_gz() -> None:
    # upstream case: unzip "generic [tar] [tar_gz]" — the gzip layer comes off, leaving hello.tar
    out = await unzip(await fixture("generic/tar/hello.tar.gz"))
    names = sorted(os.path.basename(f.path) for f in await out.list_files())
    assert names == ["hello.tar"], names
    tar = await out.get_file("hello.tar")
    await assert_md5(tar, "80c66db79a773bc87b3346035ff9593e", label="hello.tar")


tests = [test_unzip_tar_gz]
