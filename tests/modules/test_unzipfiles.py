"""Tests for flyte_bio.modules.unzipfiles.

The case and expected md5 are upstream's unzipfiles test (modules commit
9a1cc9f, p7zip 16.02).
"""

import os

from flyte_bio.modules.unzipfiles import unzipfiles
from tests.framework import assert_md5, env, fixture


@env.task
async def test_unzipfiles_tar_gz() -> None:
    # upstream case "generic [tar] [tar_gz]": the gzip layer comes off, leaving hello.tar
    result = await unzipfiles(await fixture("generic/tar/hello.tar.gz"))
    names = [os.path.basename(f.path) for f in result.files]
    assert names == ["hello.tar"], names
    await assert_md5(result.files[0], "80c66db79a773bc87b3346035ff9593e", label="hello.tar")


tests = [test_unzipfiles_tar_gz]
