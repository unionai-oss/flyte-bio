"""Tests for flyte_bio.modules.pcgr.

Upstream's pcgr/getref (modules commit 9a1cc9f) has no real test: its one
case, bundle ``20250314`` for ``grch38``, is commented out (the bundle is
about 5 GB), leaving only a stub run. This test runs that case for real (the
download is cached, so it runs once) and checks the archive is named after
the bundle version, starts with the bundle's version file, which reads
``20250314``, and holds upstream's layout, ``20250314/data/grch38/...``. It
reads only the archive's first bytes.
"""

import io
import tarfile
import zlib

from flyte.io import File

from flyte_bio.modules.pcgr import pcgr_getref
from tests.framework import env

BUNDLE_VERSION = "20250314"
GENOME = "grch38"


@env.task
async def pcgr_archive_head(archive: File, head_bytes: int = 1 << 20) -> tuple[list[str], str]:
    """The members at the start of the ``.tar.gz`` ``archive``, and its first member's text."""
    async with archive.open("rb") as fh:
        head = bytes(await fh.read(head_bytes))
    data = zlib.decompressobj(16 + zlib.MAX_WBITS).decompress(head)
    names = []
    first_text = ""
    try:
        with tarfile.open(fileobj=io.BytesIO(data), mode="r|") as tar:
            for member in tar:
                names.append(member.name)
                if len(names) == 1 and member.isfile():
                    first_text = tar.extractfile(member).read().decode()
    except (tarfile.TarError, EOFError):
        pass  # the head ends inside a member
    return names, first_text


@env.task
async def test_pcgr_getref_grch38() -> None:
    archive = await pcgr_getref(BUNDLE_VERSION, GENOME)
    assert archive.name == f"{BUNDLE_VERSION}.tar.gz", f"bundle archive is {archive.name}"
    names, version = await pcgr_archive_head(archive)
    first = f"{BUNDLE_VERSION}/data/{GENOME}/.PCGR_BUNDLE_VERSION"
    assert names[:2] == [first, BUNDLE_VERSION], f"archive starts with {names[:2]}"
    assert version.strip() == BUNDLE_VERSION, f"{first} reads {version!r}"
    prefix = f"{BUNDLE_VERSION}/data/{GENOME}/"
    folders = {f"{BUNDLE_VERSION}/data", f"{BUNDLE_VERSION}/data/{GENOME}"}
    outside = [n for n in names[2:] if n not in folders and not n.startswith(prefix)]
    assert len(names) > 2 and not outside, f"members outside {prefix}: {outside[:5]}"


tests = [
    test_pcgr_getref_grch38,
]
