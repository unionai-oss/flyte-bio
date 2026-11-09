"""Tests for flyte_bio.modules.checkv.

The cases and expected md5s are upstream's checkv/downloaddatabase and
checkv/endtoend tests (modules commit 9a1cc9f, CheckV 1.0.3).

The download case downloads the full database, as upstream; it's cached, so
it runs once. Upstream checks the database's top-level entries
(``README.txt``, ``genome_db``, ``hmm_db``); here they're the first members
of the archive, read from its first bytes only.

The end-to-end case runs the gunzipped Bacteroides fragilis contigs against
upstream's minimal test database (``checkv_minimal_db.tar``, about 47 MB),
passed as the archive it is. ``proviruses.fna`` is empty and checked by name,
as upstream. Upstream's second end-to-end case, named "stub", runs without
``-stub`` and so repeats the first; it isn't repeated here.
"""

import io
import tarfile
import zlib

from flyte.io import File

from flyte_bio.modules.checkv import checkv_downloaddatabase, checkv_endtoend
from flyte_bio.modules.gunzip import gunzip
from tests.framework import assert_md5, env, fixture

MINIMAL_DB = (
    "https://github.com/nf-core/test-datasets/raw/phageannotator/modules/nfcore/checkv/endtoend/checkv_minimal_db.tar"
)
ENDTOEND_MD5S = {
    "quality_summary": "104c579b018158b1c3f016adb0b69eb3",
    "completeness": "7b252af89a74e8ac83a7d81727f2f716",
    "contamination": "cd37c2aa2aaa60cdcd798a0d989cdac4",
    "complete_genomes": "4afbc73b59d01b7c9fac30432a343eaf",
    "viruses": "ea1d49a630943747f21d935d83ab4399",
}


@env.task
async def checkv_archive_first_members(archive: File, head_bytes: int = 1 << 20) -> list[str]:
    """The names of the members at the start of the ``.tar.gz`` ``archive``, read from its first bytes only."""
    async with archive.open("rb") as fh:
        head = bytes(await fh.read(head_bytes))
    data = zlib.decompressobj(16 + zlib.MAX_WBITS).decompress(head)
    names = []
    try:
        with tarfile.open(fileobj=io.BytesIO(data), mode="r|") as tar:
            for member in tar:
                names.append(member.name)
    except (tarfile.TarError, EOFError):
        pass  # the head ends inside a member
    return names


@env.task
async def test_checkv_downloaddatabase() -> None:
    # upstream case: checkv/downloaddatabase "No input: download only"
    db = await checkv_downloaddatabase()
    assert db.name == "checkv_db.tar.gz", f"database archive is {db.name}"
    names = await checkv_archive_first_members(db)
    assert names[:3] == ["README.txt", "genome_db", "hmm_db"], f"archive starts with {names[:3]}"


@env.task
async def test_checkv_endtoend_contigs() -> None:
    # upstream case: checkv/endtoend "bacteroides_fragilis - illumina - test1_contigs_fa_gz"
    contigs = await gunzip(
        await fixture("genomics/prokaryotes/bacteroides_fragilis/illumina/fasta/test1.contigs.fa.gz")
    )
    result = await checkv_endtoend(contigs, await fixture(MINIMAL_DB), prefix="test")
    for field, md5 in ENDTOEND_MD5S.items():
        file = getattr(result, field)
        await assert_md5(file, md5, label=file.name)
    assert result.proviruses.name == "proviruses.fna", result.proviruses.name


tests = [test_checkv_downloaddatabase, test_checkv_endtoend_contigs]
