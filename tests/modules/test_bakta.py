"""Tests for flyte_bio.modules.bakta.

The cases are upstream's bakta/baktadbdownload and bakta/bakta tests (modules
commit 9a1cc9f, bakta 1.12.0). As upstream, both download the light database
(``--type light``); the download is cached, so it runs once. Here the
database is an archive, whose first members (``version.json``, then
``bakta.db``) the download test reads from its start. Upstream checks the
database exists and bakta's outputs contain known sequences and
headers rather than md5s, as the outputs record bakta's and the database's
versions; these tests check the same.
"""

import io
import tarfile
import zlib

from flyte.io import File

from flyte_bio.modules.bakta import bakta, bakta_db_download
from tests.framework import env, fixture

LIGHT_DB_ARGS = "--type light"


async def text(file: File) -> str:
    async with file.open("rb") as fh:
        return bytes(await fh.read()).decode()


async def light_database() -> File:
    return await bakta_db_download(args=LIGHT_DB_ARGS)


@env.task
async def bakta_archive_first_members(archive: File, head_bytes: int = 1 << 20) -> list[str]:
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
async def test_bakta_db_download() -> None:
    db = await light_database()
    assert db.name == "bakta_db.tar.gz", f"database archive is {db.name}"
    names = await bakta_archive_first_members(db)
    assert names[:2] == ["version.json", "bakta.db"], f"archive starts with {names[:2]}"


@env.task
async def test_bakta_bacteroides_fragilis_genome() -> None:
    contigs = await fixture("genomics/prokaryotes/bacteroides_fragilis/illumina/fasta/test1.contigs.fa.gz")
    result = await bakta(contigs, await light_database(), prefix="test")
    expected = [
        (result.embl, '/translation="MKNTLKIAILLIAIISMGHWMPVKQVCDLNSLSLQNVEALANGET'),
        (result.faa, "MKNTLKIAILLIAIISMGHWMPVKQVCDLNSLSLQNVEALANGETPNYTFCIGAGSVDCPIQHDKVKYVSQGFSLDY"),
        (result.ffn, "ATGAAAAACACTTTAAAAATAGCTATTCTTCTTATTGCTATTATTTCTATGGGGCATTGGATGCCTGTAAAACAAGT"),
        (result.fna, "TCTTTTTACTCATAATCTACTTTTATGATGTTAATTATTTTTTCCGTGTCTCTCTTTCGG"),
        (result.gbff, '/translation="MKNTLKIAILLIAIISMGHWMPVKQVCDLNSLSLQNVEALANGET'),
        (result.gff, "##sequence-region contig_1 1 2926"),
        (result.hypotheticals_tsv, "#Annotated with Bakta"),
        (result.hypotheticals_faa, "MKNLILVLGCFFFLISCQQTEKEKLEELVKNWNGKEVLL"),
        (result.tsv, "SO:0001217, UniRef:UniRef50_A0A0I9S7A3"),
        (result.txt, "Length: 1739120"),
        (result.json, "features"),
    ]
    for file, needle in expected:
        assert needle in await text(file), f"{file.name} lacks {needle!r}"


tests = [
    test_bakta_db_download,
    test_bakta_bacteroides_fragilis_genome,
]
