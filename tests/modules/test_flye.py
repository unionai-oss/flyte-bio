"""Tests for flyte_bio.modules.flye.

The cases and expectations are upstream's flye tests (modules commit 9a1cc9f,
Flye 2.9.5), all on the same PacBio HiFi reads in different read modes, with
the module's 4 CPUs as upstream's tests. Flye's heuristics make its assemblies
vary (upstream points at Flye issue #164), and upstream snapshots no assembly
md5s. Its stated intent is to check that ``assembly_info.txt`` contains
``contig_1``, but it tests whether the file's second line contains the text
``contig_1`` (true in three modes of its snapshot, false for
``--pacbio-hifi``). That line is whichever contig Flye lists first, which
varies between runs: on this cluster it was ``contig_9`` with
``--pacbio-corr``, ``contig_6`` with ``--nano-hq`` and ``contig_12`` with
``--pacbio-hifi``, the opposite of upstream's snapshot in all three. So this
test checks the intent instead: every output is written, one
``assembly_info.txt`` row is ``contig_1``, the assembly's contigs are exactly
the rows of ``assembly_info.txt``, and ``params.json`` has upstream's md5.
"""

import gzip
import os

from flyte.io import File

from flyte_bio.modules.flye import flye
from tests.framework import assert_md5, env, fixture

HIFI = "genomics/homo_sapiens/pacbio/fastq/test_hifi.fastq.gz"
PARAMS_MD5 = "54b576cb6d4d27656878a7fd3657bde9"


async def text(file: File) -> str:
    async with file.open("rb") as fh:
        return bytes(await fh.read()).decode()


async def check(mode: str) -> None:
    result = await flye(await fixture(HIFI), mode=mode, prefix="test")
    names = [f.name or os.path.basename(f.path) for f in await result.results.list_files()]
    for name in [
        "test.assembly.fasta.gz",
        "test.assembly_graph.gfa.gz",
        "test.assembly_graph.gv.gz",
        "test.flye.log",
        "test.assembly_info.txt",
        "test.params.json",
    ]:
        assert name in names, f"flye wrote no {name} (got {names})"
    rows = [line.split("\t")[0] for line in (await text(result.txt)).splitlines()[1:] if line]
    assert "contig_1" in rows, f"assembly_info names no contig_1 (rows: {rows})"
    async with result.fasta.open("rb") as fh:
        headers = [
            line[1:].split()[0]
            for line in gzip.decompress(bytes(await fh.read())).decode().splitlines()
            if line.startswith(">")
        ]
    assert sorted(headers) == sorted(rows), f"assembly contigs {headers} differ from assembly_info rows {rows}"
    await assert_md5(result.json, PARAMS_MD5, label="params.json")


@env.task
async def test_flye_pacbio_corr() -> None:
    await check("--pacbio-corr")


@env.task
async def test_flye_pacbio_hifi() -> None:
    await check("--pacbio-hifi")


@env.task
async def test_flye_nano_corr() -> None:
    await check("--nano-corr")


@env.task
async def test_flye_nano_hq() -> None:
    await check("--nano-hq")


tests = [test_flye_pacbio_corr, test_flye_pacbio_hifi, test_flye_nano_corr, test_flye_nano_hq]
