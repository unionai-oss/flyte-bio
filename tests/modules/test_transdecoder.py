"""Tests for flyte_bio.modules.transdecoder.

The cases and expected values are upstream's transdecoder/longorf and
transdecoder/predict tests (modules commit 9a1cc9f, TransDecoder 5.7.1) on
the SARS-CoV-2 genome, plain and gzipped. LongOrfs' four files are checked by
md5 (the same for both inputs). Predict's ORFs come out in a random order,
so, as upstream, its cases check which outputs mention ORF
``GENE.MT192765.1~~MT192765.1.p7``: only the GFF3 does.

Upstream's third Predict case reruns the plain one with copy staging, a
Nextflow staging mode; inputs here are always copies, so it's the plain case.
"""

from flyte.io import File

from flyte_bio.modules.transdecoder import transdecoder_longorf, transdecoder_predict
from tests.framework import assert_md5, env, fixture

GENOME = "genomics/sarscov2/genome/genome.fasta"
ORF_ID = "ID=GENE.MT192765.1~~MT192765.1.p7"


async def text(file: File) -> str:
    async with file.open("rb") as fh:
        return bytes(await fh.read()).decode()


async def check_longorf(path: str) -> None:
    result = await transdecoder_longorf(await fixture(path))
    assert result.folder.path.rstrip("/").endswith("/genome.fasta.transdecoder_dir"), result.folder.path
    await assert_md5(result.cds, "9b24da198e1916b9ac3c33b0dc74642c", label="longest_orfs.cds")
    await assert_md5(result.dat, "5226ef383532e9b16c16b03e35ce181e", label="base_freqs.dat")
    await assert_md5(result.gff3, "82cf8837bf30feac1a2a2309e10e90db", label="longest_orfs.gff3")
    await assert_md5(result.pep, "fd301d53b23d739ad4d608d0b00b7996", label="longest_orfs.pep")


@env.task
async def test_transdecoder_longorf() -> None:
    await check_longorf(GENOME)


@env.task
async def test_transdecoder_longorf_gz_input() -> None:
    await check_longorf(GENOME + ".gz")


async def check_predict(path: str) -> None:
    fasta = await fixture(path)
    longorfs = await transdecoder_longorf(fasta)
    result = await transdecoder_predict(fasta, longorfs.folder)
    outputs = {"bed": result.bed, "cds": result.cds, "gff3": result.gff3, "pep": result.pep}
    for kind, file in outputs.items():
        assert file.name == f"genome.fasta.transdecoder.{kind}", f"{kind} output is {file.name}"
    mentions = {kind: ORF_ID in await text(file) for kind, file in outputs.items()}
    assert mentions == {"bed": False, "cds": False, "gff3": True, "pep": False}, f"{ORF_ID} in: {mentions}"


@env.task
async def test_transdecoder_predict_uncompressed() -> None:
    await check_predict(GENOME)


@env.task
async def test_transdecoder_predict_compressed() -> None:
    await check_predict(GENOME + ".gz")


tests = [
    test_transdecoder_longorf,
    test_transdecoder_longorf_gz_input,
    test_transdecoder_predict_uncompressed,
    test_transdecoder_predict_compressed,
]
