"""Tests for flyte_bio.modules.last.

The cases and expected values are upstream's last/lastdb, last/train,
last/lastal, last/split, last/mafconvert and last/dotplot tests (modules
commit 9a1cc9f, LAST 1652), with their test configs' ``ext.args``: ``-Q0``
for lastdb, ``-P1`` for lastal, and ``-B`` / ``--ploidy 1 -m`` for
maf-convert's ``bcftools mpileup`` / ``call``. Gzipped outputs are compared
decompressed, as nf-test hashes them.

Upstream's snapshots were made with 2 CPUs: last-train records its ``-P``,
and lastdb's suffix array for the FASTQ's many sequences depends on it. So
those tests pin 2 CPUs, salting the cache (which ignores resources). Upstream
extracts its fixture index into a directory named ``genome`` (the untar
step's ``meta.id``), whose path lastal records, so the lastal tests pass
``index_dir_name="genome"``.

Upstream tests each tool on fixtures (an index and a MAF made by an older
LAST); :func:`test_last_chain` also runs the chain lastdb → last-train →
lastal → last-split → maf-convert / last-dotplot from the genome FASTA.
"""

import asyncio
import gzip
import hashlib
import os
import tempfile
from pathlib import Path

import flyte
from flyte.io import Dir, File

from flyte_bio.modules.last import (
    LastTrainResult,
    last_dotplot,
    last_lastal,
    last_lastdb,
    last_lastdb_cmd,
    last_mafconvert,
    last_split,
    last_train,
    last_train_cmd,
)
from flyte_bio.modules.samtools import samtools_view
from tests.framework import assert_gunzipped_md5, assert_md5, assert_reads_md5, env, fixture, fixture_dir

GENOME = "genomics/sarscov2/genome/"
LAST = GENOME + "alignment/last/"
CONTIGS = "genomics/sarscov2/illumina/fasta/contigs.fasta"
MAF = LAST + "contigs.genome.maf.gz"
INDEX_FILES = ("bck", "des", "prj", "sds", "ssp", "suf", "tis")
MAFCONVERT_ARGS = {"args3": "-B", "args4": "--ploidy 1 -m"}
# Upstream's snapshot CPUs; the salt keeps these runs' cache entries apart from default-resource ones.
TWO_CPUS = {"resources": flyte.Resources(cpu=2, memory="8Gi"), "cache": flyte.Cache(behavior="auto", salt="cpu=2")}


async def last_lastdb_two_cpus(fastx: File, prefix: str, args: str) -> Dir:
    return await last_lastdb_cmd.override(**TWO_CPUS)(fastx=fastx, prefix=prefix, args=args)


async def last_train_two_cpus(fastx: File, index: Dir, prefix: str) -> LastTrainResult:
    task = last_train_cmd.override(**TWO_CPUS)
    results = await task(fastx=fastx, index=index, index_dir_name="lastdb", prefix=prefix, args="")
    param_file, multiqc = await asyncio.gather(
        results.get_file(f"{prefix}.train"), results.get_file(f"{prefix}.train.tsv")
    )
    assert param_file is not None and multiqc is not None, "last-train wrote no .train or .train.tsv"
    return LastTrainResult(param_file=param_file, multiqc=multiqc, results=results)


async def last_fixture_index() -> Dir:
    """Upstream's ``lastdb.tar.gz`` index of the SARS-CoV-2 genome, extracted (``lastdb/genome.*``)."""
    return await fixture_dir(LAST + "lastdb.tar.gz", "lastdb")


async def last_file_md5(file: File, gunzip: bool = False) -> str:
    async with file.open("rb") as fh:
        data = bytes(await fh.read())
    return hashlib.md5(gzip.decompress(data) if gunzip else data).hexdigest()


async def check_index(index: Dir, expected: dict[str, str], label: str) -> None:
    names = sorted(os.path.basename(f.path) for f in await index.list_files())
    assert names == sorted(expected), f"{label}: index files {names}"
    for name, md5 in expected.items():
        found = await index.get_file(name)
        assert found is not None, f"{label}: no {name}"
        await assert_md5(found, md5, label=f"{label} {name}")


async def last_text_file(name: str, text: str) -> File:
    path = Path(tempfile.mkdtemp()) / name
    path.write_text(text)
    return await File.from_local(str(path))


async def check_image(plot: File, name: str, signature: bytes) -> None:
    assert plot.path.endswith("/" + name), plot.path
    async with plot.open("rb") as fh:
        head = bytes(await fh.read(len(signature)))
    assert head == signature, f"{name}: starts with {head!r}"


@env.task
async def test_last_lastdb_fasta() -> None:
    # upstream case: last/lastdb "sarscov2 - fasta"
    index = await last_lastdb_two_cpus(await fixture(GENOME + "genome.fasta"), "test", "-Q0")
    md5s = (
        "157526d333b88523cb15ac4efe00738f",
        "3a9ea6d336e113a74d7fdca5e7b623fc",
        "4171737a2df621788d5fdc8466027205",
        "e7729db27ac7a5a109c9d48cfcdc9015",
        "53524efdea3d8989201419a29e81ec1f",
        "ef7482260705bb8146acbbbdce6c0068",
        "b7c40f06b1309dc6f37849eeb86dfd22",
    )
    await check_index(index, {f"test.{ext}": md5 for ext, md5 in zip(INDEX_FILES, md5s)}, "fasta")


@env.task
async def test_last_lastdb_fastq_gz() -> None:
    # upstream case: last/lastdb "sarscov2 - fastq gzipped"
    index = await last_lastdb_two_cpus(await fixture("genomics/sarscov2/illumina/fastq/test_1.fastq.gz"), "test", "-Q0")
    md5s = (
        "05b14d8ac418b3193d9cc921086cea05",
        "26ab49015cc572172b9efa50fc5190bc",
        "67010829ad4de5848bcd229b113e9275",
        "d3deb4c985081c9f5ad6684d405bd20b",
        "5c17139a9022b0cb97f007146fa1c6da",
        "9ac359afa86a8964d81a87a1d4f05ef0",
        "d57a3a5f7e3e036807356c15bd3aad97",
    )
    await check_index(index, {f"test.{ext}": md5 for ext, md5 in zip(INDEX_FILES, md5s)}, "fastq.gz")


@env.task
async def test_last_train() -> None:
    # upstream case: last/train "sarscov2 - genome - contigs"
    contigs, index = await asyncio.gather(fixture(CONTIGS), last_fixture_index())
    result = await last_train_two_cpus(contigs, index, "contigs")
    assert result.param_file.path.endswith("/contigs.train"), result.param_file.path
    await assert_md5(result.param_file, "dfb9d91777869f805cee030c51c9f94f", label="contigs.train")
    await assert_md5(result.multiqc, "c8721d9cabf1e2b7c0d2d3dac9ed3db8", label="contigs.train.tsv")


@env.task
async def test_last_lastal() -> None:
    # upstream case: last/lastal "sarscov2 - contigs - genome"
    contigs, index = await asyncio.gather(fixture(CONTIGS), last_fixture_index())
    result = await last_lastal(contigs, index, prefix="contigs", args="-P1", index_dir_name="genome")
    assert result.maf.path.endswith("/contigs.maf.gz"), result.maf.path
    await assert_gunzipped_md5(result.maf, "64863a9c1184079fb19cc450dccf3bd3", label="contigs.maf.gz")
    await assert_md5(result.multiqc, "0d6880667f9a7fcb3ff54372b6e558c7", label="contigs.tsv")


@env.task
async def test_last_lastal_trained() -> None:
    # upstream case: last/lastal "sarscov2 - contigs - genome - withparams"
    contigs, index, params = await asyncio.gather(
        fixture(CONTIGS), last_fixture_index(), fixture(LAST + "contigs.genome.par")
    )
    result = await last_lastal(contigs, index, params, prefix="contigs", args="-P1", index_dir_name="genome")
    await assert_gunzipped_md5(result.maf, "d780477da11631ddd8d27128a9588887", label="contigs.maf.gz")
    await assert_md5(result.multiqc, "b0dc88a92948e2f7e73a4d87cfd096d8", label="contigs.tsv")


@env.task
async def test_last_split() -> None:
    # upstream case: last/split "sarscov2 - contigs_genome"
    result = await last_split(await fixture(MAF), prefix="sarscov.contigs.genome")
    assert result.maf.path.endswith("/sarscov.contigs.genome.maf.gz"), result.maf.path
    await assert_gunzipped_md5(result.maf, "689cb18ff7098ff90eaf87017f590208", label="split maf.gz")
    await assert_md5(result.multiqc, "e73bff59577eae51459061686482955c", label="split tsv")


@env.task
async def test_last_mafconvert_psl_sam() -> None:
    # upstream cases: last/mafconvert "sarscov2 - psl" and "sarscov2 - sam"
    maf = await fixture(MAF)
    psl, sam = await asyncio.gather(
        last_mafconvert(maf, "psl", prefix="contigs.genome", **MAFCONVERT_ARGS),
        last_mafconvert(maf, "sam", prefix="contigs.genome", **MAFCONVERT_ARGS),
    )
    for result, name in ((psl, "contigs.genome.psl.gz"), (sam, "contigs.genome.sam.gz")):
        assert result.alignment.path.endswith("/" + name), result.alignment.path
        assert result.index is None and result.stats is None, result
    await assert_gunzipped_md5(psl.alignment, "515d3cff55d159309bedd38f47dd034b", label="psl.gz")
    await assert_gunzipped_md5(sam.alignment, "b5f9324545a03c6bc54566b5ccf5b35f", label="sam.gz")


@env.task
async def test_last_mafconvert_bam_cram() -> None:
    # upstream cases: last/mafconvert "sarscov2 - bam" and "sarscov2 - cram" (the reads' md5, as upstream)
    maf, fasta, fai, sequence_dict = await asyncio.gather(
        fixture(MAF),
        fixture(GENOME + "genome.fasta"),
        fixture(GENOME + "genome.fasta.fai"),
        fixture(GENOME + "genome.dict"),
    )
    bam, cram = await asyncio.gather(
        last_mafconvert(maf, "bam", prefix="contigs.genome", **MAFCONVERT_ARGS),
        last_mafconvert(
            maf, "cram", fasta=fasta, fai=fai, sequence_dict=sequence_dict, prefix="contigs.genome", **MAFCONVERT_ARGS
        ),
    )
    assert bam.alignment.path.endswith("/contigs.genome.bam"), bam.alignment.path
    assert cram.alignment.path.endswith("/contigs.genome.cram"), cram.alignment.path
    await assert_reads_md5(bam.alignment, "5027207cb5c220e7101c85617fd4a0e2", label="bam")
    sam = await samtools_view(bam=cram.alignment, fasta=fasta, fai=fai)
    h = hashlib.md5()
    async with sam.open("rb") as fh:
        for line in bytes(await fh.read()).decode().splitlines():
            h.update(line.split("\t")[9].encode())
    assert h.hexdigest().lstrip("0") == "5027207cb5c220e7101c85617fd4a0e2", f"cram reads md5 {h.hexdigest()}"


@env.task
async def test_last_mafconvert_bcf() -> None:
    # upstream case: last/mafconvert "sarscov2 - bcf" (a nonempty .bcf; upstream snapshots only versions)
    maf, fasta, fai, sequence_dict = await asyncio.gather(
        fixture(MAF),
        fixture(GENOME + "genome.fasta"),
        fixture(GENOME + "genome.fasta.fai"),
        fixture(GENOME + "genome.dict"),
    )
    result = await last_mafconvert(
        maf, "bcf", fasta=fasta, fai=fai, sequence_dict=sequence_dict, prefix="contigs.genome", **MAFCONVERT_ARGS
    )
    assert result.alignment.path.endswith("/contigs.genome.bcf"), result.alignment.path
    async with result.alignment.open("rb") as fh:
        assert len(bytes(await fh.read(1))) > 0, "empty bcf"
    assert result.stats is not None, "no bcftools stats"


@env.task
async def test_last_mafconvert_gff() -> None:
    # upstream case: last/mafconvert "sarscov2 - gff"
    maf, sizes = await asyncio.gather(fixture(MAF), fixture(GENOME + "genome.sizes"))
    result = await last_mafconvert(maf, "gff", sizes=sizes, prefix="contigs.genome", **MAFCONVERT_ARGS)
    assert result.alignment.path.endswith("/contigs.genome.gff.gz"), result.alignment.path
    await assert_gunzipped_md5(result.alignment, "f34ea8242bbc9b4d7eefa04e0fcf8c32", label="gff.gz")


@env.task
async def test_last_dotplot() -> None:
    # upstream cases: last/dotplot "... png" and "... png - filter" (the plot exists) and "... gif" (md5)
    maf, annot_a, annot_b = await asyncio.gather(
        fixture(MAF),
        last_text_file("dummy_annot_a.bed", "MT192765.1\t1000\t1010\n"),
        last_text_file("dummy_annot_b.bed", "NODE_1_length_20973_cov_191.628754\t2000\t2010\n"),
    )
    png, png_filtered, gif = await asyncio.gather(
        last_dotplot(maf, annot_a, annot_b, format="png", prefix="test"),
        last_dotplot(maf, annot_a, annot_b, format="png", filter=True, prefix="test"),
        last_dotplot(maf, format="gif", prefix="test"),
    )
    await check_image(png, "test.png", b"\x89PNG")
    await check_image(png_filtered, "test.png", b"\x89PNG")
    await check_image(gif, "test.gif", b"GIF8")
    await assert_md5(gif, "cec730f1d060a82553cefe2311708bad", label="test.gif")


@env.task
async def test_last_chain() -> None:
    # Not an upstream case: upstream's chained tools read fixtures. Here an index built from the genome
    # must train and align contigs exactly as upstream's fixture index does, and the alignments go on
    # through last-split to PSL and a dot plot.
    fasta, contigs, fixture_index = await asyncio.gather(
        fixture(GENOME + "genome.fasta"), fixture(CONTIGS), last_fixture_index()
    )
    index = await last_lastdb(fasta, prefix="genome")
    built, upstream = await asyncio.gather(
        last_train(contigs, index, prefix="contigs"), last_train(contigs, fixture_index, prefix="contigs")
    )
    for name in ("param_file", "multiqc"):
        a, b = await asyncio.gather(last_file_md5(getattr(built, name)), last_file_md5(getattr(upstream, name)))
        assert a == b, f"train {name}: built index {a} != fixture index {b}"
    aligned = await last_lastal(contigs, index, built.param_file, prefix="contigs.genome", args="-P1")
    split = await last_split(aligned.maf, prefix="contigs.genome.split")
    psl, plot = await asyncio.gather(
        last_mafconvert(split.maf, "psl", prefix="contigs.genome", **MAFCONVERT_ARGS),
        last_dotplot(split.maf, prefix="contigs.genome"),
    )
    async with psl.alignment.open("rb") as fh:
        rows = [line.split("\t") for line in gzip.decompress(bytes(await fh.read())).decode().splitlines()]
    assert rows and all(len(r) == 21 and r[13] == "MT192765.1" for r in rows), rows[:2]
    await check_image(plot, "contigs.genome.png", b"\x89PNG")
    async with split.multiqc.open("rb") as fh:
        table = bytes(await fh.read()).decode().splitlines()
    assert table[0].startswith("Sample\tTotalAlignmentLength"), table
    assert table[1].split("\t")[0] == "contigs.genome.split" and table[1].split("\t")[-4:-2] == ["1", "29829"], table


tests = [
    test_last_lastdb_fasta,
    test_last_lastdb_fastq_gz,
    test_last_train,
    test_last_lastal,
    test_last_lastal_trained,
    test_last_split,
    test_last_mafconvert_psl_sam,
    test_last_mafconvert_bam_cram,
    test_last_mafconvert_bcf,
    test_last_mafconvert_gff,
    test_last_dotplot,
    test_last_chain,
]
