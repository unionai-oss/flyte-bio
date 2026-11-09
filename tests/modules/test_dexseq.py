"""Tests for flyte_bio.modules.dexseq.

The cases and expected values are upstream's dexseq/prepareannotation,
dexseq/count, dexseq/exon and dexseq/dtu tests (modules commit 9a1cc9f,
DEXSeq 1.56.0, HTSeq 2.1.2), with their test configs' ``ext.args``. As
upstream, the count tests count against the annotation flattened (with hash
seed 0) from the same GTF, and without a seed only the GFF's name and its
aggregate genes are checked.

The DEU and DTU fits depend on the OpenBLAS kernel picked for the CPU, so,
as upstream, the exon tests check the IDs, counts and significant calls of
the results, and the DTU test the results rounded to 6 decimals.
"""

import hashlib
import re
from decimal import ROUND_HALF_UP, Context, Decimal

from flyte.io import File

from flyte_bio.modules.dexseq import (
    DexseqResult,
    dexseq_count,
    dexseq_dtu,
    dexseq_exon,
    dexseq_prepareannotation,
)
from tests.framework import assert_md5, env, fixture

GTF = "genomics/homo_sapiens/genome/genome.gtf"
BAM = "genomics/homo_sapiens/illumina/bam/test.rna.paired_end.sorted.bam"
EXON = "delete_me/dexseq/exon/"
DTU = "delete_me/dexseq/dtu/"
EXON_SAMPLES = ["ERR188383", "ERR188428", "ERR188454", "ERR204916"]
AGGREGATE_GENE = re.compile(r'gene_id "[^"]*\+[^"]*"')
EXON_HEADER = [
    "groupID",
    "featureID",
    "exonBaseMean",
    "dispersion",
    "stat",
    "pvalue",
    "padj",
    "YRI",
    "GBR",
    "log2fold_GBR_YRI",
    "countData.1",
    "countData.2",
    "countData.3",
    "countData.4",
]
EXON_SIGNIFICANT = [
    "ENSG00000147403+ENSG00000207165:E024",
    "ENSG00000147403+ENSG00000207165:E032",
    "ENSG00000147403+ENSG00000207165:E033",
]
# A decimal number as Groovy's String.isNumber() (BigDecimal) accepts it.
NUMBER = re.compile(r"[+-]?(\d+(\.\d*)?|\.\d+)([eE][+-]?\d+)?")


async def text(file: File) -> str:
    async with file.open("rb") as fh:
        return bytes(await fh.read()).decode()


async def table(file: File) -> list[list[str]]:
    """``file``'s lines split on tabs, keeping empty fields (Groovy's ``readLines()*.split('\\t', -1)``)."""
    return [line.split("\t") for line in (await text(file)).splitlines()]


def md5(lines: list[str]) -> str:
    return hashlib.md5("\n".join(lines).encode()).hexdigest()


def significant(value: str) -> bool:
    return value != "NA" and float(value) < 0.05


def rounded(value: str) -> str:
    """Java's ``String.format('%.6f', value as double)`` for numbers, else ``value``.

    Java rounds the double's shortest decimal form half up.
    """
    if not NUMBER.fullmatch(value):
        return value
    return format(
        Decimal(repr(float(value))).quantize(Decimal("0.000001"), rounding=ROUND_HALF_UP, context=Context(prec=400)),
        "f",
    )


async def genome_gff() -> File:
    """The test GTF flattened with hash seed 0, as the count tests' setup step."""
    return await dexseq_prepareannotation(await fixture(GTF), prefix="genome", random_seed=0)


@env.task
async def test_dexseq_prepareannotation() -> None:
    # upstream case: dexseq/prepareannotation "homo_sapiens - gtf"
    gff = await dexseq_prepareannotation(await fixture(GTF), prefix="test", random_seed=0)
    assert AGGREGATE_GENE.search(await text(gff)), "no aggregate genes"
    await assert_md5(gff, "e685a7bede4b9edd028f9ace4b7ef69a", label="gff")


@env.task
async def test_dexseq_prepareannotation_no_aggregation() -> None:
    # upstream case: dexseq/prepareannotation "homo_sapiens - gtf - no aggregation"
    gff = await dexseq_prepareannotation(await fixture(GTF), prefix="test", random_seed=0, args="-r no")
    assert not AGGREGATE_GENE.search(await text(gff)), "aggregate genes despite -r no"
    await assert_md5(gff, "30288b2948e496452f33abed5ed6666d", label="gff")


@env.task
async def test_dexseq_prepareannotation_no_seed() -> None:
    # upstream case: dexseq/prepareannotation "homo_sapiens - gtf - no seed"
    gff = await dexseq_prepareannotation(await fixture(GTF), prefix="test")
    assert gff.name == "test.gff", f"unexpected name {gff.name}"
    assert AGGREGATE_GENE.search(await text(gff)), "no aggregate genes"


@env.task
async def test_dexseq_count_paired_end() -> None:
    # upstream case: dexseq/count "homo_sapiens - paired_end - bam"
    counts = await dexseq_count(await fixture(BAM), await genome_gff(), prefix="test", args="--stranded no")
    assert "ENSG00000239435:001\t71" in (await text(counts)).splitlines(), "no ENSG00000239435:001 count of 71"
    await assert_md5(counts, "502433dfddcd443eb3aaab62bdc9391e", label="counts")


@env.task
async def test_dexseq_count_paired_end_reverse() -> None:
    # upstream case: dexseq/count "homo_sapiens - paired_end - bam - reverse stranded"
    counts = await dexseq_count(
        await fixture(BAM), await genome_gff(), prefix="test", args="--stranded reverse --minaqual 10"
    )
    await assert_md5(counts, "aac2e2d89746199fc1f0a54baf0c7cc3", label="counts")


@env.task
async def test_dexseq_count_single_end() -> None:
    # upstream case: dexseq/count "homo_sapiens - single_end - bam" (upstream's paired-end BAM, counted as single-end)
    counts = await dexseq_count(
        await fixture(BAM), await genome_gff(), prefix="test", single_end=True, args="--stranded no"
    )
    assert "ENSG00000239435:001\t132" in (await text(counts)).splitlines(), "no ENSG00000239435:001 count of 132"
    await assert_md5(counts, "3da54cd3b59319df4c37326aa56b231b", label="counts")


async def run_exon(args: str) -> DexseqResult:
    return await dexseq_exon(
        samplesheet=await fixture(EXON + "chrX.samplesheet.csv"),
        counts=[await fixture(f"{EXON}{sample}.txt") for sample in EXON_SAMPLES],
        gff=await fixture(EXON + "chrX.DEXSeq.gff.gz"),
        contrast_variable="condition",
        reference="YRI",
        target="GBR",
        prefix="GBR_vs_YRI",
        args=args,
    )


async def check_exon(result: DexseqResult) -> None:
    """Upstream's snapshot: the header, the IDs and counts, the significant exons and genes."""
    results = await table(result.results)
    qvalues = await table(result.gene_qvalues)
    header = results[0]
    assert header == EXON_HEADER, f"unexpected results header {header}"
    stable = [i for i, column in enumerate(header) if re.fullmatch(r"groupID|featureID|countData\..+", column)]
    padj = header.index("padj")
    stable_md5 = md5(["\t".join(row[i] for i in stable) for row in results])
    assert stable_md5 == "1ab928f072f369098af5646bd13593dc", f"IDs and counts md5 {stable_md5}"
    calls = [f"{row[0]}:{row[1]}" for row in results[1:] if significant(row[padj])]
    assert calls == EXON_SIGNIFICANT, f"significant exons {calls}"
    ids_md5 = md5([row[0] for row in qvalues])
    assert ids_md5 == "5a61555c5923c3100519f4634e387040", f"gene IDs md5 {ids_md5}"
    genes = [row[0] for row in qvalues[1:] if significant(row[1])]
    assert genes == ["ENSG00000147403+ENSG00000207165"], f"significant genes {genes}"
    await assert_md5(result.versions, "469726af854d0c52e1b4dac0968f9bf2", label="versions.yml")


@env.task
async def test_dexseq_exon() -> None:
    # upstream case: dexseq/exon "homo_sapiens - chrX - GBR vs YRI"
    result = await run_exon("")
    await check_exon(result)
    assert result.plots is not None and result.plots.name == "GBR_vs_YRI.dexseq.plots.pdf", "no plots"


@env.task
async def test_dexseq_exon_no_plots() -> None:
    # upstream case: dexseq/exon "homo_sapiens - chrX - GBR vs YRI - no plots"
    result = await run_exon("--n_plot 0")
    await check_exon(result)
    assert result.plots is None, "plots despite --n_plot 0"


@env.task
async def test_dexseq_dtu() -> None:
    # upstream case: dexseq/dtu "homo_sapiens - chrX - GBR vs YRI"
    result = await dexseq_dtu(
        samplesheet=await fixture(DTU + "chrX.samplesheet.csv"),
        counts=await fixture(DTU + "chrX.transcript_counts.tsv"),
        contrast_variable="condition",
        reference="YRI",
        target="GBR",
        prefix="GBR_vs_YRI",
        args="--feature_id_col tx",
    )
    results = await table(result.results)
    qvalues = await table(result.gene_qvalues)
    assert "log2fold_GBR_YRI" in results[0], f"unexpected results header {results[0]}"
    n_significant = sum(significant(row[1]) for row in qvalues[1:])
    assert n_significant == 23, f"{n_significant} significant genes, expected 23"
    for rows, expected, label in [
        (results, "a1d1eb62c416e082965e1a1e672f1e70", "results"),
        (qvalues, "e5787d3fd004cd02f812dd41358374cf", "gene q-values"),
    ]:
        actual = md5(["\t".join(rounded(v) for v in row) for row in rows])
        assert actual == expected, f"rounded {label} md5 {actual}, expected {expected}"
    await assert_md5(result.versions, "415dfc79a4c7dde0f048b1ea6c1dc1c2", label="versions.yml")


tests = [
    test_dexseq_prepareannotation,
    test_dexseq_prepareannotation_no_aggregation,
    test_dexseq_prepareannotation_no_seed,
    test_dexseq_count_paired_end,
    test_dexseq_count_paired_end_reverse,
    test_dexseq_count_single_end,
    test_dexseq_exon,
    test_dexseq_exon_no_plots,
    test_dexseq_dtu,
]
