"""Tests for flyte_bio.modules.cnvkit.

The cases and expected values are upstream's cnvkit/target and cnvkit/batch
tests (modules commit 9a1cc9f, CNVkit 0.9.14), with each batch case's
``ext.args``. As upstream, outputs that aren't byte-stable (the segments, and
for most cases the target/antitarget BEDs) are checked by name only; the rest
by md5.
"""

import asyncio

from flyte.io import File

from flyte_bio.modules.cnvkit import CnvkitBatchResult, cnvkit_batch, cnvkit_target
from tests.framework import assert_md5, env, fixture

HS = "genomics/homo_sapiens/"
SC2 = "genomics/sarscov2/"
CHR21 = HS + "genome/chr21/sequence/"


@env.task
async def test_cnvkit_target() -> None:
    # upstream cases: cnvkit/target "test-cnvkit-target" and "test-cnvkit-target-with-gff3"
    baits, gff3 = await asyncio.gather(fixture(HS + "genome/genome.bed"), fixture(HS + "genome/genome.gff3"))
    plain, annotated = await asyncio.gather(
        cnvkit_target(baits, prefix="test"), cnvkit_target(baits, prefix="test", annotation=gff3)
    )
    for bed, label in ((plain, "target"), (annotated, "target with gff3")):
        assert bed.name == "test.bed", f"{label}: {bed.name}"
        await assert_md5(bed, "c206fd9dc2a001ff3ca3e88e5a8e9d1b", label=label)


# Expected outputs per category: name -> md5, or None where upstream records only the name.
Expected = dict[str, dict[str, str | None]]


async def check_batch(result: CnvkitBatchResult, expected: Expected, label: str) -> None:
    for category in ("bed", "cnn", "cnr", "cns", "pdf", "png"):
        files: list[File] = getattr(result, category)
        want = expected.get(category, {})
        names = [f.name for f in files]
        assert names == sorted(want), f"{label}: {category} files {names}, expected {sorted(want)}"
        for f in files:
            if want[f.name] is not None:
                await assert_md5(f, want[f.name], label=f"{label} {f.name}")


def segments(sample: str) -> dict[str, str | None]:
    return {f"{sample}.bintest.cns": None, f"{sample}.call.cns": None, f"{sample}.cns": None}


# Upstream's antitarget coverage md5 for every one of these BAMs.
EMPTY_ANTITARGET = "203caf8cef6935bb50b4138097955cb8"

HYBRID: Expected = {
    "bed": {"baits.antitarget.bed": None, "baits.target.bed": None},
    "cnn": {
        "reference.cnn": "ce4ef720a8ca89976842c31bd41977fa",
        "test.paired_end.sorted.antitargetcoverage.cnn": EMPTY_ANTITARGET,
        "test.paired_end.sorted.targetcoverage.cnn": "3fe80b6013ffc3e9968345e810158215",
        "test.single_end.sorted.antitargetcoverage.cnn": EMPTY_ANTITARGET,
        "test.single_end.sorted.targetcoverage.cnn": "aa8a018b1d4d1e688c9f9f6ae01bf4d7",
    },
    "cnr": {"test.paired_end.sorted.cnr": "7e37d73ab604dbc3fe4ebb56aca9bdc3"},
    "cns": segments("test.paired_end.sorted"),
}

WGS: Expected = {
    "bed": {"genome.antitarget.bed": None, "genome.bed": None, "genome.target.bed": None},
    "cnn": {
        "reference.cnn": "005e2e68209ee787f1af4840a10f1fe0",
        "test.paired_end.sorted.antitargetcoverage.cnn": EMPTY_ANTITARGET,
        "test.paired_end.sorted.targetcoverage.cnn": "ff526714696aa49bdc1dc8d00d965266",
        "test2.paired_end.sorted.antitargetcoverage.cnn": EMPTY_ANTITARGET,
        "test2.paired_end.sorted.targetcoverage.cnn": "6ae6b3fce7299eedca6133d911c38fe1",
    },
    "cnr": {"test2.paired_end.sorted.cnr": "80318d06c6b095945a0fb0e85e887cbc"},
    "cns": segments("test2.paired_end.sorted"),
}

TUMOUR_ONLY: Expected = {
    "bed": {
        "reference_chr21.antitarget-tmp.bed": "3d4d20f9f23b39970865d29ef239d20b",
        "reference_chr21.target-tmp.bed": "657b25dbda8516624efa8cb2cf3716ca",
    },
    "cnn": {
        "test2.paired_end.recalibrated.sorted.antitargetcoverage.cnn": "067115082c4af4b64d58c0dc3a3642e4",
        "test2.paired_end.recalibrated.sorted.targetcoverage.cnn": "e6d0190c1c37ce6e41f76ca5b24ccca3",
    },
    "cnr": {"test2.paired_end.recalibrated.sorted.cnr": "d9bdb71ce807051369577ee7f807a40c"},
    "cns": segments("test2.paired_end.recalibrated.sorted"),
}

TUMOUR_ONLY_FLAT: Expected = {
    "bed": {"genome.antitarget.bed": None, "genome.target.bed": None},
    "cnn": {
        "reference.cnn": "baada08f38a3e1fca613886a9cd2c557",
        "test2.paired_end.sorted.antitargetcoverage.cnn": EMPTY_ANTITARGET,
        "test2.paired_end.sorted.targetcoverage.cnn": "8268feb252c0f5c47b57e4d8ce51e652",
    },
    "cnr": {"test2.paired_end.sorted.cnr": "2b63cfc469a4744e736a31d1112bea1f"},
    "cns": segments("test2.paired_end.sorted"),
}

GERMLINE_CRAM: Expected = {
    "bed": {"multi_intervals.antitarget.bed": None, "multi_intervals.target.bed": None},
    "cnn": {
        "reference.cnn": "a09ee4be5dda1cf0f68073bdb3aad8ec",
        "test2.paired_end.recalibrated.sorted.antitargetcoverage.cnn": "067115082c4af4b64d58c0dc3a3642e4",
        "test2.paired_end.recalibrated.sorted.targetcoverage.cnn": "b4a49faf170e436ec32dcc21ccc3ce8f",
    },
    "cnr": {"test2.paired_end.recalibrated.sorted.cnr": "3b4fc0cc73be78f978cfe2422470753e"},
    "cns": segments("test2.paired_end.recalibrated.sorted"),
}

PANEL_OF_NORMALS: Expected = {
    "bed": {"genome.antitarget.bed": None, "genome.bed": None, "genome.target.bed": None},
    "cnn": {
        "panel_of_normals.cnn": "ae12d3c464e55e784a83fb24fef9bb84",
        "test.paired_end.sorted.antitargetcoverage.cnn": EMPTY_ANTITARGET,
        "test.paired_end.sorted.targetcoverage.cnn": "ff526714696aa49bdc1dc8d00d965266",
        "test2.paired_end.sorted.antitargetcoverage.cnn": EMPTY_ANTITARGET,
        "test2.paired_end.sorted.targetcoverage.cnn": "6ae6b3fce7299eedca6133d911c38fe1",
    },
}


async def indexed(path: str, index: str = ".bai") -> tuple[File, File]:
    return await asyncio.gather(fixture(path), fixture(path + index))


@env.task
async def test_cnvkit_batch_hybrid() -> None:
    # upstream cases: cnvkit/batch "cnvkit batch hybrid mode - bam" (--output-reference reference.cnn) and
    # "cnvkit batch germline hybrid mode - bam" (no args; the same inputs, so upstream's same outputs)
    (tumor, tumor_bai), (normal, normal_bai), fasta, baits = await asyncio.gather(
        indexed(SC2 + "illumina/bam/test.paired_end.sorted.bam"),
        indexed(SC2 + "illumina/bam/test.single_end.sorted.bam"),
        fixture(SC2 + "genome/genome.fasta"),
        fixture(SC2 + "genome/bed/baits.bed"),
    )
    common = {
        "tumor": tumor,
        "tumor_index": tumor_bai,
        "normals": [normal],
        "normal_indices": [normal_bai],
        "fasta": fasta,
        "targets": baits,
    }
    hybrid, germline_hybrid = await asyncio.gather(
        cnvkit_batch(**common, args="--output-reference reference.cnn"), cnvkit_batch(**common)
    )
    await check_batch(hybrid, HYBRID, "hybrid - bam")
    await check_batch(germline_hybrid, HYBRID, "germline hybrid - bam")


@env.task
async def test_cnvkit_batch_wgs() -> None:
    # upstream cases: cnvkit/batch "cnvkit batch wgs - bam" and "cnvkit batch wgs - cram"
    # (--output-reference reference.cnn --method wgs)
    (tumor, tumor_bai), (normal, normal_bai), tumor_cram, normal_cram, fasta, fai = await asyncio.gather(
        indexed(HS + "illumina/bam/test2.paired_end.sorted.bam"),
        indexed(HS + "illumina/bam/test.paired_end.sorted.bam"),
        fixture(HS + "illumina/cram/test2.paired_end.sorted.cram"),
        fixture(HS + "illumina/cram/test.paired_end.sorted.cram"),
        fixture(HS + "genome/genome.fasta"),
        fixture(HS + "genome/genome.fasta.fai"),
    )
    args = "--output-reference reference.cnn --method wgs"
    from_bam, from_cram = await asyncio.gather(
        cnvkit_batch(
            tumor=tumor, tumor_index=tumor_bai, normals=[normal], normal_indices=[normal_bai], fasta=fasta, args=args
        ),
        cnvkit_batch(tumor=tumor_cram, normals=[normal_cram], fasta=fasta, fai=fai, args=args),
    )
    await check_batch(from_bam, WGS, "wgs - bam")
    await check_batch(from_cram, WGS, "wgs - cram")


@env.task
async def test_cnvkit_batch_tumour_only() -> None:
    # upstream cases: cnvkit/batch "cnvkit batch tumouronly mode - bam", "... - cram" (with a reference) and
    # "cnvkit batch tumour only mode, no reference - bam" (--method wgs)
    (bam, bai), cram, (plain, plain_bai), fasta, fai, reference, targets = await asyncio.gather(
        indexed(HS + "illumina/bam/test2.paired_end.recalibrated.sorted.bam"),
        fixture(HS + "illumina/cram/test2.paired_end.recalibrated.sorted.cram"),
        indexed(HS + "illumina/bam/test2.paired_end.sorted.bam"),
        fixture(CHR21 + "genome.fasta"),
        fixture(CHR21 + "genome.fasta.fai"),
        fixture(CHR21 + "reference_chr21.cnn"),
        fixture(HS + "genome/genome.bed"),
    )
    args = "--method wgs"
    from_bam, from_cram, flat = await asyncio.gather(
        cnvkit_batch(tumor=bam, tumor_index=bai, reference=reference, args=args),
        cnvkit_batch(tumor=cram, fasta=fasta, fai=fai, reference=reference, args=args),
        cnvkit_batch(tumor=plain, tumor_index=plain_bai, targets=targets, args=args),
    )
    await check_batch(from_bam, TUMOUR_ONLY, "tumouronly - bam")
    await check_batch(from_cram, TUMOUR_ONLY, "tumouronly - cram")
    await check_batch(flat, TUMOUR_ONLY_FLAT, "tumour only, no reference - bam")


@env.task
async def test_cnvkit_batch_germline() -> None:
    # upstream case: cnvkit/batch "cnvkit batch germline mode - cram" (no args)
    cram, fasta, fai, targets = await asyncio.gather(
        fixture(HS + "illumina/cram/test2.paired_end.recalibrated.sorted.cram"),
        fixture(CHR21 + "genome.fasta"),
        fixture(CHR21 + "genome.fasta.fai"),
        fixture(CHR21 + "multi_intervals.bed"),
    )
    result = await cnvkit_batch(normals=[cram], fasta=fasta, fai=fai, targets=targets)
    await check_batch(result, GERMLINE_CRAM, "germline - cram")


@env.task
async def test_cnvkit_batch_panel_of_normals() -> None:
    # upstream case: cnvkit/batch "cnvkit batch pon mode - bam" (--method wgs --output-reference panel_of_normals.cnn)
    (one, one_bai), (two, two_bai), fasta = await asyncio.gather(
        indexed(HS + "illumina/bam/test.paired_end.sorted.bam"),
        indexed(HS + "illumina/bam/test2.paired_end.sorted.bam"),
        fixture(HS + "genome/genome.fasta"),
    )
    result = await cnvkit_batch(
        normals=[one, two],
        normal_indices=[one_bai, two_bai],
        fasta=fasta,
        panel_of_normals=True,
        args="--method wgs --output-reference panel_of_normals.cnn",
    )
    await check_batch(result, PANEL_OF_NORMALS, "pon - bam")


tests = [
    test_cnvkit_target,
    test_cnvkit_batch_hybrid,
    test_cnvkit_batch_wgs,
    test_cnvkit_batch_tumour_only,
    test_cnvkit_batch_germline,
    test_cnvkit_batch_panel_of_normals,
]
