"""Tests for flyte_bio.modules.annotsv.

Upstream's annotsv/annotsv and annotsv/installannotations tests (modules
commit 9a1cc9f, AnnotSV 3.5.10) are all stubs: they only check that the
outputs are created. Two of them check that installannotations refuses one
version without the other; those are ported as checks of the wrapper.

:func:`test_annotsv_annotate` is not an upstream case: it installs the real
annotations (human 3.5; Exomiser 2406, which the container's script fails to
fetch) and annotates upstream's stub inputs, then reports the table with knotAnnotSV, as a rare-disease
pipeline chains them. AnnotSV's ranking depends on the annotation release,
so it checks the outputs' structure: every input SV of at least AnnotSV's
minimum size (50 bp) is annotated (a ``full`` line) or counted as
unannotated, split lines belong to annotated SVs, and the HTML report is
written. Upstream's test config adds ``-vcf 1``; it is left out because the
container can't write the VCF (its variantconvert fails under its pandas 3;
see the module notes), and the module now fails on it. The install downloads 5.3 GB (about two hours from the
annotation server); its result is cached.
"""

import asyncio
import gzip
import tempfile
from pathlib import Path

from flyte.io import File

from flyte_bio.modules.annotsv import annotsv_annotsv, annotsv_installannotations
from flyte_bio.modules.knotannotsv import knotannotsv
from tests.framework import env, fixture

ILLUMINA = "genomics/homo_sapiens/illumina/"
CALLS = ILLUMINA + "gatk/haplotypecaller_calls/"


async def annotsv_text(file: File) -> str:
    async with file.open("rb") as fh:
        data = bytes(await fh.read())
    return (gzip.decompress(data) if data[:2] == b"\x1f\x8b" else data).decode()


async def annotsv_text_file(name: str, text: str) -> File:
    path = Path(tempfile.mkdtemp()) / name
    path.write_text(text)
    return await File.from_local(str(path))


@env.task
async def test_annotsv_installannotations_versions() -> None:
    # upstream cases: annotsv/installannotations "annotsv version empty" and "exomiser version empty" (both fail)
    for versions in ((None, "2406"), ("3.5", None)):
        try:
            await annotsv_installannotations(*versions)
        except ValueError:
            continue
        raise AssertionError(f"annotsv_installannotations accepted versions {versions}")


@env.task
async def test_annotsv_annotate() -> None:
    # Not an upstream case (upstream's are stubs): upstream's stub inputs, annotated against the real annotations.
    annotations, sv_vcf, sv_tbi, small, false_positives, genes, transcripts = await asyncio.gather(
        annotsv_installannotations("3.5", "2406"),
        fixture(ILLUMINA + "vcf/sv_query.vcf.gz"),
        fixture(ILLUMINA + "vcf/sv_query.vcf.gz.tbi"),
        fixture(CALLS + "test2_haplotc.vcf.gz"),
        fixture(CALLS + "test_haplotcaller.cnn.vcf.gz"),
        # As upstream's collectFile, without newlines.
        annotsv_text_file("gene_candidates.txt", "GENE1GENE2GENE3"),
        annotsv_text_file("gene_transcripts.txt", "GENE1 GENE2 GENE3"),
    )
    result = await annotsv_annotsv(
        sv_vcf,
        annotations,
        prefix="test",
        sv_vcf_index=sv_tbi,
        candidate_small_variants=small,
        candidate_genes=genes,
        false_positive_snv=false_positives,
        gene_transcripts=transcripts,
    )
    assert result.tsv.path.endswith("/test.tsv"), f"tsv path {result.tsv.path}"

    records = [line.split("\t") for line in (await annotsv_text(sv_vcf)).splitlines() if line and line[0] != "#"]
    all_ids = {r[2] for r in records}
    # AnnotSV skips SVs under its minimum size (-SVminSize, 50 bp), measured on sequence alleles.
    small = {r[2] for r in records if not r[4].startswith("<") and abs(len(r[4]) - len(r[3])) < 50}
    sv_ids = all_ids - small
    rows = [line.split("\t") for line in (await annotsv_text(result.tsv)).splitlines()]
    header, body = rows[0], rows[1:]
    assert header[0] == "AnnotSV_ID" and "Annotation_mode" in header and "ID" in header, f"tsv header {header[:10]}"
    mode, sv_id = header.index("Annotation_mode"), header.index("ID")
    annotated = {r[sv_id] for r in body if r[mode] == "full"}
    split = {r[sv_id] for r in body if r[mode] == "split"}
    unannotated = 0
    if result.unannotated_tsv is not None:
        unannotated = len([line for line in (await annotsv_text(result.unannotated_tsv)).splitlines()[1:] if line])
    assert annotated, f"no annotated SVs; modes {sorted({r[mode] for r in body})}, {len(body)} rows"
    assert split <= annotated, f"split lines of unannotated SVs: {sorted(split - annotated)}"
    assert annotated <= all_ids, f"annotated SVs not in the input: {sorted(annotated - all_ids)}"
    missing = sorted(sv_ids - annotated)
    assert len(annotated) + unannotated >= len(sv_ids), f"SVs neither annotated nor unannotated: {missing}"

    assert result.vcf is None, f"unexpected VCF {result.vcf.path}"

    report = await knotannotsv(result.tsv, prefix="test", args="--genomeBuild hg38")
    html = await annotsv_text(report)
    assert "<html" in html.lower(), f"report {report.path} is not HTML: {html[:200]!r}"
    # The report names each SV by its AnnotSV_ID.
    reported = {r[0] for r in body if r[mode] == "full"}
    absent = sorted(i for i in reported if i not in html)
    assert not absent, f"report {report.path} ({len(html)} chars) lacks annotated SVs {absent[:5]}"


tests = [test_annotsv_installannotations_versions, test_annotsv_annotate]
