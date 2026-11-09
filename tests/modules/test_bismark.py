"""Tests for flyte_bio.modules.bismark.

The cases and expected values are upstream's bismark/genomepreparation,
bismark/align, bismark/deduplicate, bismark/methylationextractor,
bismark/coverage2cytosine, bismark/report and bismark/summary tests (modules
commit 9a1cc9f, Bismark 3.1.0). As upstream, the tests chain the tools: each
builds its index with genome preparation (with the aligner options the
alignment uses, as upstream's test config passes them to both), and the
report and summary tests run alignment, deduplication and extraction first.
Gzipped outputs are compared decompressed, as nf-test hashes them; BAMs by
their reads' md5.

Upstream's tests run with 4 CPUs, the alignment's default here: its
threading cases check the per-instance ``-p`` it records in its report.
"""

from flyte.io import Dir, File

from flyte_bio.modules.bismark import (
    BismarkAlignResult,
    BismarkMethylationExtractorResult,
    bismark_align,
    bismark_coverage2cytosine,
    bismark_deduplicate,
    bismark_genomepreparation,
    bismark_methylationextractor,
    bismark_report,
    bismark_summary,
)
from tests.framework import assert_gunzipped_md5, assert_md5, assert_reads_md5, env, fixture

GENOME = "genomics/sarscov2/genome/genome.fasta"
FASTQ = "genomics/sarscov2/illumina/fastq/"
METHYLATED_BAM = "genomics/sarscov2/illumina/bam/test.paired_end.methylated.bam"
UNIQUE_HITS = "Number of alignments with a unique best hit from the different alignments:\t"


async def bismark_text(file: File) -> str:
    async with file.open("rb") as fh:
        return bytes(await fh.read()).decode()


async def bismark_index(args: str = "") -> Dir:
    return await bismark_genomepreparation(await fixture(GENOME), args=args)


async def bismark_reads(paired: bool) -> list[File]:
    mates = ["1", "2"] if paired else ["1"]
    return [await fixture(f"{FASTQ}test.methylated_{m}.fastq.gz") for m in mates]


async def bismark_check_align(result: BismarkAlignResult, reads_md5: str, unique_hits: str, found: bool) -> None:
    await assert_reads_md5(result.bam, reads_md5, label=result.bam.name)
    lines = (await bismark_text(result.report)).splitlines()
    wanted = UNIQUE_HITS + unique_hits
    assert (wanted in lines) == found, f"{result.report.name}: {wanted!r} {'missing' if found else 'present'}"
    assert result.unmapped == [], f"unexpected unmapped reads {[f.name for f in result.unmapped]}"


@env.task
async def test_bismark_genomepreparation() -> None:
    # upstream case: "sarscov2 genome [fasta]"
    index = await bismark_index()
    root = index.path.rstrip("/") + "/"
    files = {f.path.removeprefix(root): f async for f in index.walk()}
    expected = {
        "Bisulfite_Genome/CT_conversion/BS_CT.1.bt2": "10a8854e5eeb95629f46d24b0d8da70c",
        "Bisulfite_Genome/CT_conversion/BS_CT.2.bt2": "bc5d407e6fce3e1cdb4b23ab2a20f707",
        "Bisulfite_Genome/CT_conversion/BS_CT.3.bt2": "4ed93abba181d8dfab2e303e33114777",
        "Bisulfite_Genome/CT_conversion/BS_CT.4.bt2": "6c4f549575e5882ff3a97ae10ae2e7be",
        "Bisulfite_Genome/CT_conversion/BS_CT.rev.1.bt2": "d37d586b41b6a36efb25839961d242cb",
        "Bisulfite_Genome/CT_conversion/BS_CT.rev.2.bt2": "48bf400de85cd7324fbf6bf38fe2dc95",
        "Bisulfite_Genome/CT_conversion/genome_mfa.CT_conversion.fa": "903b9f357eea4a5f36e21e78e0fe1dfa",
        "Bisulfite_Genome/GA_conversion/BS_GA.1.bt2": "a6fa4068ed10872568f32568c66cb600",
        "Bisulfite_Genome/GA_conversion/BS_GA.2.bt2": "fc5b634e2f5137801e3801d62ad05f74",
        "Bisulfite_Genome/GA_conversion/BS_GA.3.bt2": "4ed93abba181d8dfab2e303e33114777",
        "Bisulfite_Genome/GA_conversion/BS_GA.4.bt2": "b1e855685ed634daa7df8e180c9b2fc4",
        "Bisulfite_Genome/GA_conversion/BS_GA.rev.1.bt2": "17f51e65f44477962226a83b369b2bc4",
        "Bisulfite_Genome/GA_conversion/BS_GA.rev.2.bt2": "d36044c5dfe37af0411f60d446585bc6",
        "Bisulfite_Genome/GA_conversion/genome_mfa.GA_conversion.fa": "e529375ec2039112161465b5169fcd59",
        "genome.fasta": "6e9fe4042a72f2345f644f239272b7e6",
    }
    assert sorted(files) == sorted(expected), f"index files differ: {sorted(files)}"
    for name, md5 in expected.items():
        await assert_md5(files[name], md5, label=name)


@env.task
async def test_bismark_align_bowtie2_single_end() -> None:
    # upstream case: "bowtie2 | single-end | sarscov2 genome [fasta]"
    result = await bismark_align(await bismark_reads(False), await bismark_index("--bowtie2"), args="--bowtie2")
    await bismark_check_align(result, "936c0d5ce713130113e99e09a5b53afd", "5009", True)


@env.task
async def test_bismark_align_bowtie2_paired_end() -> None:
    # upstream case: "bowtie2 | paired-end | sarscov2 genome [fasta]"
    result = await bismark_align(await bismark_reads(True), await bismark_index("--bowtie2"), args="--bowtie2")
    await bismark_check_align(result, "458c5350c3f96a8dde710985bfb62827", "5009", False)


@env.task
async def test_bismark_align_hisat2_single_end() -> None:
    # upstream case: "hisat2 | single-end | sarscov2 genome [fasta]"
    result = await bismark_align(await bismark_reads(False), await bismark_index("--hisat2"), args="--hisat2")
    await bismark_check_align(result, "936c0d5ce713130113e99e09a5b53afd", "5009", True)


@env.task
async def test_bismark_align_hisat2_paired_end() -> None:
    # upstream case: "hisat2 | paired-end | sarscov2 genome [fasta]"
    result = await bismark_align(await bismark_reads(True), await bismark_index("--hisat2"), args="--hisat2")
    await bismark_check_align(result, "fb9d284aec4b2c727af3c5ffd77c6381", "5009", False)


@env.task
async def test_bismark_align_minimap2_single_end() -> None:
    # upstream case: "minimap2 | single-end | sarscov2 genome [fasta]"
    result = await bismark_align(await bismark_reads(False), await bismark_index("--minimap2"), args="--minimap2")
    await bismark_check_align(result, "d354e6c6d6a06e36630590de8f27031a", "5047", True)


async def bismark_report_options(paired: bool, index_args: str, args: str) -> str:
    result = await bismark_align(await bismark_reads(paired), await bismark_index(index_args), args=args)
    return await bismark_text(result.report)


@env.task
async def test_bismark_align_bowtie2_threads() -> None:
    # upstream case: "bowtie2 | intra-instance -p threading | paired-end | 4 cpus"
    report = await bismark_report_options(True, "--bowtie2", "--bowtie2")
    assert "-p 2 --reorder" in report, "directional alignment with 4 CPUs should run 2 instances with -p 2"


@env.task
async def test_bismark_align_bowtie2_non_directional_threads() -> None:
    # upstream case: "bowtie2 | non_directional | low-cpu omits -p (4 cpus / 4 instances < 2)"
    report = await bismark_report_options(False, "--bowtie2", "--bowtie2 --non_directional")
    assert "-p " not in report, "non-directional alignment with 4 CPUs should omit -p"


@env.task
async def test_bismark_align_hisat2_threads() -> None:
    # upstream case: "hisat2 | intra-instance -p threading | paired-end | 4 cpus"
    report = await bismark_report_options(True, "--hisat2", "--hisat2")
    assert "-p 2 --reorder" in report, "directional HISAT2 alignment with 4 CPUs should run -p 2"


@env.task
async def test_bismark_align_bowtie2_combined_index_threads() -> None:
    # upstream case: "bowtie2 | combined_index | intra-instance -p threading | single-end | 4 cpus"
    report = await bismark_report_options(False, "--bowtie2 --combined_genome", "--bowtie2 --combined_index")
    assert "-p 4 --reorder" in report, "combined-index alignment with 4 CPUs should run one instance with -p 4"


@env.task
async def test_bismark_deduplicate() -> None:
    # upstream case: "sarscov2 paired-end methylated [bam]"
    result = await bismark_deduplicate(await fixture(METHYLATED_BAM), single_end=False)
    await assert_reads_md5(result.bam, "8036b4b920354438ddd1871c88f5088", label=result.bam.name)
    assert result.report.name == "test.paired_end.methylated.deduplication_report.txt", result.report.name
    await assert_md5(result.report, "86aaf87bac9ccdc3cc83628f20b902dc", label=result.report.name)


async def bismark_extract_methylated() -> BismarkMethylationExtractorResult:
    return await bismark_methylationextractor(await fixture(METHYLATED_BAM), await bismark_index(), single_end=False)


@env.task
async def test_bismark_methylationextractor() -> None:
    # upstream case: "sarscov2 paired-end methylated [bam] | sarscov2 genome [fasta]"
    result = await bismark_extract_methylated()
    name = "test.paired_end.methylated"
    await assert_gunzipped_md5(result.bedgraph, "3e461daa0b89ee4eeb76e336ad6aad7a", label=result.bedgraph.name)
    await assert_gunzipped_md5(result.coverage, "1524da1472ed1ad0c0663a9676c142e0", label=result.coverage.name)
    assert result.coverage.name == f"{name}.bismark.cov.gz", result.coverage.name
    await assert_md5(result.mbias, "0b100924d46b3c35115f1206f34c4a59", label=result.mbias.name)
    await assert_md5(result.report, "e09d4e72fb8752df4ad46a9ed2894943", label=result.report.name)
    calls = {
        f"CHG_OB_{name}.txt.gz": "dfd299b7d7661fc6e331bd4fe60c7c07",
        f"CHG_OT_{name}.txt.gz": "9a85dc0b1650f06ecda55643e081a66e",
        f"CHH_OB_{name}.txt.gz": "20d937df40f2fa9a1bd878873335d3c3",
        f"CHH_OT_{name}.txt.gz": "d3705b17290e017eb1244808dae2ef23",
        f"CpG_OB_{name}.txt.gz": "4c3619590c86e408dd60d0075fe08d1a",
        f"CpG_OT_{name}.txt.gz": "90b9b96895fcdae9785a7c35d829d379",
    }
    found = {f.name: f for f in result.methylation_calls}
    assert sorted(found) == sorted(calls), f"methylation calls differ: {sorted(found)}"
    for call_name, md5 in calls.items():
        await assert_gunzipped_md5(found[call_name], md5, label=call_name)


@env.task
async def test_bismark_coverage2cytosine() -> None:
    # upstream case: "sarscov2 genome [fasta] | sarscov2 paired-end methylated [bam]"
    extracted = await bismark_extract_methylated()
    result = await bismark_coverage2cytosine(extracted.coverage, await bismark_index(), prefix="test")
    assert result.coverage == [], f"unexpected coverage files {[f.name for f in result.coverage]}"
    assert result.report.name == "test.CpG_report.txt.gz", result.report.name
    await assert_gunzipped_md5(result.report, "0b83c760c7b3df21a320f054bbb865fb", label=result.report.name)
    await assert_md5(result.summary, "5ed486bfa844a5e476a7ab054a7a8d6b", label=result.summary.name)


async def bismark_single_end_chain() -> tuple[BismarkAlignResult, File, BismarkMethylationExtractorResult]:
    """Upstream's report and summary setup: align single-end reads, deduplicate and extract."""
    index = await bismark_index()
    aligned = await bismark_align(await bismark_reads(False), index)
    deduplicated = await bismark_deduplicate(aligned.bam, single_end=True)
    extracted = await bismark_methylationextractor(aligned.bam, index, single_end=True)
    return aligned, deduplicated.report, extracted


@env.task
async def test_bismark_report() -> None:
    # upstream case: "sarscov2 genome [fasta] | sarscov2 single-end methylated [fastq]"
    aligned, dedup_report, extracted = await bismark_single_end_chain()
    results = await bismark_report(aligned.report, dedup_report, extracted.report, extracted.mbias)
    names = sorted([f.name async for f in results.walk()])
    assert names == ["test.methylated_1_bismark_bt2_SE_report.html"], f"bismark2report wrote {names}"


@env.task
async def test_bismark_summary() -> None:
    # upstream case: "sarscov2 genome [fasta] | sarscov2 single-end methylated [fastq]"
    aligned, dedup_report, extracted = await bismark_single_end_chain()
    result = await bismark_summary(
        [aligned.bam], [aligned.report], [dedup_report], [extracted.report], [extracted.mbias]
    )
    names = sorted([f.name async for f in result.results.walk()])
    assert names == ["bismark_summary_report.html", "bismark_summary_report.txt"], f"bismark2summary wrote {names}"
    assert "CpGs" in (await bismark_text(result.txt)).splitlines()[0], "summary table's header lacks CpGs"
    assert "</html>" in (await bismark_text(result.html)).splitlines()[-1], "summary HTML doesn't end with </html>"


tests = [
    test_bismark_genomepreparation,
    test_bismark_align_bowtie2_single_end,
    test_bismark_align_bowtie2_paired_end,
    test_bismark_align_hisat2_single_end,
    test_bismark_align_hisat2_paired_end,
    test_bismark_align_minimap2_single_end,
    test_bismark_align_bowtie2_threads,
    test_bismark_align_bowtie2_non_directional_threads,
    test_bismark_align_hisat2_threads,
    test_bismark_align_bowtie2_combined_index_threads,
    test_bismark_deduplicate,
    test_bismark_methylationextractor,
    test_bismark_coverage2cytosine,
    test_bismark_report,
    test_bismark_summary,
]
