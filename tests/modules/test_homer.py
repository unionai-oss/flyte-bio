"""Tests for flyte_bio.modules.homer.

The cases and expected md5s are upstream's homer/maketagdirectory,
homer/findpeaks, homer/makeucscfile, homer/pos2bed and homer/annotatepeaks
tests (modules commit 9a1cc9f, HOMER 4.11), with the ``args`` each test's
config sets. As upstream, the downstream tests run on the tag directory of
the two sarscov2 BEDs (``-format bed``); HOMER finds no peaks in it, so the
peak files are headers only, which record the tag directory's name and the
command line. Upstream treats the bedGraph and the annotation table as
unstable and snapshots only their names; here their content is checked
instead (the track line, the annotation header and its rows).
"""

import gzip
import tempfile
import zipfile
from pathlib import Path

from flyte.io import Dir, File

from flyte_bio.modules.homer import (
    HomerTagDirectory,
    homer_annotatepeaks,
    homer_findpeaks,
    homer_maketagdirectory,
    homer_makeucscfile,
    homer_pos2bed,
)
from tests.framework import assert_md5, env, fixture

GENOME = "genomics/sarscov2/genome/"
BAM = "genomics/sarscov2/illumina/bam/"
UNIQMAP_ZIP = "https://raw.githubusercontent.com/nf-core/test-datasets/nascent/reference/uniqmap.GRCh38_chr21.50nt.zip"


async def text(file: File) -> str:
    async with file.open("rb") as fh:
        return bytes(await fh.read()).decode()


async def bed_tagdir() -> HomerTagDirectory:
    """The tag directory of upstream's setup step: the two sarscov2 BEDs, ``-format bed``."""
    beds = [await fixture(GENOME + "bed/test.bed"), await fixture(GENOME + "bed/test2.bed")]
    return await homer_maketagdirectory(beds, await fixture(GENOME + "genome.fasta"), prefix="test", args="-format bed")


async def check_tagdir(result: HomerTagDirectory, md5s: dict[str, str]) -> None:
    assert result.tagdir.name == "test_tagdir", result.tagdir.name
    files = {f.name: f for f in await result.tagdir.list_files()}
    assert sorted(files) == sorted(md5s), f"tag directory holds {sorted(files)}"
    for name, md5 in md5s.items():
        await assert_md5(files[name], md5, label=name)
    await assert_md5(result.taginfo, md5s["tagInfo.txt"], label="taginfo")


@env.task
async def test_homer_maketagdirectory_bed() -> None:
    # upstream case: homer/maketagdirectory "sarscov2 - [bed, bed], fasta"
    await check_tagdir(
        await bed_tagdir(),
        {
            "MT192765.1.tags.tsv": "e29522171ca2169b57396495f8b97485",
            "tagAutocorrelation.txt": "62b107c4971b94126fb89a0bc2800455",
            "tagCountDistribution.txt": "fd4ee7ce7c5dfd7c9d739534b8180578",
            "tagInfo.txt": "c9bb2ca53bb101d74c1ec92d2b0ad26e",
            "tagLengthDistribution.txt": "e5aa2b9843ca9c04ace297280aed6af4",
        },
    )


@env.task
async def test_homer_maketagdirectory_bam() -> None:
    # upstream case: homer/maketagdirectory "sarscov2 - [bam, bam], fasta"
    bams = [await fixture(BAM + "test.single_end.sorted.bam"), await fixture(BAM + "test.paired_end.sorted.bam")]
    result = await homer_maketagdirectory(
        bams, await fixture(GENOME + "genome.fasta"), prefix="test", args="-format sam"
    )
    await check_tagdir(
        result,
        {
            "MT192765.1.tags.tsv": "365808c4751ef6dd7085ac52037a22bc",
            "tagAutocorrelation.txt": "8b396f2aef1cdd3af4fab57b142d3250",
            "tagCountDistribution.txt": "afc6d007096c3872bbe84c9dc8edb832",
            "tagInfo.txt": "aebf6ff15fd0a238ee6a94d623c578ca",
            "tagLengthDistribution.txt": "44f231adb2a705ae81950808c55cf248",
        },
    )


@env.task
async def test_homer_findpeaks_tagdir() -> None:
    # upstream case: homer/findpeaks "sarscov2 - tagdir"
    tagdir = await bed_tagdir()
    peaks = await homer_findpeaks(tagdir.tagdir, prefix="test", args="-style factor")
    assert peaks.name == "test.peaks.txt", peaks.name
    await assert_md5(peaks, "86e15beaa4b439585786478e58418c0c", label="peaks")


async def uniqmap_dir() -> Dir:
    """upstream's UNZIP of the uniqmap archive: a directory named after the archive, holding its contents."""
    archive = await fixture(UNIQMAP_ZIP)
    root = Path(tempfile.mkdtemp()) / "uniqmap.GRCh38_chr21.50nt"
    with zipfile.ZipFile(await archive.download()) as zf:
        zf.extractall(root)
    return await Dir.from_local(str(root))


@env.task
async def test_homer_findpeaks_uniqmap() -> None:
    # upstream case: homer/findpeaks "sarscov2 - tagdir - uniqmap"
    tagdir = await bed_tagdir()
    uniqmap = await uniqmap_dir()
    assert uniqmap.name == "uniqmap.GRCh38_chr21.50nt", uniqmap.name
    peaks = await homer_findpeaks(
        tagdir.tagdir, uniqmap=uniqmap, prefix="test-uniqmap.GRCh38_chr21", args="-style factor"
    )
    assert peaks.name == "test-uniqmap.GRCh38_chr21.peaks.txt", peaks.name
    await assert_md5(peaks, "d96e7bba4f4fe827a18acae200b422a4", label="peaks")


@env.task
async def test_homer_makeucscfile_tagdir() -> None:
    # upstream case: homer/makeucscfile "sarscov2 - tagdir" (upstream snapshots only the name)
    tagdir = await bed_tagdir()
    bedgraph = await homer_makeucscfile(tagdir.tagdir, prefix="test")
    assert bedgraph.name == "test.bedGraph.gz", bedgraph.name
    async with bedgraph.open("rb") as fh:
        lines = gzip.decompress(bytes(await fh.read())).decode().splitlines()
    track = 'track type=bedGraph name="test_tagdir Total Tags = 7.00e+00, normalized to 1.00e+07"'
    assert lines[0].startswith(track), f"bedGraph track line: {lines[0]!r}"
    assert len(lines) > 1 and all(line.startswith("MT192765.1\t") for line in lines[1:]), "bedGraph body"


@env.task
async def test_homer_pos2bed_peaks() -> None:
    # upstream case: homer/pos2bed "sarscov2 - peaks"
    tagdir = await bed_tagdir()
    peaks = await homer_findpeaks(tagdir.tagdir, prefix="test")
    bed = await homer_pos2bed(peaks, prefix="test")
    assert bed.name == "test.bed", bed.name
    await assert_md5(bed, "5d6ddd9c7e621a66f6f045b9b5abecb4", label="bed")


@env.task
async def test_homer_annotatepeaks() -> None:
    # upstream case: homer/annotatepeaks "test-homer-annotatepeaks" (the table is unstable upstream)
    result = await homer_annotatepeaks(
        await fixture(GENOME + "bed/test.bed"),
        await fixture(GENOME + "genome.fasta"),
        await fixture(GENOME + "genome.gtf"),
        prefix="test",
        args="-annStats test.annStats.txt",
    )
    assert result.txt.name == "test.annotatePeaks.txt", result.txt.name
    assert result.stats is not None and result.stats.name == "test.annStats.txt"
    await assert_md5(result.stats, "eba336bce28fb5ec5d4ea3215b502dc5", label="annStats")
    lines = (await text(result.txt)).splitlines()
    cmd = "PeakID (cmd=annotatePeaks.pl test.bed genome.fasta -annStats test.annStats.txt -gtf genome.gtf -cpu "
    assert lines[0].startswith(cmd), f"annotation header: {lines[0][:200]!r}"
    peaks = sorted(line.split("\t")[0] for line in lines[1:])
    expected = ["nCoV-2019_5_LEFT", "nCoV-2019_5_RIGHT", "nCoV-2019_6_LEFT", "nCoV-2019_6_RIGHT"]
    assert peaks == expected, f"annotated peaks: {peaks}"


tests = [
    test_homer_maketagdirectory_bed,
    test_homer_maketagdirectory_bam,
    test_homer_findpeaks_tagdir,
    test_homer_findpeaks_uniqmap,
    test_homer_makeucscfile_tagdir,
    test_homer_pos2bed_peaks,
    test_homer_annotatepeaks,
]
