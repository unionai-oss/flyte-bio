"""Tests for flyte_bio.modules.starfusion.

The cases are upstream's starfusion/build and starfusion/detect tests
(modules commit 9a1cc9f, STAR-Fusion 1.15.1), with upstream's inputs: the
minigenome, its GTF, the mini fusion annotation library, Pfam and a small
Dfam, and the Broad's ``AnnotFilterRule.pm``.

- build: upstream snapshots the library's file names (its contents are
  ``unstableKeys``); this checks the same names, read from the archive.
  ctatsplicing/prepgenomelib's test checks the contents' md5s.
- detect: the md5s of the predictions and the abridged predictions, and the
  first column of the coding-effect table, with ``--examine_coding_effect``
  and upstream's ``<id>.starfusion`` prefix, for paired- and single-end reads.

The library comes back as a ``.tar.gz`` (see the module docstring), so these
tests list and read it with :func:`genome_lib_members` and
:func:`genome_lib_md5s`. Upstream's tests run with 4 CPUs, as the module does
by default.
"""

import gzip
import hashlib
import tarfile

from flyte.io import File

from flyte_bio.modules.starfusion import StarFusionResult, starfusion_build, starfusion_detect
from tests.framework import assert_md5, env, fixture

HUMAN = "genomics/homo_sapiens/genome/"
ANNOT_FILTER_RULE = "https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/AnnotFilterRule.pm"
LIB_DIR = "minigenome_fasta_genome_lib_build_dir"

CHECKPOINTS = [
    "_fusion_annot_lib.idx.ok",
    "_prot_info_db.ok",
    "annotfiltrule_cp.ok",
    "blast_pairs.idx.ok",
    "cp_gene_blast_pairs.ok",
    "cp_pfam_dat.ok",
    "cp_ref_annot_cdna.ok",
    "fusion_annot_lib.cp.ok",
    "index_pfam_hits.ok",
    "index_ref_annot_cdna.ok",
    "makeblastdb.ok",
    "mm2.splice_bed.ok",
    "mm2_genome_idx.ok",
    "ref_annot.cdsplus.dfam_masked.fa.cp.ok",
    "ref_annot.cdsplus.dfam_masked.fa.idx.ok",
    "ref_annot.gtf.gene_spans.ok",
    "ref_annot.gtf.mini.sortu.ok",
    "ref_annot.gtf.ok",
    "ref_genome.fa.ok",
    "ref_genome_fai.ok",
    "trans.blast.dat.cp.ok",
    "trans.blast.dat.index.ok",
    "validate_ctat_genome_lib.ok",
]
STAR_INDEX = [
    "Genome",
    "Log.out",
    "SA",
    "SAindex",
    "build.ok",
    "chrLength.txt",
    "chrName.txt",
    "chrNameLength.txt",
    "chrStart.txt",
    "exonGeTrInfo.tab",
    "exonInfo.tab",
    "geneInfo.tab",
    "genomeParameters.txt",
    "sjdbInfo.txt",
    "sjdbList.fromGTF.out.tab",
    "sjdbList.out.tab",
    "transcriptInfo.tab",
]
LIB_FILES = [
    "AnnotFilterRule.pm",
    "PFAM.domtblout.dat.gz",
    "blast_pairs.dat.gz",
    "blast_pairs.idx",
    "fusion_annot_lib.gz",
    "fusion_annot_lib.idx",
    "pfam_domains.dbm",
    "ref_annot.cdna.fa",
    "ref_annot.cdna.fa.idx",
    "ref_annot.cds",
    "ref_annot.cdsplus.fa",
    "ref_annot.cdsplus.fa.idx",
    "ref_annot.gtf",
    "ref_annot.gtf.gene_spans",
    "ref_annot.gtf.mini.sortu",
    "ref_annot.gtf.mm2.splice.bed",
    "ref_annot.pep",
    "ref_annot.prot_info.dbm",
    "ref_genome.fa",
    "ref_genome.fa.fai",
    "ref_genome.fa.mm2",
    "ref_genome.fa.ndb",
    "ref_genome.fa.nhr",
    "ref_genome.fa.nin",
    "ref_genome.fa.njs",
    "ref_genome.fa.not",
    "ref_genome.fa.nsq",
    "ref_genome.fa.ntf",
    "ref_genome.fa.nto",
    "trans.blast.align_coords.align_coords.dat",
    "trans.blast.align_coords.align_coords.dbm",
    "trans.blast.dat.gz",
]
# The build's files, by path within the library directory.
BUILD_FILES = sorted(
    LIB_FILES + [f"__chkpts/{name}" for name in CHECKPOINTS] + [f"ref_genome.fa.star.idx/{n}" for n in STAR_INDEX]
)


async def minigenome_lib() -> File:
    """Upstream's test library: STAR-Fusion's build over the minigenome."""
    dfam = [await fixture(HUMAN + "test_starfusion_dfam.hmm" + ext) for ext in ("", ".h3f", ".h3i", ".h3m", ".h3p")]
    return await starfusion_build(
        fasta=await fixture(HUMAN + "minigenome.fa"),
        gtf=await fixture(HUMAN + "minigenome.gtf"),
        fusion_annot_lib=await fixture(HUMAN + "CTAT_HumanFusionLib.mini.dat.gz"),
        pfam=await fixture(HUMAN + "Pfam-A.hmm.gz"),
        dfam=dfam,
        annot_filter_rule=await fixture(ANNOT_FILTER_RULE),
        prefix="minigenome_fasta",
    )


async def genome_lib_members(archive: File) -> tuple[list[str], list[str]]:
    """The archive's (directories, files), as sorted paths; the top-level directory is ``""``."""
    local = await archive.download()
    dirs, files = [], []
    with tarfile.open(local) as tar:
        for member in tar:
            parts = member.name.strip("/").split("/", 1)
            path = parts[1] if len(parts) > 1 else ""
            (dirs if member.isdir() else files).append(path)
    return sorted(dirs), sorted(files)


async def genome_lib_md5s(archive: File, paths: set[str]) -> dict[str, str]:
    """md5s of the archive's files at ``paths`` (within the library directory); ``.gz`` files decompressed."""
    local = await archive.download()
    md5s = {}
    with tarfile.open(local) as tar:
        for member in tar:
            parts = member.name.strip("/").split("/", 1)
            if len(parts) < 2 or parts[1] not in paths or not member.isfile():
                continue
            fh = tar.extractfile(member)
            assert fh is not None
            stream = gzip.GzipFile(fileobj=fh) if member.name.endswith(".gz") else fh
            h = hashlib.md5()
            while chunk := stream.read(1 << 20):
                h.update(chunk)
            md5s[parts[1]] = h.hexdigest()
    return md5s


@env.task
async def test_starfusion_build() -> None:
    # upstream case: starfusion/build "STARFUSION_BUILD - human - minigenome"
    lib = await minigenome_lib()
    assert lib.name == f"{LIB_DIR}.tar.gz", lib.name
    dirs, files = await genome_lib_members(lib)
    assert files == BUILD_FILES, (
        f"missing {sorted(set(BUILD_FILES) - set(files))}, extra {sorted(set(files) - set(BUILD_FILES))}"
    )
    assert dirs == ["", "__chkpts", "ref_genome.fa.star.idx"], dirs


async def first_column(file: File) -> list[str]:
    async with file.open("rb") as fh:
        return [line.split("\t")[0] for line in bytes(await fh.read()).decode().splitlines()]


async def check_detect(result: StarFusionResult, fusions_md5: str, abridged_md5: str, coding: list[str]) -> None:
    assert result.fusions.name == "test.starfusion.fusion_predictions.tsv", result.fusions.name
    assert result.abridged.name == "test.starfusion.abridged.tsv", result.abridged.name
    await assert_md5(result.fusions, fusions_md5, label="fusion predictions")
    await assert_md5(result.abridged, abridged_md5, label="abridged predictions")
    coding_effect = result.coding_effect
    assert coding_effect is not None and coding_effect.name == "test.starfusion.abridged.coding_effect.tsv", (
        coding_effect
    )
    names = await first_column(coding_effect)
    assert names == coding, names


PAIRED_CODING = [
    "#FusionName",
    "TATDN1--GSDMB",
    "ACACA--STAC2",
    "TATDN1--GSDMB",
    "CCDC6--RET",
    "TOB1--SYNRG",
    "TATDN1--GSDMB",
    "BCAS4--BCAS3",
    "RPS6KB1--SNF8",
    "BAG6--SLC44A4",
    "BCR--ABL1",
    "VAPB--IKZF3",
    "CYTH1--EIF3H",
    "RAB22A--MYO9B",
    "RARA--PKIA",
    "ARFGEF2--SULF2",
    "SUMF1--LRRFIP2",
    "ANKHD1--PCDH1",
    "MED1--ACSF2",
    "AHCTF1--NAAA",
    "MED1--STXBP4",
    "BCAS4--BCAS3",
    "ZMYND8--CEP250",
    "AC007562.1--AC012005.2",
    "VAPB--IKZF3",
    "SUMF1--LRRFIP2",
    "ZMYND8--CEP250",
    "RARA--PKIA",
    "CYTH1--EIF3H",
    "MED13--BCAS3",
    "VAPB--IKZF3",
]
SINGLE_CODING = [
    "#FusionName",
    "CCDC6--RET",
    "TATDN1--GSDMB",
    "TOB1--SYNRG",
    "ACACA--STAC2",
    "TATDN1--GSDMB",
    "BAG6--SLC44A4",
    "TATDN1--GSDMB",
    "RARA--PKIA",
    "AC007562.1--AC012005.2",
    "SUMF1--LRRFIP2",
    "RAB22A--MYO9B",
    "BCR--ABL1",
]


@env.task
async def test_starfusion_detect_paired_end() -> None:
    # upstream case: starfusion/detect "Paired end - Should run without failures"
    lib = await minigenome_lib()
    r1 = await fixture(HUMAN + "test_starfusion_rnaseq_1.fastq.gz")
    r2 = await fixture(HUMAN + "test_starfusion_rnaseq_2.fastq.gz")
    result = await starfusion_detect(
        lib, prefix="test.starfusion", reads_1=r1, reads_2=r2, args="--examine_coding_effect"
    )
    await check_detect(result, "82834fffed743afe07da82bd56d50c99", "d6d20fdd4b5cba21b9c0ebf8e0ea19ff", PAIRED_CODING)


@env.task
async def test_starfusion_detect_single_end() -> None:
    # upstream case: starfusion/detect "Single end - Should run without failures"
    lib = await minigenome_lib()
    r1 = await fixture(HUMAN + "test_starfusion_rnaseq_1.fastq.gz")
    result = await starfusion_detect(lib, prefix="test.starfusion", reads_1=r1, args="--examine_coding_effect")
    await check_detect(result, "62edb2285932b5afd88f6aba7d4ebdd1", "4f76c050e3b012320b5e526f7076ffb8", SINGLE_CODING)


tests = [test_starfusion_build, test_starfusion_detect_paired_end, test_starfusion_detect_single_end]
