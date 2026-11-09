"""Tests for flyte_bio.modules.ctatsplicing.

The cases are upstream's ctatsplicing/prepgenomelib and
ctatsplicing/startocancerintrons tests (modules commit 9a1cc9f, CTAT-Splicing
0.0.3), on upstream's library: STAR-Fusion's build over the minigenome (see
test_starfusion), with the Broad's GRCh38 cancer introns added.

- prepgenomelib: upstream snapshots the names of every file and directory in
  the library and the md5s of its non-empty files but those it lists as
  unstable (``.gz`` files decompressed, as nf-test hashes them). These tests
  check the same, reading the library archive.
- startocancerintrons: upstream aligns ``test_rnaseq_{1,2}.fastq.gz`` as one
  single-end sample (two comma-separated files) with STAR against the
  library's own index and GTF, indexes the BAM and runs with
  ``--min_total_reads 1 --vis``; it snapshots the md5s of
  ``test.cancer.introns`` and ``test.cancer.introns.prelim`` and the name of
  ``test.introns``. STAR reads comma-separated files one after the other, so
  here the two files are concatenated into one (cat_fastq) instead. The STAR
  index and GTF are unpacked from the library archive in the test.
"""

import fnmatch
import tarfile
import tempfile
from pathlib import Path

from flyte.io import Dir, File

from flyte_bio.modules.cat import cat_fastq
from flyte_bio.modules.ctatsplicing import ctatsplicing_prepgenomelib, ctatsplicing_startocancerintrons
from flyte_bio.modules.samtools import samtools_index
from flyte_bio.modules.star import star_align
from tests.framework import assert_md5, env, fixture
from tests.modules.test_starfusion import genome_lib_md5s, minigenome_lib

CANCER_INTRONS = (
    "https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/CANCER_SPLICING_LIB_SUPPLEMENT/"
    "cancer_introns.GRCh38.Jun232020.tsv.gz"
)

# Upstream's snapshot: the name of every file and directory of the prepared library, sorted.
PREPARED_NAMES = sorted(
    [
        "minigenome_fasta_genome_lib_build_dir",
        "AnnotFilterRule.pm",
        "PFAM.domtblout.dat.gz",
        "__chkpts",
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
        "blast_pairs.dat.gz",
        "blast_pairs.idx",
        "cancer_splicing_lib",
        "cancer_splicing.idx",
        "fusion_annot_lib.gz",
        "fusion_annot_lib.idx",
        "pfam_domains.dbm",
        "refGene.bed",
        "refGene.sort.bed.gz",
        "refGene.sort.bed.gz.tbi",
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
        "ref_genome.fa.star.idx",
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
        "trans.blast.align_coords.align_coords.dat",
        "trans.blast.align_coords.align_coords.dbm",
        "trans.blast.dat.gz",
    ]
)
# Upstream's snapshot md5s, by path within the library directory.
STAR_IDX = "ref_genome.fa.star.idx/"
PREPARED_MD5S = {
    "AnnotFilterRule.pm": "f94966013cd0df9624a6dda0b75fefa0",
    "fusion_annot_lib.gz": "23d82a5da81f91ca4e1ecd6481992a12",
    "refGene.sort.bed.gz": "01bf4093e3e9dcad8ed3a0f45a381c95",
    "ref_annot.gtf": "5ce8afe99ef3940a877a04caeacf9181",
    "ref_annot.gtf.mini.sortu": "1d29ccecdbb7b40a99c84a02d6c2c1be",
    "ref_annot.gtf.mm2.splice.bed": "340585ea1843bf06bf555575ddecf28c",
    "ref_genome.fa": "ad699c56ed38566c7d3e9579486b1706",
    "ref_genome.fa.fai": "e3f74a27219b33ae80fd5de5cbeaf32b",
    "ref_genome.fa.mm2": "ce50979ea284748eb9f84ae88cfd930e",
    "ref_genome.fa.ndb": "6ea574753b557610f62f6e4ab79e19f5",
    "ref_genome.fa.nhr": "50f28dae71683c4394bfaf94a1ef4392",
    "ref_genome.fa.not": "1e53e9d08f1d23af0299cfa87478a7bb",
    "ref_genome.fa.nsq": "d2361e7871ce4cf51181c112a48f191b",
    "ref_genome.fa.ntf": "de1250813f0c7affc6d12dac9d0fb6bb",
    "ref_genome.fa.nto": "33cdeccccebe80329f1fdbee7f5874cb",
    STAR_IDX + "Genome": "9e3efdd0901cabb5a2d589664a63b372",
    STAR_IDX + "SA": "7dd9083264be9c6a2194d990bc10d237",
    STAR_IDX + "SAindex": "ac4711df685109e04356db9e9cb9fb7f",
    STAR_IDX + "chrLength.txt": "e02cd536b7281b894246863b160d5d06",
    STAR_IDX + "chrName.txt": "f4d0d6595f423084e6b9472e40dfe6e8",
    STAR_IDX + "chrNameLength.txt": "07a67d7ac441d7d30d80840b0927e717",
    STAR_IDX + "chrStart.txt": "e2031239a74fe5ee9051e9364e4f608a",
    STAR_IDX + "exonGeTrInfo.tab": "3c35618d07a8e35a0f9108699fcdda42",
    STAR_IDX + "exonInfo.tab": "bcbb3f32fa31fe504cc737f337ad341c",
    STAR_IDX + "geneInfo.tab": "db5db4b6e003904e9908fce7c05f0125",
    STAR_IDX + "sjdbInfo.txt": "e4cc1bbf8bd687cfc3d7c2c702e6def7",
    STAR_IDX + "sjdbList.fromGTF.out.tab": "8f3e8604b00d4067e4eb80aa476a8113",
    STAR_IDX + "sjdbList.out.tab": "5d78dd49d5db24ca2c056b7ebe5c2059",
    STAR_IDX + "transcriptInfo.tab": "b758c0ccaddcf0453bab5905b3cec4a2",
    "trans.blast.align_coords.align_coords.dat": "9f6b7a75aea03a9671190be25ecdd4c2",
    "trans.blast.dat.gz": "85ba5ea96c566f751ad83a3e4b8ab128",
}
# The files upstream leaves out of its md5s, as its (top-level or any-depth) glob patterns.
UNSTABLE = [
    "refGene.sort.bed.gz.tbi",
    "refGene.bed",
    "*/cancer_splicing.idx",
    "ref_annot.c*",
    "PFAM.domtblout.dat.gz",
    "blast_pairs*",
    "fusion_annot_lib.idx",
    "pfam_domains.dbm",
    "ref_annot.gtf.gene_spans",
    "ref_annot.prot_info.dbm",
    "ref_annot.pep",
    "ref_genome.fa.nin",
    "ref_genome.fa.njs",
    "*/Log.out",
    "*/genomeParameters.txt",
    "trans.blast.align_coords.align_coords.dbm",
]


async def prepared_minigenome_lib() -> File:
    """Upstream's prepared test library: the minigenome library with the GRCh38 cancer introns."""
    return await ctatsplicing_prepgenomelib(await minigenome_lib(), await fixture(CANCER_INTRONS))


@env.task
async def test_ctatsplicing_prepgenomelib() -> None:
    # upstream case: ctatsplicing/prepgenomelib "homo_sapiens - reference"
    lib = await prepared_minigenome_lib()
    assert lib.name == "minigenome_fasta_genome_lib_build_dir.tar.gz", lib.name
    local = await lib.download()
    names, stable = [], set()
    with tarfile.open(local) as tar:
        for member in tar:
            parts = member.name.strip("/").split("/", 1)
            names.append(Path(member.name).name)
            if len(parts) < 2 or not member.isfile() or member.size == 0:
                continue
            if not any(fnmatch.fnmatchcase(parts[1], pattern) for pattern in UNSTABLE):
                stable.add(parts[1])
    assert sorted(names) == PREPARED_NAMES, (
        f"missing {sorted(set(PREPARED_NAMES) - set(names))}, extra {sorted(set(names) - set(PREPARED_NAMES))}"
    )
    assert stable == set(PREPARED_MD5S), f"stable files {sorted(stable ^ set(PREPARED_MD5S))} differ"
    md5s = await genome_lib_md5s(lib, set(PREPARED_MD5S))
    wrong = {path: md5 for path, md5 in md5s.items() if md5 != PREPARED_MD5S[path]}
    assert not wrong, f"md5 mismatches: {wrong}"


async def star_index_and_gtf(lib: File) -> tuple[Dir, File]:
    """The library's STAR index (``ref_genome.fa.star.idx``) and GTF (``ref_annot.gtf``), unpacked."""
    local = await lib.download()
    root = Path(tempfile.mkdtemp())
    with tarfile.open(local) as tar:
        top = tar.next()
        assert top is not None
        members = [m for m in tar if "/ref_genome.fa.star.idx/" in m.name or m.name.endswith("/ref_annot.gtf")]
        tar.extractall(root, members=members, filter="data")
    lib_dir = root / top.name.strip("/").split("/")[0]
    return await Dir.from_local(str(lib_dir / "ref_genome.fa.star.idx")), await File.from_local(
        str(lib_dir / "ref_annot.gtf")
    )


# Upstream's test config for STAR_ALIGN, plus the read group upstream's module adds.
STAR_ARGS = (
    "--readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outReadsUnmapped None --twopassMode None "
    "--outSAMstrandField intronMotif --outSAMunmapped Within --chimSegmentMin 12 --chimJunctionOverhangMin 8 "
    "--chimOutJunctionFormat 1 --alignSJDBoverhangMin 10 --alignMatesGapMax 100000 --alignIntronMax 100000 "
    "--alignSJstitchMismatchNmax 5 -1 5 5 --chimMultimapScoreRange 3 --chimScoreJunctionNonGTAG -4 "
    "--chimMultimapNmax 20 --chimNonchimScoreDropMin 10 --peOverlapNbasesMin 12 --peOverlapMMp 0.1 "
    "--alignInsertionFlush Right --alignSplicedMateMapLminOverLmate 0 --alignSplicedMateMapLmin 30 "
    "--outSAMattrRGline ID:test SM:test"
)


async def star_output(results: Dir, name: str) -> File:
    found = await results.get_file(name)
    assert found is not None, f"STAR wrote no {name}"
    return found


@env.task
async def test_ctatsplicing_startocancerintrons() -> None:
    # upstream case: ctatsplicing/startocancerintrons "homo_sapiens - bam - junctions"
    lib = await prepared_minigenome_lib()
    index, gtf = await star_index_and_gtf(lib)
    reads = await cat_fastq(
        [
            await fixture("genomics/homo_sapiens/illumina/fastq/test_rnaseq_1.fastq.gz"),
            await fixture("genomics/homo_sapiens/illumina/fastq/test_rnaseq_2.fastq.gz"),
        ],
        prefix="test",
    )
    aligned = await star_align(reads_1=reads, reads_2=None, index=index, gtf=gtf, args=STAR_ARGS)
    bam = await star_output(aligned, "Aligned.sortedByCoord.out.bam")
    bai = await samtools_index(bam=bam)
    result = await ctatsplicing_startocancerintrons(
        split_junction=await star_output(aligned, "SJ.out.tab"),
        junction=await star_output(aligned, "Chimeric.out.junction"),
        genome_lib=lib,
        prefix="test",
        bam=bam,
        bai=bai,
        args="--min_total_reads 1 --vis",
    )
    await assert_md5(result.cancer_introns, "c0eb94972d845218c2e470b85c55871b", label="test.cancer.introns")
    await assert_md5(
        result.cancer_introns_prelim, "c0eb94972d845218c2e470b85c55871b", label="test.cancer.introns.prelim"
    )
    assert result.introns.name == "test.introns", result.introns.name


tests = [test_ctatsplicing_prepgenomelib, test_ctatsplicing_startocancerintrons]
