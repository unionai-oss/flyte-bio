"""Tests for flyte_bio.modules.ultra.

The cases and expected md5s are upstream's ultra/index and ultra/align tests
(modules commit 9a1cc9f, uLTRA 0.1, samtools 1.17), with their setup chains:
the GTF sorted by ``sort -k1,1 -k4,4n`` (prefix ``<id>_sorted`` for the index
case, ``<id>`` for align's), indexed with ``--disable_infer``, and the HiFi
reads gunzipped before alignment. Upstream hashes ``database.db`` and only
names the pickles, whose bytes aren't stable; the aligned BAM is compared by
its reads' md5.
"""

import asyncio
import os

from flyte.io import File

from flyte_bio.modules.gnusort import gnu_sort
from flyte_bio.modules.gunzip import gunzip
from flyte_bio.modules.ultra import UltraIndexResult, ultra_align, ultra_index
from tests.framework import assert_md5, assert_reads_md5, env, fixture

GENOME = "genomics/homo_sapiens/genome/"
PICKLES = [
    "all_splice_pairs_annotations",
    "all_splice_sites_annotations",
    "chr_to_id",
    "exon_choordinates_to_id",
    "flank_choordinates",
    "gene_to_small_segments",
    "id_to_chr",
    "max_intron_chr",
    "parts_to_segments",
    "ref_exon_sequences",
    "ref_flank_sequences",
    "ref_part_sequences",
    "ref_segment_sequences",
    "refs_id_lengths",
    "refs_lengths",
    "segment_id_to_choordinates",
    "segment_to_gene",
    "segment_to_ref",
    "splices_to_transcripts",
    "transcripts_to_splices",
]


async def index_genome(gtf_prefix: str) -> tuple[File, UltraIndexResult]:
    fasta, gtf = await asyncio.gather(fixture(GENOME + "genome.fasta"), fixture(GENOME + "genome.gtf"))
    sorted_gtf = await gnu_sort(gtf, "gtf", prefix=gtf_prefix, args="-k1,1 -k4,4n")
    return fasta, await ultra_index(fasta, sorted_gtf, args="--disable_infer")


@env.task
async def test_ultra_index() -> None:
    # upstream case: ultra/index "ultra_index"
    result = (await index_genome("test_sorted"))[1]
    await assert_md5(result.database, "f9b780bf5a0a5d23256afb99df7137d9", label="database.db")
    names = sorted(os.path.basename(f.path) for f in await result.index.list_files())
    expected = sorted(["database.db", *(f"{p}.pickle" for p in PICKLES)])
    assert names == expected, f"index files: expected {expected}, got {names}"


@env.task
async def test_ultra_align() -> None:
    # upstream case: ultra/align "ultra_align"
    (fasta, index), reads = await asyncio.gather(
        index_genome("test"), fixture("genomics/homo_sapiens/pacbio/fastq/test_hifi.fastq.gz")
    )
    bam = await ultra_align(await gunzip(reads), fasta, index.index, prefix="test")
    assert bam.name == "test.bam", bam.name
    await assert_reads_md5(bam, "4a12eec345cffd63b51bd3d29d89f13c", label="test.bam")


tests = [test_ultra_index, test_ultra_align]
