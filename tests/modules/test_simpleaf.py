"""Tests for flyte_bio.modules.simpleaf.

Upstream's module tests (and snapshots) only check that the expected index
and quantification files exist, since piscem and alevin-fry outputs carry
timestamps and run paths. These do the same, on the same data: the
homo_sapiens chr22 genome + GTF and a 10x v3 sample, filtered with ``knee``.
"""

from flyte.io import Dir

from flyte_bio.modules.simpleaf import simpleaf_index, simpleaf_quant
from tests.framework import env, fixture

FASTQS = "genomics/homo_sapiens/10xgenomics/cellranger/5k_cmvpos_tcells/fastqs/gex_1/"
READS = "subsampled_5k_human_antiCMV_T_TBNK_connect_GEX_1_S1_L001_R{}_001.fastq.gz"


async def names(d: Dir) -> set[str]:
    """Paths of every file under ``d``, relative to it."""
    root = d.path.rstrip("/") + "/"
    return {f.path.removeprefix(root) async for f in d.walk()}


async def build_index():
    return await simpleaf_index(
        prefix="simpleaf_index",
        genome_fasta=await fixture("genomics/homo_sapiens/genome/genome.fasta"),
        gtf=await fixture("genomics/homo_sapiens/genome/genome.gtf"),
    )


@env.task
async def test_simpleaf_index() -> None:
    index = await build_index()
    files = await names(index.index)
    piscem = ("_cfish.json", ".ctab", ".ectab", ".json", ".refinfo", ".sshash")
    for name in ["simpleaf_index.json", *(f"piscem_idx{ext}" for ext in piscem)]:
        assert name in files, sorted(files)
    ref = await names(index.ref)
    for name in ("roers_ref.fa", "t2g_3col.tsv", "gene_id_to_name.tsv", "roers_make-ref.json"):
        assert name in ref, sorted(ref)
    assert index.t2g is not None and index.t2g.path.endswith("/t2g_3col.tsv"), index.t2g


@env.task
async def test_simpleaf_quant() -> None:
    index = await build_index()
    quant = await simpleaf_quant(
        reads_1=[await fixture(FASTQS + READS.format(1))],
        reads_2=[await fixture(FASTQS + READS.format(2))],
        index=index.index,
        t2g=index.t2g,
        chemistry="10xv3",
        prefix="test",
        cell_filter="knee",
    )
    mapped = await names(quant.map)
    assert {"map.rad", "map_info.json", "unmapped_bc_count.bin"} <= mapped, sorted(mapped)
    quantified = await names(quant.quant)
    expected = {
        "permit_map.bin",
        "collate.json",
        "generate_permit_list.json",
        "quant.json",
        "featureDump.txt",
        "permit_freq.bin",
        "unmapped_bc_count_collated.bin",
        "alevin/quants_mat.mtx",
        "alevin/quants_mat_cols.txt",
        "alevin/quants_mat_rows.txt",
        "alevin/quants.h5ad",
    }
    assert expected <= quantified, sorted(quantified)


tests = [test_simpleaf_index, test_simpleaf_quant]
