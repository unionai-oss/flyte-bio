"""Tests for flyte_bio.modules.kallistobustools.

The cases and expected md5s are upstream's kallistobustools tests (modules
commit 9a1cc9f, kb-python 0.28.2):

- ref builds the test genome's index with the ``standard`` and ``nac``
  workflows and checks every output but the index (which upstream also
  leaves unhashed);
- count runs the subsampled 10x CellPlex gene-expression reads
  (``10XV3``, ``standard``, upstream's ``--cellranger -m 1``) against the
  ``standard`` reference and checks the matrix, the files of
  ``test.count`` upstream snapshots, and the lines of ``kb_info.json`` (the
  commands run) and ``run_info.json`` upstream snapshots.

The commands in ``kb_info.json`` record the threads and memory, so count runs
as upstream's test does: 2 CPUs and 4 GiB (``-t 2 -m 3G``), salted apart in
the cache.
"""

import hashlib

import flyte
from flyte.io import File

from flyte_bio.modules.kallistobustools import (
    kallistobustools_count_cmd,
    kallistobustools_ref,
    read_count_results,
)
from tests.framework import assert_md5, env, fixture

GENOME = "genomics/homo_sapiens/genome/"
GEX = (
    "genomics/homo_sapiens/10xgenomics/cellranger/10k_pbmc_cmo/fastqs/gex_1/"
    "subsampled_SC3_v3_NextGem_DI_CellPlex_Human_PBMC_10K_1_gex_S2_L001_"
)


async def kb_reference(workflow_mode: str):
    return await kallistobustools_ref(
        await fixture(GENOME + "genome.fasta"), await fixture(GENOME + "genome.gtf"), workflow_mode
    )


@env.task
async def test_kallistobustools_ref_standard() -> None:
    reference = await kb_reference("standard")
    await assert_md5(reference.cdna, "7bca59288fab822451de250d2eee48dc", label="cdna.fa")
    await assert_md5(reference.t2g, "d1a8a22c59b9cb0bda39c0c9bb3f6afe", label="t2g.txt")
    assert reference.index.name == "kb_ref_out.idx", reference.index.name
    assert reference.intron is None and reference.cdna_t2c is None and reference.intron_t2c is None


@env.task
async def test_kallistobustools_ref_nac() -> None:
    reference = await kb_reference("nac")
    await assert_md5(reference.cdna, "7bca59288fab822451de250d2eee48dc", label="cdna.fa")
    await assert_md5(reference.t2g, "58591306b33bb948bac7b40f346d0cd7", label="t2g.txt")
    assert reference.intron is not None and reference.cdna_t2c is not None and reference.intron_t2c is not None
    await assert_md5(reference.intron, "1aad4e3f5d006f495cc6647fa0bbf6ff", label="intron.fa")
    await assert_md5(reference.cdna_t2c, "23861cf43033e7c596e6989a88a3a373", label="cdna_t2c.txt")
    await assert_md5(reference.intron_t2c, "fe6d5501923867b514a0447aa4b4995f", label="intron_t2c.txt")
    assert reference.index.name == "kb_ref_out.idx", reference.index.name


async def kb_lines(file: File) -> list[str]:
    async with file.open("rb") as fh:
        return bytes(await fh.read()).decode().splitlines()


async def kb_md5(file: File) -> str:
    h = hashlib.md5()
    async with file.open("rb") as fh:
        while chunk := bytes(await fh.read(1 << 20)):
            h.update(chunk)
    return h.hexdigest()


COUNT_MD5S = {
    "10x_version3_whitelist.txt": "3d36d0a4021fd292b265e2b5e72aaaf3",
    "counts_unfiltered/cells_x_genes.barcodes.txt": "41f7adaf43b60f2f4f62d6a7073688de",
    "counts_unfiltered/cells_x_genes.genes.names.txt": "b29afa75be300c7f24fbd0740a66689b",
    "counts_unfiltered/cells_x_genes.genes.txt": "fe6d5501923867b514a0447aa4b4995f",
    "counts_unfiltered/cells_x_genes.mtx": "e12a45e7f7d6527f698dd9cb2e99ecd1",
    "inspect.json": "bafb47a58ac1bbf9be953f21c361d266",
    "matrix.ec": "31a4c1a3e8e0c562b12f6569ffbf5459",
    "output.bus": "d6fa0612a4a16eaf8a3e08bdc13ff49c",
    "output.unfiltered.bus": "bf899b967657f612ba864188868d58cc",
    "transcripts.txt": "23861cf43033e7c596e6989a88a3a373",
}
READS = (
    "subsampled_SC3_v3_NextGem_DI_CellPlex_Human_PBMC_10K_1_gex_S2_L001_R1_001.fastq.gz "
    "subsampled_SC3_v3_NextGem_DI_CellPlex_Human_PBMC_10K_1_gex_S2_L001_R2_001.fastq.gz"
)
KB_INFO_LINES = [
    '    "commands": [',
    f'        "kallisto bus -i kb_ref_out.idx -o test.count -x 10XV3 -t 2 {READS}",',
    '        "bustools sort -o test.count/tmp/output.s.bus -T test.count/tmp -t 2 -m 3G test.count/output.bus",',
    (
        '        "bustools inspect -o test.count/inspect.json -w test.count/10x_version3_whitelist.txt '
        'test.count/tmp/output.s.bus",'
    ),
    (
        '        "bustools correct -o test.count/tmp/output.s.c.bus -w test.count/10x_version3_whitelist.txt '
        'test.count/tmp/output.s.bus",'
    ),
    (
        '        "bustools sort -o test.count/output.unfiltered.bus -T test.count/tmp -t 2 -m 3G '
        'test.count/tmp/output.s.c.bus",'
    ),
    (
        '        "bustools count -o test.count/counts_unfiltered/cells_x_genes -g t2g.txt -e test.count/matrix.ec '
        '-t test.count/transcripts.txt --genecounts --umi-gene test.count/output.unfiltered.bus"'
    ),
    "    ],",
]
RUN_INFO_LINES = [
    "{",
    '\t"n_targets": 12,',
    '\t"n_bootstraps": 0,',
    '\t"n_processed": 10000,',
    '\t"n_pseudoaligned": 26,',
    '\t"n_unique": 26,',
    '\t"p_pseudoaligned": 0.3,',
    '\t"p_unique": 0.3,',
    '\t"kallisto_version": "0.50.1",',
    '\t"index_version": 13,',
]


@env.task
async def test_kallistobustools_count() -> None:
    reference = await kb_reference("standard")
    reads = [await fixture(GEX + "R1_001.fastq.gz"), await fixture(GEX + "R2_001.fastq.gz")]
    # The cache key ignores resources, so salt it to keep the 2-CPU run apart from a default one.
    task = kallistobustools_count_cmd.override(
        resources=flyte.Resources(cpu=2, memory="4Gi"), cache=flyte.Cache(behavior="auto", salt="cpu=2,mem=4Gi")
    )
    results = await task(
        reads=reads,
        index=reference.index,
        t2g=reference.t2g,
        t1c=None,
        t2c=None,
        technology="10XV3",
        workflow_mode="standard",
        prefix="test",
        args="--cellranger -m 1",
    )
    result = await read_count_results(results, "test")
    names = [f.name for f in result.matrices]
    assert names == ["cells_x_genes.mtx"], names
    await assert_md5(result.matrices[0], "e12a45e7f7d6527f698dd9cb2e99ecd1", label="cells_x_genes.mtx")
    root = result.count.path.rstrip("/") + "/"
    files = {f.path[len(root) :]: f async for f in result.count.walk()}
    mismatches = []
    for path, md5 in COUNT_MD5S.items():
        if path not in files:
            mismatches.append(f"{path}: missing (have {sorted(files)})")
        elif (actual := await kb_md5(files[path])) != md5:
            mismatches.append(f"{path}: expected {md5}, got {actual}")
    assert not mismatches, "; ".join(mismatches)
    kb_info = (await kb_lines(files["kb_info.json"]))[15:23]
    assert kb_info == KB_INFO_LINES, f"kb_info.json lines 16-23: {kb_info}"
    run_info = (await kb_lines(files["run_info.json"]))[0:10]
    assert run_info == RUN_INFO_LINES, f"run_info.json lines 1-10: {run_info}"


tests = [test_kallistobustools_ref_standard, test_kallistobustools_ref_nac, test_kallistobustools_count]
