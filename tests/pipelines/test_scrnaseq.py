"""End-to-end tests for flyte_bio.pipelines.scrnaseq on the upstream test profile.

The test profile: two 10x v2 samples (Sample_Y with two lanes) against
mouse chr19 (GRCm38.p6 genome, GENCODE vM19 annotation).
"""

import asyncio
import gzip
import json
import math

import anndata
import flyte
from flyte.io import Dir, File

from flyte_bio.modules.simpleaf import SimpleafQuant
from flyte_bio.pipelines.scrnaseq import (
    convert_matrices,
    get_barcode_whitelist,
    get_protocol,
    multiqc_report,
    prepare_index,
    prepare_reference,
    quantify,
    read_qc,
    remove_background,
    samples_from_samplesheet,
    scrnaseq,
)
from flyte_bio.samplesheet import read_samplesheet
from tests.framework import assert_nonempty, env, fixture, stage_samplesheet

DATA = "https://github.com/nf-core/test-datasets/raw/scrnaseq/"
SAMPLESHEET = DATA + "samplesheet-2-0.csv"
FASTA = DATA + "reference/GRCm38.p6.genome.chr19.fa"
GTF = DATA + "reference/gencode.vM19.annotation.chr19.gtf"


async def text(f: File) -> str:
    async with f.open("rb") as fh:
        return bytes(await fh.read()).decode()


async def staged_samples():
    sheet = await stage_samplesheet(SAMPLESHEET, ("fastq_1", "fastq_2"))
    return sheet, samples_from_samplesheet(await read_samplesheet(sheet))


@env.task
async def test_prepare_reference() -> None:
    fasta, gtf = await asyncio.gather(fixture(FASTA), fixture(GTF))
    ref = await prepare_reference(fasta, gtf)
    assert ref.fasta is fasta
    assert ref.gtf.path.endswith("/GRCm38.p6.genome.chr19_genes.gtf"), ref.gtf.path
    # Every annotation line is on chr19, the genome's only sequence, so all are kept.
    original = [line for line in (await text(gtf)).splitlines() if line]
    kept = [line for line in (await text(ref.gtf)).splitlines() if line]
    assert kept == [line for line in original if line.split("\t")[0] == "chr19"], (len(kept), len(original))


@env.task
async def test_read_qc() -> None:
    _, samples = await staged_samples()
    assert [(s.id, len(s.fastq_1), s.expected_cells, s.seq_center) for s in samples] == [
        ("Sample_X", 1, 5000, "Broad Institute"),
        ("Sample_Y", 2, 5000, "CRG Barcelona"),
    ], samples
    reports = await read_qc(samples[1])
    # Two lanes of R1 + R2: four reports, numbered in lane order.
    for n in (1, 2, 3, 4):
        assert await reports.get_file(f"Sample_Y_{n}_fastqc.html") is not None, f"no Sample_Y_{n}_fastqc.html"


async def read_count(f: File) -> int:
    local = await f.download()
    with gzip.open(local, "rb") as fh:
        return sum(1 for _ in fh) // 4


@env.task
async def test_quantify() -> None:
    _, samples = await staged_samples()
    sample = samples[1]  # Sample_Y: two lanes
    fasta, gtf = await asyncio.gather(fixture(FASTA), fixture(GTF))
    proto = get_protocol("10XV2")
    reference, whitelist = await asyncio.gather(prepare_reference(fasta, gtf), get_barcode_whitelist(proto))
    assert whitelist is not None and whitelist.path.endswith("10x_V2_barcode_whitelist.txt.gz"), whitelist
    index = await prepare_index(reference)
    assert index.t2g is not None and index.t2g.path.endswith("/t2g_3col.tsv"), index.t2g
    quant = await quantify(sample, index, proto, whitelist)
    for name in ("quants_mat.mtx", "quants_mat_cols.txt", "quants_mat_rows.txt", "quants.h5ad"):
        assert await quant.quant.get_file(f"alevin/{name}") is not None, f"no alevin/{name}"
    # Both lanes were mapped: the mapper saw every read pair of the sample.
    info = json.loads(await text(await quant.map.get_file("map_info.json")))
    reads = sum(await asyncio.gather(*(read_count(f) for f in sample.fastq_1)))
    assert reads in info.values(), (reads, info)


async def read_h5ad(f: File):
    return anndata.read_h5ad(await f.download())


async def simpleaf_counts(quant: SimpleafQuant) -> tuple[int, int, float]:
    """(cells, genes, total count) of simpleaf's own h5ad export in ``alevin/``.

    Not ``quants_mat.mtx``: that keeps spliced, unspliced and ambiguous counts
    as separate columns per gene, while the h5ad (the conversion's input) has
    one column per gene.
    """
    m = await read_h5ad(await quant.quant.get_file("alevin/quants.h5ad"))
    return m.n_obs, m.n_vars, float(m.X.sum())


@env.task
async def test_convert_matrices() -> None:
    # Upstream's snapshot pins these files' md5s, but they record each tool's command
    # line (with its file paths) in uns, so they can't match byte for byte here: the
    # checks are on content instead.
    sheet, samples = await staged_samples()
    fasta, gtf = await asyncio.gather(fixture(FASTA), fixture(GTF))
    proto = get_protocol("10XV2")
    reference, whitelist = await asyncio.gather(prepare_reference(fasta, gtf), get_barcode_whitelist(proto))
    index = await prepare_index(reference)
    quants = await asyncio.gather(*(quantify(s, index, proto, whitelist) for s in samples))
    matrices = await convert_matrices({s.id: q for s, q in zip(samples, quants)}, sheet)
    assert list(matrices.samples) == ["Sample_X", "Sample_Y"], list(matrices.samples)

    total_cells = 0
    for (sample_id, matrix), quant in zip(matrices.samples.items(), quants):
        assert matrix.h5ad.path.endswith(f"/{sample_id}_raw_matrix.h5ad"), matrix.h5ad.path
        adata = await read_h5ad(matrix.h5ad)
        cells, genes, total = await simpleaf_counts(quant)
        assert adata.shape == (cells, genes), (sample_id, adata.shape, cells, genes)
        assert math.isclose(float(adata.X.sum()), total, rel_tol=1e-6), (sample_id, float(adata.X.sum()), total)
        assert list(adata.obs_names) == sorted(adata.obs_names), sample_id
        assert set(adata.obs["sample"]) == {sample_id}, set(adata.obs["sample"])
        assert not any("." in g for g in adata.var_names), sample_id
        assert "runtime_seconds" not in json.loads(adata.uns["simpleaf_map_info"]), sample_id
        total_cells += cells
        for f in (matrix.rds.sce, matrix.rds.seurat):
            await assert_nonempty(f, label=f.path)

    combined = matrices.combined
    assert combined.h5ad.path.endswith("/combined_raw_matrix.h5ad"), combined.h5ad.path
    adata = await read_h5ad(combined.h5ad)
    # Upstream labels each cell with its file's name minus "_matrix.h5ad".
    assert adata.n_obs == total_cells, (adata.n_obs, total_cells)
    assert set(adata.obs["sample"]) == {"Sample_X_raw", "Sample_Y_raw"}, set(adata.obs["sample"])
    for f in (combined.rds.sce, combined.rds.seurat):
        await assert_nonempty(f, label=f.path)


@env.task
async def test_multiqc_report() -> None:
    _, samples = await staged_samples()
    fasta, gtf = await asyncio.gather(fixture(FASTA), fixture(GTF))
    proto = get_protocol("10XV2")
    reference, whitelist = await asyncio.gather(prepare_reference(fasta, gtf), get_barcode_whitelist(proto))
    index = await prepare_index(reference)
    fastqc, quants = await asyncio.gather(
        asyncio.gather(*(read_qc(s) for s in samples)),
        asyncio.gather(*(quantify(s, index, proto, whitelist) for s in samples)),
    )
    ids = [s.id for s in samples]
    report = await multiqc_report(dict(zip(ids, fastqc)), dict(zip(ids, quants)))
    await assert_nonempty(report.report, label="multiqc_report.html")

    files = [f async for f in report.results.walk()]
    table = next((f for f in files if f.path.endswith("_data/multiqc_fastqc.txt")), None)
    assert table is not None, sorted(f.path.rsplit("/", 1)[-1] for f in files)
    reported = {line.split("\t")[0] for line in (await text(table)).splitlines()[1:]}
    # Every lane's R1 and R2: Sample_X has one lane, Sample_Y two.
    expected = {"Sample_X_1", "Sample_X_2", *(f"Sample_Y_{n}" for n in (1, 2, 3, 4))}
    assert reported == expected, reported
    # Upstream's config exports the plots too.
    assert any("_plots/" in f.path for f in files), "no exported plots"


@env.task
async def test_remove_background() -> None:
    # Upstream's pipeline test skips CellBender, which can't train on the toy reads' few
    # cells, so this runs it on the module test's raw PBMC matrix for 5 epochs instead.
    sheet, _ = await staged_samples()
    raw = await fixture("genomics/homo_sapiens/scrnaseq/h5ad/SRR28679759_raw_matrix.h5ad")
    result = await remove_background({"pbmc": raw}, sheet, cellbender_args="--epochs 5")
    assert list(result.cellbender) == ["pbmc"], list(result.cellbender)

    matrix = result.matrices.samples["pbmc"]
    assert matrix.h5ad.path.endswith("/pbmc_cellbender_filter_matrix.h5ad"), matrix.h5ad.path
    cells = (await text(result.cellbender["pbmc"].barcodes)).split()
    adata, unfiltered = await asyncio.gather(read_h5ad(matrix.h5ad), read_h5ad(raw))
    assert 0 < len(cells) < unfiltered.n_obs, (len(cells), unfiltered.n_obs)
    # The raw counts of CellBender's cells, in its order.
    assert list(adata.obs_names) == cells, (adata.obs_names[:3], cells[:3])
    assert list(adata.var_names) == list(unfiltered.var_names)
    assert math.isclose(float(adata.X.sum()), float(unfiltered[cells].X.sum()), rel_tol=1e-6)
    for f in (matrix.rds.sce, matrix.rds.seurat):
        await assert_nonempty(f, label=f.path)

    combined = result.matrices.combined
    assert combined.h5ad.path.endswith("/combined_cellbender_filter_matrix.h5ad"), combined.h5ad.path
    adata = await read_h5ad(combined.h5ad)
    assert adata.n_obs == len(cells), (adata.n_obs, len(cells))
    assert set(adata.obs["sample"]) == {"pbmc_cellbender_filter"}, set(adata.obs["sample"])
    for f in (combined.rds.sce, combined.rds.seurat):
        await assert_nonempty(f, label=f.path)


# Upstream's simpleaf test-profile snapshot (tests/main_pipeline_simpleaf.nf.test.snap): every published
# file outside multiqc/ and pipeline_info/.
AF_MAP = ["map.rad", "map_info.json", "unmapped_bc_count.bin"]
AF_QUANT = [
    "alevin/quants.h5ad",
    "alevin/quants_mat.mtx",
    "alevin/quants_mat_cols.txt",
    "alevin/quants_mat_rows.txt",
    "collate.json",
    "featureDump.txt",
    "gene_id_to_name.tsv",
    "generate_permit_list.json",
    "map.collated.rad",
    "permit_freq.bin",
    "permit_map.bin",
    "quant.json",
    "simpleaf_map_info.json",
    "unmapped_bc_count_collated.bin",
]
LANE_READS = {"Sample_X": 2, "Sample_Y": 4}  # FastQC reports: R1 and R2 per lane


def matrices(prefix: str, input_type: str) -> list[str]:
    return [f"{prefix}_{input_type}_matrix.{ext}" for ext in ("h5ad", "sce.rds", "seurat.rds")]


def upstream_paths(fastqc: bool = True, af_map: bool = True) -> set[str]:
    paths = {f"simpleaf/mtx_conversions/{f}" for f in matrices("combined", "raw")}
    for sid, reads in LANE_READS.items():
        if fastqc:
            paths |= {f"fastqc/{sid}_{n}_fastqc.{ext}" for n in range(1, reads + 1) for ext in ("html", "zip")}
        quant = f"simpleaf/{sid}/simpleaf_quant"
        if af_map:
            paths |= {f"{quant}/af_map/{f}" for f in AF_MAP}
        paths |= {f"{quant}/af_quant/{f}" for f in AF_QUANT}
        paths |= {f"simpleaf/mtx_conversions/{sid}/{f}" for f in matrices(sid, "raw")}
    return paths


async def published(tree: Dir) -> set[str]:
    base = tree.path.rstrip("/")
    return {f.path[len(base) :].lstrip("/") async for f in tree.walk()}


async def profile_inputs() -> dict:
    """scrnaseq() keyword arguments for upstream's test profile, which skips QCatch and CellBender."""
    sheet, fasta, gtf = await asyncio.gather(
        stage_samplesheet(SAMPLESHEET, ("fastq_1", "fastq_2")), fixture(FASTA), fixture(GTF)
    )
    return {
        "samplesheet": sheet,
        "fasta": fasta,
        "gtf": gtf,
        "protocol": "10XV2",
        "skip_qcatch": True,
        "skip_cellbender": True,
    }


@env.task
async def test_scrnaseq() -> None:
    # Upstream's simpleaf pipeline test, published in upstream's --outdir layout.
    result = await scrnaseq(**await profile_inputs(), publish_results=True)
    assert [r.sample.id for r in result.samples] == ["Sample_X", "Sample_Y"]
    assert result.cellbender is None and result.multiqc is not None and result.outdir is not None
    paths = await published(result.outdir)
    others = {p for p in paths if not p.startswith("multiqc/")}
    assert others == upstream_paths(), (sorted(others - upstream_paths()), sorted(upstream_paths() - others))
    for path in ("multiqc/multiqc_report.html", "multiqc/multiqc_data/multiqc_fastqc.txt"):
        assert path in paths, f"not published: {path}"
    assert any(p.startswith("multiqc/multiqc_plots/") for p in paths), "no MultiQC plots published"


@env.task
async def test_scrnaseq_options() -> None:
    # With FastQC and MultiQC skipped, the reference saved and the mapping left out, published to an
    # explicit outdir.
    outdir = flyte.ctx().raw_data_path.get_random_remote_path("scrnaseq_outdir")
    result = await scrnaseq(
        **await profile_inputs(),
        skip_fastqc=True,
        skip_multiqc=True,
        save_reference=True,
        save_align_intermeds=False,
        outdir=outdir,
    )
    assert result.multiqc is None and all(r.fastqc is None for r in result.samples)
    assert result.outdir is not None and result.outdir.path.rstrip("/") == outdir.rstrip("/")
    paths = await published(result.outdir)
    reference = {p for p in paths if p.startswith(("gtf_filter/", "simpleaf/simpleaf_index/"))}
    assert paths - reference == upstream_paths(fastqc=False, af_map=False), sorted(paths - reference)
    assert "gtf_filter/GRCm38.p6.genome.chr19_genes.gtf" in reference, sorted(reference)
    assert any(p.startswith("simpleaf/simpleaf_index/index/") for p in reference), sorted(reference)
    assert any(p.startswith("simpleaf/simpleaf_index/ref/") for p in reference), sorted(reference)


@env.task
async def test_protocols() -> None:
    assert get_protocol("10XV3").chemistry == "10xv3"
    # Anything else goes to simpleaf as it is, as upstream: a chemistry name or custom geometry.
    custom = get_protocol("1{b[16]u[12]x:}2{r:}")
    assert (custom.chemistry, custom.whitelist, custom.qcatch) == ("1{b[16]u[12]x:}2{r:}", None, None)
    try:
        get_protocol("auto")
    except ValueError:
        pass
    else:
        raise AssertionError("protocol 'auto' was accepted")


tests = [
    test_prepare_reference,
    test_read_qc,
    test_quantify,
    test_convert_matrices,
    test_multiqc_report,
    test_remove_background,
    test_scrnaseq,
    test_scrnaseq_options,
    test_protocols,
]
