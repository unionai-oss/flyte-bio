"""scrnaseq — droplet single-cell RNA-seq quantification with simpleaf (alevin-fry).

Ports the default ``simpleaf`` route of the upstream single-cell RNA-seq
pipeline (4.2.0):

- :func:`samples_from_samplesheet` — group samplesheet rows (``sample,
  fastq_1, fastq_2``, optionally ``expected_cells`` and ``seq_center``; one
  per lane) into samples.
- :func:`prepare_reference` — decompress the genome and annotation and keep
  the GTF lines on the genome's sequences, as upstream.
- :func:`read_qc` — FastQC over all of a sample's FASTQs.
- :func:`prepare_index` — the simpleaf (piscem) index of the spliced +
  intronic transcriptome, unless one is given.
- :func:`quantify` — simpleaf quant per sample with the protocol's
  chemistry and barcode whitelist.
- :func:`cell_qc` — QCatch cell calling and QC report on a sample's quantification.
- :func:`convert_matrices` — each sample's counts as h5ad, all of them
  combined, and every h5ad as Seurat and SingleCellExperiment RDS files.
- :func:`remove_background` — CellBender cell calling (CPU or GPU) on the
  raw matrices, which are restricted to its cells and converted likewise.
- :func:`multiqc_report` — one MultiQC report over the FastQC reports and
  simpleaf quantifications.
- :func:`scrnaseq` — the whole route from a samplesheet, optionally
  published in upstream's ``--outdir`` layout (:func:`scrnaseq_layout`).

These are plain async functions, not tasks: they run inside the caller's
task and fan out to the module tasks, so the caller's
:class:`flyte.TaskEnvironment` must ``depends_on`` :data:`flyte_bio.modules.env`.
"""

import asyncio
import re
import tempfile
import urllib.request
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path

import flyte
from flyte.io import Dir, File

from flyte_bio.modules.anndata import anndata_barcodes
from flyte_bio.modules.anndatar import RdsObjects, anndatar_convert
from flyte_bio.modules.cellbender import CellbenderResult, cellbender_removebackground
from flyte_bio.modules.fastqc import fastqc_reads
from flyte_bio.modules.gtfgenefilter import gtf_gene_filter
from flyte_bio.modules.gunzip import gunzip
from flyte_bio.modules.multiqc import MultiQCResult, multiqc
from flyte_bio.modules.qcatch import QcatchResult, qcatch
from flyte_bio.modules.scanpy import concat_h5ad, mtx_to_h5ad_simpleaf
from flyte_bio.modules.simpleaf import SimpleafIndex, SimpleafQuant, simpleaf_index, simpleaf_quant
from flyte_bio.publish import expand, name, publish
from flyte_bio.samplesheet import Row, read_samplesheet
from flyte_bio.scripts import path as asset


@dataclass
class Sample:
    """One sample: its lanes' paired FASTQs (R1 barcodes + UMIs, R2 cDNA), in samplesheet order."""

    id: str
    fastq_1: list[File]
    fastq_2: list[File]
    expected_cells: int | None = None
    seq_center: str = ""


def samples_from_samplesheet(rows: list[Row]) -> list[Sample]:
    """Group samplesheet rows (one per lane) into Samples, in file order, as upstream.

    Upstream's columns: ``sample``, ``fastq_1`` and ``fastq_2`` (both
    required: R1 carries the cell barcode and UMI), and optionally
    ``expected_cells`` and ``seq_center``. A sample's rows must agree on the
    optional columns.
    """
    groups: dict[str, list[Row]] = {}
    for row in rows:
        groups.setdefault(row["sample"], []).append(row)
    samples = []
    for sample_id, sample_rows in groups.items():
        lines = [r.line for r in sample_rows]
        if any(c.isspace() for c in sample_id):
            raise ValueError(f"sample {sample_id!r} (lines {lines}): sample IDs must not contain whitespace")
        for r in sample_rows:
            if not r["fastq_1"] or not r["fastq_2"]:
                raise ValueError(f"sample {sample_id!r} (line {r.line}): needs both fastq_1 and fastq_2")
        for column in ("expected_cells", "seq_center"):
            values = {r.get(column, "") for r in sample_rows}
            if len(values) > 1:
                raise ValueError(f"sample {sample_id!r} (lines {lines}): rows disagree on {column} {sorted(values)}")
        expected = sample_rows[0].get("expected_cells", "")
        samples.append(
            Sample(
                id=sample_id,
                fastq_1=[File.from_existing_remote(r["fastq_1"]) for r in sample_rows],
                fastq_2=[File.from_existing_remote(r["fastq_2"]) for r in sample_rows],
                expected_cells=int(expected) if expected else None,
                seq_center=sample_rows[0].get("seq_center", ""),
            )
        )
    return samples


@dataclass
class Reference:
    fasta: File
    gtf: File  # filtered to the genome's sequences: <genome name>_genes.gtf


async def maybe_gunzip(file: File) -> File:
    """Decompress ``file`` when its name ends in ``.gz``, as upstream decides."""
    return await gunzip(archive=file) if file.path.endswith(".gz") else file


async def prepare_reference(fasta: File, gtf: File) -> Reference:
    """Upstream's genome prep: decompress the FASTA and GTF, keep the GTF lines on the genome's sequences."""
    fasta, gtf = await asyncio.gather(maybe_gunzip(fasta), maybe_gunzip(gtf))
    return Reference(fasta=fasta, gtf=await gtf_gene_filter(fasta, gtf))


async def read_qc(sample: Sample) -> Dir:
    """FastQC over all of a sample's FASTQs, lane by lane (R1, R2, ...), as ``<sample>_<n>_fastqc.*``."""
    reads = [f for pair in zip(sample.fastq_1, sample.fastq_2) for f in pair]
    return await fastqc_reads(reads=reads, prefix=sample.id, args="--quiet")


ASSETS = "https://raw.githubusercontent.com/nf-core/scrnaseq/4.2.0/assets/"


@dataclass
class Protocol:
    chemistry: str  # simpleaf's --chemistry
    whitelist: str | None  # barcode whitelist under ASSETS
    qcatch: str | None  # QCatch's --chemistry


# Upstream's simpleaf protocols.
PROTOCOLS = {
    "10XV1": Protocol("10xv1", "whitelist/10x_V1_barcode_whitelist.txt.gz", None),
    "10XV2": Protocol("10xv2", "whitelist/10x_V2_barcode_whitelist.txt.gz", "10X_3p_v2"),
    "10XV3": Protocol("10xv3", "whitelist/10x_V3_barcode_whitelist.txt.gz", "10X_3p_v3"),
    "10XV4": Protocol("10xv4-3p", "whitelist/10x_V4_barcode_whitelist.txt.gz", "10X_3p_v4"),
    "dropseq": Protocol("dropseq", None, None),
}


def get_protocol(name: str) -> Protocol:
    """Upstream's protocol for ``name``; any other name goes to simpleaf as its chemistry, unmodified, as upstream.

    That covers simpleaf's own chemistry names and custom read geometries.
    Upstream's ``auto`` is cellranger-only, so it's refused here.
    """
    if name == "auto":
        raise ValueError(f"protocol 'auto' is only supported by cellranger: give simpleaf's, e.g. {sorted(PROTOCOLS)}")
    if name not in PROTOCOLS:
        flyte.logger.warning(f"Protocol {name!r} not recognized: passing it on to simpleaf unmodified")
        return Protocol(name, None, None)
    return PROTOCOLS[name]


async def get_barcode_whitelist(proto: Protocol, override: File | None = None) -> File | None:
    """``override``, else the protocol's whitelist from upstream's assets (None for protocols without one)."""
    if override is not None or proto.whitelist is None:
        return override
    local = Path(tempfile.mkdtemp()) / proto.whitelist.rsplit("/", 1)[-1]

    def download() -> None:
        with urllib.request.urlopen(ASSETS + proto.whitelist, timeout=300) as response, local.open("wb") as out:
            while chunk := response.read(1 << 20):
                out.write(chunk)

    await asyncio.to_thread(download)
    return await File.from_local(str(local))


async def prepare_index(
    reference: Reference,
    index: Dir | None = None,
    txp2gene: File | None = None,
    transcript_fasta: File | None = None,
) -> SimpleafIndex:
    """Upstream's index step: ``index`` as given, else built from ``transcript_fasta`` or the genome."""
    if index is not None:
        return SimpleafIndex(index=index, ref=None, t2g=txp2gene)
    if transcript_fasta is not None:
        if txp2gene is None:
            raise ValueError("txp2gene is required when building the index from transcript_fasta")
        built = await simpleaf_index(prefix="simpleaf_index", transcript_fasta=transcript_fasta)
        return SimpleafIndex(index=built.index, ref=None, t2g=txp2gene)
    built = await simpleaf_index(prefix="simpleaf_index", genome_fasta=reference.fasta, gtf=reference.gtf)
    return SimpleafIndex(index=built.index, ref=built.ref, t2g=txp2gene or built.t2g)


async def quantify(
    sample: Sample,
    index: SimpleafIndex,
    proto: Protocol,
    whitelist: File | None,
    resolution: str = "cr-like",
) -> SimpleafQuant:
    """simpleaf quant of one sample, unfiltered against the barcode whitelist, as upstream."""
    return await simpleaf_quant(
        reads_1=sample.fastq_1,
        reads_2=sample.fastq_2,
        index=index.index,
        t2g=index.t2g,
        chemistry=proto.chemistry,
        prefix=sample.id,
        resolution=resolution,
        cell_filter="unfiltered-pl",
        permit_list=whitelist,
    )


async def cell_qc(
    sample: Sample,
    quant: SimpleafQuant,
    proto: Protocol,
    remove_doublets: bool = False,
    n_partitions: int | None = None,
) -> QcatchResult:
    """QCatch over a sample's simpleaf quantification, with the protocol's QCatch chemistry, as upstream.

    QCatch's cell calling needs a realistic number of barcodes: on a toy
    dataset (such as upstream's test profile, which skips QCatch) it fails.
    """
    return await qcatch(
        quant.quant,
        prefix=sample.id,
        chemistry=proto.qcatch,
        remove_doublets=remove_doublets,
        n_partitions=n_partitions,
    )


@dataclass
class Matrix:
    """One count matrix: the h5ad and its R conversions, named ``<id>_<input_type>_matrix.*``."""

    h5ad: File
    rds: RdsObjects


@dataclass
class Matrices:
    samples: dict[str, Matrix]  # by sample ID
    combined: Matrix  # every sample's cells, annotated from the samplesheet


async def convert_matrices(quants: dict[str, SimpleafQuant], samplesheet: File) -> Matrices:
    """Upstream's h5ad conversion of the simpleaf counts (``raw``: not cell-filtered).

    Each sample's ``af_quant`` becomes ``<sample>_raw_matrix.h5ad``; those are
    concatenated into ``combined_raw_matrix.h5ad``, each cell annotated with
    its sample's samplesheet columns; and every h5ad is also written as
    ``.seurat.rds`` and ``.sce.rds``.
    """
    ids = list(quants)
    h5ads = await asyncio.gather(*(mtx_to_h5ad_simpleaf(quants[s].quant, s, "raw") for s in ids))
    return await matrices_from_h5ads(dict(zip(ids, h5ads)), samplesheet, "raw")


async def matrices_from_h5ads(h5ads: dict[str, File], samplesheet: File, input_type: str) -> Matrices:
    """Per-sample ``<sample>_<input_type>_matrix.h5ad`` files, combined and each converted to RDS."""
    ids = list(h5ads)
    combined = await concat_h5ad(list(h5ads.values()), samplesheet, prefix=f"combined_{input_type}")
    all_rds = await asyncio.gather(
        *(anndatar_convert(h5ads[s], f"{s}_{input_type}_matrix") for s in ids),
        anndatar_convert(combined, f"combined_{input_type}_matrix"),
    )
    return Matrices(
        samples={s: Matrix(h5ad=h5ads[s], rds=r) for s, r in zip(ids, all_rds)},
        combined=Matrix(h5ad=combined, rds=all_rds[-1]),
    )


async def multiqc_report(fastqc: dict[str, Dir], quants: dict[str, SimpleafQuant]) -> MultiQCResult:
    """One MultiQC report over the samples' FastQC reports and simpleaf quantifications, with upstream's config.

    As upstream stages them: every FastQC zip and HTML, and each sample's
    ``af_quant`` in a directory of its own. Of ``af_quant``, only its
    top-level JSON summaries are passed, not the count matrices.
    """
    root = Path(tempfile.mkdtemp(prefix="multiqc_"))
    downloads = []

    async def put_dir(d: Dir, rel: str, keep: Callable[[str], bool]) -> None:
        base = d.path.rstrip("/")
        async for f in d.walk():
            sub = f.path[len(base) :].lstrip("/")
            if keep(sub):
                (root / rel / sub).parent.mkdir(parents=True, exist_ok=True)
                downloads.append(f.download(str(root / rel / sub)))

    def fastqc_report(sub: str) -> bool:
        return sub.endswith(("_fastqc.zip", "_fastqc.html"))

    def summary_json(sub: str) -> bool:
        return "/" not in sub and sub.endswith(".json")

    for sample_id, reports in fastqc.items():
        await put_dir(reports, f"fastqc/{sample_id}", fastqc_report)
    for sample_id, quant in quants.items():
        await put_dir(quant.quant, f"simpleaf/{sample_id}/af_quant", summary_json)
    await asyncio.gather(*downloads)
    data = await Dir.from_local(str(root))
    config = await File.from_local(str(asset("multiqc_config_scrnaseq.yml")))
    return await multiqc(data, configs=[config], prefix="multiqc_report")


@dataclass
class BackgroundRemoval:
    cellbender: dict[str, CellbenderResult]  # by sample ID
    matrices: Matrices  # the raw h5ads restricted to CellBender's cells: <sample>_cellbender_filter_matrix.*


async def remove_background(
    raw: dict[str, File],
    samplesheet: File,
    gpu: bool = False,
    cellbender_args: str = "",
) -> BackgroundRemoval:
    """Upstream's CellBender step on each sample's raw h5ad (by sample ID).

    CellBender calls the cells (on a GPU when ``gpu``); each raw h5ad is then
    restricted to them as ``<sample>_cellbender_filter_matrix.h5ad``, which
    are combined and converted to RDS like the raw matrices. The counts stay
    raw: CellBender's denoised counts are only in its ``.h5`` outputs.
    """
    ids = list(raw)
    results = await asyncio.gather(
        *(cellbender_removebackground(raw[s], s, gpu=gpu, args=cellbender_args) for s in ids)
    )
    filtered = await asyncio.gather(
        *(anndata_barcodes(raw[s], r.barcodes, f"{s}_cellbender_filter_matrix") for s, r in zip(ids, results))
    )
    return BackgroundRemoval(
        cellbender=dict(zip(ids, results)),
        matrices=await matrices_from_h5ads(dict(zip(ids, filtered)), samplesheet, "cellbender_filter"),
    )


@dataclass
class SampleResult:
    sample: Sample
    fastqc: Dir | None  # None when skip_fastqc
    quant: SimpleafQuant
    qcatch: QcatchResult | None  # None when skip_qcatch


@dataclass
class ScrnaseqResult:
    reference: Reference
    index: SimpleafIndex
    samples: list[SampleResult]
    raw: Matrices  # <sample>_raw_matrix.* and combined_raw_matrix.*
    cellbender: BackgroundRemoval | None = None  # None when skip_cellbender
    multiqc: MultiQCResult | None = None  # None when skip_multiqc
    outdir: Dir | None = None  # the published results tree, when outdir/publish_results was asked for


async def scrnaseq(
    samplesheet: File,
    fasta: File,
    gtf: File,
    protocol: str,
    barcode_whitelist: File | None = None,
    simpleaf_index: Dir | None = None,
    txp2gene: File | None = None,
    transcript_fasta: File | None = None,
    simpleaf_umi_resolution: str = "cr-like",
    skip_fastqc: bool = False,
    skip_qcatch: bool = False,
    remove_doublets: bool = False,
    qcatch_n_partitions: int | None = None,
    skip_cellbender: bool = False,
    cellbender_gpu: bool = False,
    cellbender_args: str = "",
    skip_multiqc: bool = False,
    save_reference: bool = False,
    save_align_intermeds: bool = True,
    outdir: str | None = None,
    publish_results: bool = False,
) -> ScrnaseqResult:
    """Droplet scRNA-seq from FASTQ to count matrices with simpleaf, as upstream's default (``simpleaf``) route.

    ``samplesheet`` is upstream's ``--input``: a CSV with columns ``sample``,
    ``fastq_1`` (cell barcode + UMI) and ``fastq_2`` (cDNA), optionally
    ``expected_cells`` and ``seq_center``, one row per lane. It is read
    before anything runs. ``protocol`` is required, as upstream's ``auto``
    is cellranger-only: one of upstream's (``10XV1``..``10XV4``,
    ``dropseq``), or else passed to simpleaf as its chemistry.

    The reference is prepared and indexed once (or ``simpleaf_index`` is
    used, with ``txp2gene``; ``transcript_fasta`` + ``txp2gene`` build the
    index from transcripts instead of the genome). Every sample is then
    FastQC'd and quantified against the protocol's barcode whitelist (or
    ``barcode_whitelist``) in parallel, and cell-QC'd with QCatch. The
    raw counts are converted to h5ad and RDS, per sample and combined;
    CellBender (``cellbender_gpu`` runs it on a GPU) calls cells on them,
    and its cells' counts are converted likewise. MultiQC reports on all of
    it.

    QCatch and CellBender need realistic data: upstream's test profile skips
    both (``skip_qcatch``, ``skip_cellbender``).

    Outputs: everything is returned in the :class:`ScrnaseqResult`. Like
    upstream's ``--outdir``, ``outdir`` (e.g. ``s3://my-bucket/scrnaseq``)
    also lays the results out there in upstream's folder structure
    (:func:`scrnaseq_layout`); ``publish_results=True`` does the same in
    Flyte's own storage. ``save_reference`` also publishes the filtered GTF
    and the built index, and ``save_align_intermeds=False`` leaves out
    simpleaf's mapping (``af_map``), as upstream.
    """
    samples = samples_from_samplesheet(await read_samplesheet(samplesheet))
    if not samples:
        raise ValueError("the samplesheet has no samples")
    proto = get_protocol(protocol)

    reference = await prepare_reference(fasta, gtf)
    index, whitelist = await asyncio.gather(
        prepare_index(reference, index=simpleaf_index, txp2gene=txp2gene, transcript_fasta=transcript_fasta),
        get_barcode_whitelist(proto, barcode_whitelist),
    )

    async def fastqc(sample: Sample) -> Dir | None:
        return None if skip_fastqc else await read_qc(sample)

    async def run_sample(sample: Sample) -> SampleResult:
        reports, quant = await asyncio.gather(
            fastqc(sample),
            quantify(sample, index, proto, whitelist, resolution=simpleaf_umi_resolution),
        )
        qc = None
        if not skip_qcatch:
            qc = await cell_qc(sample, quant, proto, remove_doublets=remove_doublets, n_partitions=qcatch_n_partitions)
        return SampleResult(sample, reports, quant, qc)

    results = list(await asyncio.gather(*(run_sample(s) for s in samples)))
    ids = [r.sample.id for r in results]
    raw = await convert_matrices({r.sample.id: r.quant for r in results}, samplesheet)
    result = ScrnaseqResult(reference=reference, index=index, samples=results, raw=raw)

    async def background() -> None:
        if not skip_cellbender:
            result.cellbender = await remove_background(
                {s: raw.samples[s].h5ad for s in ids}, samplesheet, gpu=cellbender_gpu, cellbender_args=cellbender_args
            )

    async def report() -> None:
        if not skip_multiqc:
            fastqc = {r.sample.id: r.fastqc for r in results if r.fastqc is not None}
            result.multiqc = await multiqc_report(fastqc, {r.sample.id: r.quant for r in results})

    await asyncio.gather(background(), report())
    if outdir is not None or publish_results:
        layout = await scrnaseq_layout(result, save_reference=save_reference, save_align_intermeds=save_align_intermeds)
        result.outdir = await publish(layout, outdir)
    return result


async def scrnaseq_layout(
    result: ScrnaseqResult, save_reference: bool = False, save_align_intermeds: bool = True
) -> dict[str, File | Dir | None]:
    """Upstream's results layout for a run: relative path -> output, for :func:`~flyte_bio.publish.publish`.

    Follows upstream's default ``publishDir`` paths: ``fastqc/``,
    ``simpleaf/<sample>/`` (``simpleaf_quant/``, ``qcatch/``,
    ``cellbender_removebackground/``), ``simpleaf/mtx_conversions/`` and
    ``multiqc/``. ``save_reference`` adds ``gtf_filter/`` and the built
    ``simpleaf/simpleaf_index/``; without ``save_align_intermeds``,
    ``simpleaf_quant/af_map/`` is left out.
    """
    out: dict[str, File | Dir | None] = {}
    if save_reference:
        out[f"gtf_filter/{name(result.reference.gtf)}"] = result.reference.gtf
        if result.index.ref is not None:  # built here, not given
            out.update(await expand("simpleaf/simpleaf_index/index", result.index.index))
            out.update(await expand("simpleaf/simpleaf_index/ref", result.index.ref))
    for r in result.samples:
        sid = r.sample.id
        if r.fastqc is not None:
            out.update(await expand("fastqc", r.fastqc))
        quant_dir = f"simpleaf/{sid}/simpleaf_quant"
        if save_align_intermeds:
            out.update(await expand(f"{quant_dir}/af_map", r.quant.map))
        out.update(await expand(f"{quant_dir}/af_quant", r.quant.quant))
        if r.qcatch is not None:
            for f in (r.qcatch.report, r.qcatch.filtered_h5ad, r.qcatch.metrics_summary):
                out[f"simpleaf/{sid}/qcatch/{name(f)}"] = f

    def put_matrices(m: Matrices) -> None:
        for sid, matrix in m.samples.items():
            for f in (matrix.h5ad, matrix.rds.sce, matrix.rds.seurat):
                out[f"simpleaf/mtx_conversions/{sid}/{name(f)}"] = f
        for f in (m.combined.h5ad, m.combined.rds.sce, m.combined.rds.seurat):
            out[f"simpleaf/mtx_conversions/{name(f)}"] = f

    put_matrices(result.raw)
    if result.cellbender is not None:
        for sid, cb in result.cellbender.cellbender.items():
            for f in (
                cb.h5,
                cb.filtered_h5,
                cb.posterior_h5,
                cb.barcodes,
                cb.metrics,
                cb.report,
                cb.pdf,
                cb.log,
                cb.checkpoint,
            ):
                if f is not None:
                    out[f"simpleaf/{sid}/cellbender_removebackground/{name(f)}"] = f
        put_matrices(result.cellbender.matrices)
    if result.multiqc is not None:
        # Upstream runs MultiQC without --filename, so its folders get MultiQC's default names.
        for path, f in (await expand("multiqc", result.multiqc.results)).items():
            out[re.sub(r"^multiqc/multiqc_report_(data|plots)/", r"multiqc/multiqc_\1/", path)] = f
    return out
