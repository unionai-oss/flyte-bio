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

These are plain async functions, not tasks: they run inside the caller's
task and fan out to the module tasks, so the caller's
:class:`flyte.TaskEnvironment` must ``depends_on`` :data:`flyte_bio.modules.env`.
"""

import asyncio
import tempfile
import urllib.request
from dataclasses import dataclass
from pathlib import Path

from flyte.io import Dir, File

from flyte_bio.modules.anndatar import RdsObjects, anndatar_convert
from flyte_bio.modules.fastqc import fastqc_reads
from flyte_bio.modules.gtfgenefilter import gtf_gene_filter
from flyte_bio.modules.gunzip import gunzip
from flyte_bio.modules.qcatch import QcatchResult, qcatch
from flyte_bio.modules.scanpy import concat_h5ad, mtx_to_h5ad_simpleaf
from flyte_bio.modules.simpleaf import SimpleafIndex, SimpleafQuant, simpleaf_index, simpleaf_quant
from flyte_bio.samplesheet import Row


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


def protocol(name: str) -> Protocol:
    if name not in PROTOCOLS:
        raise ValueError(f"protocol must be one of {sorted(PROTOCOLS)} for simpleaf, got {name!r}")
    return PROTOCOLS[name]


async def barcode_whitelist(proto: Protocol, override: File | None = None) -> File | None:
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
    combined = await concat_h5ad(list(h5ads), samplesheet, prefix="combined_raw")
    all_rds = await asyncio.gather(
        *(anndatar_convert(h, f"{s}_raw_matrix") for s, h in zip(ids, h5ads)),
        anndatar_convert(combined, "combined_raw_matrix"),
    )
    return Matrices(
        samples={s: Matrix(h5ad=h, rds=r) for s, h, r in zip(ids, h5ads, all_rds)},
        combined=Matrix(h5ad=combined, rds=all_rds[-1]),
    )
