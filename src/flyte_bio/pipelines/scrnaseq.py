"""scrnaseq — droplet single-cell RNA-seq quantification with simpleaf (alevin-fry).

Ports the default ``simpleaf`` route of the upstream single-cell RNA-seq
pipeline (4.2.0):

- :func:`samples_from_samplesheet` — group samplesheet rows (``sample,
  fastq_1, fastq_2``, optionally ``expected_cells`` and ``seq_center``; one
  per lane) into samples.
- :func:`prepare_reference` — decompress the genome and annotation and keep
  the GTF lines on the genome's sequences, as upstream.
- :func:`read_qc` — FastQC over all of a sample's FASTQs.

These are plain async functions, not tasks: they run inside the caller's
task and fan out to the module tasks, so the caller's
:class:`flyte.TaskEnvironment` must ``depends_on`` :data:`flyte_bio.modules.env`.
"""

import asyncio
from dataclasses import dataclass

from flyte.io import Dir, File

from flyte_bio.modules.fastqc import fastqc_reads
from flyte_bio.modules.gtfgenefilter import gtf_gene_filter
from flyte_bio.modules.gunzip import gunzip
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
