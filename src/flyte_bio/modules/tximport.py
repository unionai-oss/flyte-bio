"""tximport — summarize transcript quantifications into count/TPM matrices.

Imports every sample's quantification output (salmon, kallisto or RSEM) with
tximport and writes transcript- and gene-level count, TPM and length tables,
plus the tx2gene table actually used. The R helper is vendored and staged
into the stock tximeta biocontainer as a File input.

Sample names come from the directory holding each quant file
(``<sample>/quant.sf``), so :func:`tximport` takes one Dir laid out that
way. Shell tasks can't take a ``list[Dir]``/``dict[str, Dir]``, so
:func:`collect_quants` builds that Dir from per-sample result Dirs.
"""

import asyncio
import shutil
import tempfile
from dataclasses import dataclass
from pathlib import Path

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.scripts import path

TXIMPORT_IMAGE = "quay.io/biocontainers/bioconductor-tximeta:1.20.1--r43hdfd78af_0"

DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory="6Gi")

tximport_cmd = shell.create(
    name="tximport",
    image=TXIMPORT_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"script": File, "quants": Dir, "tx2gene": File, "quant_type": str, "prefix": str},
    outputs={"results": Dir},
    script=(
        "Rscript {inputs.script} --quant_type {inputs.quant_type} --tx2gene {inputs.tx2gene} "
        "--quants {inputs.quants} --prefix {outputs.results}/{inputs.prefix}\n"
    ),
)


env = flyte.TaskEnvironment.from_task(
    "tximport",
    tximport_cmd.as_task(),
)


@dataclass
class TximportResult:
    tpm_gene: File
    counts_gene: File
    counts_gene_length_scaled: File
    counts_gene_scaled: File
    lengths_gene: File
    tpm_transcript: File
    counts_transcript: File
    lengths_transcript: File
    tx2gene_augmented: File


# TximportResult field -> output file suffix (after "<prefix>.").
OUTPUT_SUFFIXES = {
    "tpm_gene": "gene_tpm.tsv",
    "counts_gene": "gene_counts.tsv",
    "counts_gene_length_scaled": "gene_counts_length_scaled.tsv",
    "counts_gene_scaled": "gene_counts_scaled.tsv",
    "lengths_gene": "gene_lengths.tsv",
    "tpm_transcript": "transcript_tpm.tsv",
    "counts_transcript": "transcript_counts.tsv",
    "lengths_transcript": "transcript_lengths.tsv",
    "tx2gene_augmented": "tx2gene_augmented.tsv",
}


async def collect_quants(results: dict[str, Dir]) -> Dir:
    """Lay per-sample quant result Dirs out as one ``<sample>/...`` Dir.

    Runs in the caller's pod: downloads each sample's results and uploads
    them together. Quant outputs are small (a few MB per sample).
    """
    for sample in results:
        if not sample or "/" in sample or sample in (".", ".."):
            raise ValueError(f"invalid sample name {sample!r}")
    root = Path(tempfile.mkdtemp(prefix="quants_"))
    try:
        await asyncio.gather(*(d.download(root / sample) for sample, d in results.items()))
        return await Dir.from_local(root)
    finally:
        shutil.rmtree(root, ignore_errors=True)


async def tximport(
    quants: Dir, tx2gene: File, quant_type: str = "salmon", prefix: str = "all_samples"
) -> TximportResult:
    """Run tximport over ``quants`` (one ``<sample>/`` subdir per sample)."""
    script = await File.from_local(str(path("tximport.r")))
    results = await tximport_cmd(script=script, quants=quants, tx2gene=tx2gene, quant_type=quant_type, prefix=prefix)

    async def pick(suffix: str) -> File:
        found = await results.get_file(f"{prefix}.{suffix}")
        if found is None:
            raise FileNotFoundError(f"tximport wrote no {prefix}.{suffix}")
        return found

    files = await asyncio.gather(*(pick(s) for s in OUTPUT_SUFFIXES.values()))
    return TximportResult(**dict(zip(OUTPUT_SUFFIXES, files)))
