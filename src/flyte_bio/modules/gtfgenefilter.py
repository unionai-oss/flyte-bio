"""gtfgenefilter — keep the GTF lines on a genome's sequences.

Exposes :func:`gtf_gene_filter`, which keeps every GTF line whose sequence
is in the FASTA, writing ``<fasta name without extension>_genes.gtf``, as the
upstream scRNA-seq pipeline does before building its index. Unlike
:mod:`flyte_bio.modules.gtffilter` (the RNA-seq pipeline's filter) it keeps
lines without a ``transcript_id``, such as gene records.

The script is upstream's ``filter_gtf_for_genes_in_genome.py``, vendored
unchanged and staged into a stock Python image as a File input.
"""

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.scripts import path

# Upstream names the bare `biocontainers/python:3.9--1` (Docker Hub); use the quay.io mirror.
PYTHON_IMAGE = "quay.io/biocontainers/python:3.9--1"

DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory="4Gi")

gtf_gene_filter_cmd = shell.create(
    name="gtf_gene_filter",
    image=PYTHON_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"script": File, "fasta": File, "gtf": File},
    outputs={"results": Dir},
    script=r"""
        S=({inputs.script}); FA=({inputs.fasta}); G=({inputs.gtf})
        NAME=$(basename "${FA[0]}")
        python "${S[0]}" --gtf "${G[0]}" --fasta "${FA[0]}" -o {outputs.results}/"${NAME%.*}_genes.gtf"
    """,
)


env = flyte.TaskEnvironment.from_task(
    "gtfgenefilter",
    gtf_gene_filter_cmd.as_task(),
)


async def gtf_gene_filter(fasta: File, gtf: File) -> File:
    """``<fasta name>_genes.gtf``: the lines of ``gtf`` on sequences in ``fasta`` (both uncompressed)."""
    script = await File.from_local(str(path("filter_gtf_for_genes_in_genome.py")))
    results = await gtf_gene_filter_cmd(script=script, fasta=fasta, gtf=gtf)
    found = [f async for f in results.walk() if f.path.endswith("_genes.gtf")]
    if len(found) != 1:
        raise FileNotFoundError(f"gtf_gene_filter wrote {len(found)} *_genes.gtf files, expected one")
    return found[0]
