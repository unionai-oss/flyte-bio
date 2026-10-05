"""tx2gene — map transcripts to genes from a GTF.

Reads the top transcripts from one sample's quantification output to
discover which GTF attribute holds the transcript IDs, then writes a
``transcript → gene (+ extra attributes)`` table. Pure Python, run as a
shell task on a stock python biocontainer with the vendored script staged in
as a File input — so the tool pod needs no ``flyte_bio`` install.

The table is written under its natural ``<prefix>.tx2gene.tsv`` name inside
an output Dir: downstream readers (summarizedexperiment) pick the delimiter
from the file extension, which a bare output path wouldn't have.
"""

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.scripts import path

TX2GENE_IMAGE = "quay.io/biocontainers/python:3.10.4"

DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory="6Gi")

tx2gene_cmd = shell.create(
    name="tx2gene",
    image=TX2GENE_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "script": File,
        "gtf": File,
        "quants": Dir,
        "quant_type": str,
        "gene_id": str,
        "extra": str,
        "filename": str,
    },
    outputs={"results": Dir},
    script=(
        "python {inputs.script} --quant-type {inputs.quant_type} --gtf {inputs.gtf} "
        "--quants {inputs.quants} --id {inputs.gene_id} --extra {inputs.extra} "
        "--output {outputs.results}/{inputs.filename}\n"
    ),
)


env = flyte.TaskEnvironment.from_task(
    "tx2gene",
    tx2gene_cmd.as_task(),
)


async def tx2gene(
    gtf: File,
    quants: Dir,
    quant_type: str = "salmon",
    gene_id: str = "gene_id",
    extra: str = "gene_name",
    prefix: str = "",
) -> File:
    """Return the ``transcript → gene`` TSV.

    ``quants`` holds one sample's quant file (``quant.sf`` for salmon) or
    per-sample subdirectories of them; any one sample suffices.
    ``extra`` may list several attributes, comma-separated.
    """
    filename = f"{prefix}.tx2gene.tsv" if prefix else "tx2gene.tsv"
    script = await File.from_local(str(path("tx2gene.py")))
    results = await tx2gene_cmd(
        script=script,
        gtf=gtf,
        quants=quants,
        quant_type=quant_type,
        gene_id=gene_id,
        extra=extra,
        filename=filename,
    )
    table = await results.get_file(filename)
    if table is None:
        raise FileNotFoundError(f"tx2gene wrote no {filename}")
    return table
