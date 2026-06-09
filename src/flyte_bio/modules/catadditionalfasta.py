"""catadditionalfasta — append an extra FASTA onto a genome reference.

Concatenates an additional FASTA (e.g. spike-ins / transgenes) onto the
genome FASTA, and a generated GTF describing it onto the genome GTF. Pure
Python, but it runs as a shell task on a stock python biocontainer with the
vendored script staged in as a File input — so the tool pod is
self-contained and needs no ``flyte_bio`` install.
"""



import flyte
from flyte.extras import shell
from flyte.io import File

from flyte_bio.scripts import path

# Pinned biocontainer URI (a bare Python interpreter).
CATADDITIONALFASTA_IMAGE = "quay.io/biocontainers/python:3.12"

DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory="6Gi")

cat_additional_fasta_cmd = shell.create(
    name="cat_additional_fasta",
    image=CATADDITIONALFASTA_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"script": File, "fasta": File, "gtf": File, "add_fasta": File, "biotype": str},
    outputs={"out_fasta": File, "out_gtf": File},
    script=(
        "python {inputs.script} --fasta {inputs.fasta} --gtf {inputs.gtf} "
        "--add-fasta {inputs.add_fasta} --biotype {inputs.biotype} "
        "--out-fasta {outputs.out_fasta} --out-gtf {outputs.out_gtf}\n"
    ),
)


env = flyte.TaskEnvironment.from_task(
    "catadditionalfasta",
    cat_additional_fasta_cmd.as_task(),
)


async def cat_additional_fasta(fasta: File, gtf: File, add_fasta: File, biotype: str = "") -> tuple[File, File]:
    """Return ``(genome+add FASTA, genome+generated GTF)``."""
    script = await File.from_local(str(path("catadditionalfasta.py")))
    return await cat_additional_fasta_cmd(
        script=script, fasta=fasta, gtf=gtf, add_fasta=add_fasta, biotype=biotype
    )
