"""gtffilter — restrict a GTF to a genome's sequences.

Keeps only annotation records whose sequence name is present in the genome
FASTA (and, unless skipped, records carrying a ``transcript_id``). The
filtering is pure-stdlib Python, but it runs as a shell task on a stock
python biocontainer with the vendored script staged in as a File input —
so the tool pod is self-contained and needs no ``flyte_bio`` install. The
vendored script rides Flyte's default code bundle via the stub-module trick
in :mod:`flyte_bio.scripts`.
"""



import flyte
from flyte.extras import shell
from flyte.io import File

from flyte_bio.scripts import path

# Pinned biocontainer URI (a bare Python interpreter).
GTFFILTER_IMAGE = "quay.io/biocontainers/python:3.12"

DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory="6Gi")

gtf_filter_cmd = shell.create(
    name="gtf_filter",
    image=GTFFILTER_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"script": File, "gtf": File, "fasta": File},
    outputs={"gtf_filtered": File},
    script="python {inputs.script} --gtf {inputs.gtf} --fasta {inputs.fasta} --output {outputs.gtf_filtered}\n",
)


env = flyte.TaskEnvironment.from_task(
    "gtffilter",
    gtf_filter_cmd.as_task(),
)


async def gtf_filter(gtf: File, fasta: File) -> File:
    """Filter ``gtf`` to the sequences present in ``fasta``."""
    script = await File.from_local(str(path("gtffilter.py")))
    return await gtf_filter_cmd(script=script, gtf=gtf, fasta=fasta)
