"""catadditionalfasta — append an extra FASTA onto a genome reference.

Concatenates an additional FASTA (e.g. spike-ins / transgenes) onto the
genome FASTA, and a generated GTF describing it onto the genome GTF. Pure
Python, so it's a native Flyte task importing the vendored logic directly —
no biocontainer, no shell.
"""



from pathlib import Path

import flyte
from flyte.io import File

from flyte_bio.scripts.catadditionalfasta import cat_additional_fasta as cat_additional_fasta_impl

DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory="6Gi")

env = flyte.TaskEnvironment(
    name="catadditionalfasta",
    image=flyte.Image.from_debian_base(),
    resources=DEFAULT_RESOURCES,
)


@env.task
async def cat_additional_fasta(fasta: File, gtf: File, add_fasta: File, biotype: str = "") -> tuple[File, File]:
    """Return ``(genome+add FASTA, genome+generated GTF)``."""
    fasta_local = await fasta.download(Path("in") / (fasta.name or "genome.fasta"))
    gtf_local = await gtf.download(Path("in") / (gtf.name or "genome.gtf"))
    add_local = await add_fasta.download(Path("in") / (add_fasta.name or "add.fasta"))

    cat_additional_fasta_impl(fasta_local, gtf_local, add_local, biotype, "out.fasta", "out.gtf")
    return await File.from_local("out.fasta"), await File.from_local("out.gtf")
