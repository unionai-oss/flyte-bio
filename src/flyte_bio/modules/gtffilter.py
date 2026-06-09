"""gtffilter — restrict a GTF to a genome's sequences.

Keeps only annotation records whose sequence name is present in the genome
FASTA (and, unless skipped, records carrying a ``transcript_id``). The
filtering is pure-stdlib Python, so this is a native Flyte task that imports
the vendored logic and calls it directly — no biocontainer, no shell. The
task runs on a plain Flyte base image; ``flyte_bio`` (and the helper) ride
in via the code bundle.
"""



from pathlib import Path

import flyte
from flyte.io import File

from flyte_bio.scripts.gtffilter import filter_gtf

DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory="6Gi")

env = flyte.TaskEnvironment(
    name="gtffilter",
    image=flyte.Image.from_debian_base(),
    resources=DEFAULT_RESOURCES,
)


@env.task
async def gtf_filter(gtf: File, fasta: File) -> File:
    """Filter ``gtf`` to the sequences present in ``fasta``."""
    # Download under the original names so filter_gtf's `.gz` sniffing
    # (it keys off the path suffix) still works for gzipped inputs.
    gtf_local = await gtf.download(Path("in") / (gtf.name or "annotation.gtf"))
    fasta_local = await fasta.download(Path("in") / (fasta.name or "genome.fasta"))

    suffix = ".gtf.gz" if gtf_local.endswith(".gz") else ".gtf"
    out = f"filtered{suffix}"
    filter_gtf(fasta_local, gtf_local, out, False)
    return await File.from_local(out)
