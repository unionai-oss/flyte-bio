"""gffread — validate, filter, and convert GFF/GTF annotations.

Exposes the two conversions the RNA-seq genome-prep stage needs:

- :data:`gffread_gff_to_gtf` — convert a GFF3 annotation to GTF (``-T``).
- :data:`gffread_transcripts_fasta` — extract spliced transcript sequences
  from a genome FASTA guided by an annotation (``-w``).
"""



import flyte
from flyte.extras import shell
from flyte.io import File

# Pinned biocontainer URI.
GFFREAD_IMAGE = "quay.io/biocontainers/gffread:0.12.7--hdcf5f25_4"

DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory="6Gi")


gffread_gff_to_gtf = shell.create(
    name="gffread_gff_to_gtf",
    image=GFFREAD_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"gff": File},
    outputs={"gtf": File},
    script=r"""
        gffread {inputs.gff} -T -o {outputs.gtf}
    """,
)


# ``-w`` writes spliced exon (transcript) sequences; ``-g`` supplies the
# genome they're spliced out of. The annotation may be GFF or GTF.
gffread_transcripts_fasta = shell.create(
    name="gffread_transcripts_fasta",
    image=GFFREAD_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"annotation": File, "fasta": File},
    outputs={"transcripts": File},
    script=r"""
        gffread {inputs.annotation} -g {inputs.fasta} -w {outputs.transcripts}
    """,
)


env = flyte.TaskEnvironment.from_task(
    "gffread",
    gffread_gff_to_gtf.as_task(),
    gffread_transcripts_fasta.as_task(),
)
