"""salmon — transcript-level RNA-seq quantification.

Currently exposes the index builder:

- :data:`salmon_index` — build a salmon index from a transcript FASTA,
  optionally decoy-aware when a genome FASTA is supplied.

The quantifier (``salmon_quant``) lands in a later batch. All salmon tasks
share one biocontainer.
"""



import flyte
from flyte.extras import shell
from flyte.io import Dir, File

# Pinned biocontainer URI. Upstream's main.nf names the bare
# `biocontainers/salmon:...` (Docker Hub, not anonymously pullable) — use the
# quay.io biocontainers mirror.
SALMON_IMAGE = "quay.io/biocontainers/salmon:1.10.3--h6dccd9a_2"

INDEX_RESOURCES = flyte.Resources(cpu=2, memory="12Gi")


# When a genome FASTA is supplied, build a decoy-aware index: the genome
# sequence names become decoys and the transcript+genome concatenation
# ("gentrome") is indexed. Inputs are already decompressed upstream, so no
# gzip handling here.
salmon_index = shell.create(
    name="salmon_index",
    image=SALMON_IMAGE,
    resources=INDEX_RESOURCES,
    inputs={"transcript_fasta": File, "genome_fasta": File | None},
    outputs={"index": Dir},
    script=r"""
        if [ -e {inputs.genome_fasta} ]; then
            grep '^>' {inputs.genome_fasta} | cut -d ' ' -f 1 | sed 's/>//g' > decoys.txt
            cat {inputs.transcript_fasta} {inputs.genome_fasta} > gentrome.fa
            salmon index --threads 2 -t gentrome.fa -d decoys.txt -i {outputs.index}
        else
            salmon index --threads 2 -t {inputs.transcript_fasta} -i {outputs.index}
        fi
    """,
)


env = flyte.TaskEnvironment.from_task(
    "salmon",
    salmon_index.as_task(),
)
