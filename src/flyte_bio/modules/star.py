"""STAR — spliced transcript aligner for RNA-seq.

Currently exposes the index builder:

- :data:`star_genome_generate` — build a STAR genome index from a FASTA +
  GTF, ready for alignment.

The per-sample aligner (``star_align``) lands in a later batch. All STAR
tasks share one biocontainer (STAR + samtools + gawk).
"""



import flyte
from flyte.extras import shell
from flyte.io import Dir, File

# Pinned biocontainer URI (STAR + samtools + gawk).
STAR_IMAGE = "community.wave.seqera.io/library/htslib_samtools_star_gawk:ae438e9a604351a4"

# Indexing is the memory-heavy step; STAR is also thread-scalable.
INDEX_RESOURCES = flyte.Resources(cpu=4, memory="16Gi")


# For small genomes STAR needs a reduced --genomeSAindexNbases, computed from
# the total sequence length (same formula upstream uses). samtools writes the
# .fai next to its input, so symlink the FASTA into the writable work dir
# first. gawk is single-quoted so the shell leaves $-fields and {} blocks to
# gawk (no `{x.y}` form, so the task renderer leaves them alone too).
star_genome_generate = shell.create(
    name="star_genome_generate",
    image=STAR_IMAGE,
    resources=INDEX_RESOURCES,
    inputs={"fasta": File, "gtf": File},
    outputs={"index": Dir},
    script=r"""
        ln -s {inputs.fasta} genome.fa
        samtools faidx genome.fa
        NUM_BASES=$(gawk '
            {sum = sum + $2}
            END {
                x = (log(sum) / log(2)) / 2 - 1
                if (x > 14) printf "%.0f", 14; else printf "%.0f", x
            }
        ' genome.fa.fai)
        STAR \
            --runMode genomeGenerate \
            --genomeDir {outputs.index} \
            --genomeFastaFiles genome.fa \
            --sjdbGTFfile {inputs.gtf} \
            --runThreadN 4 \
            --genomeSAindexNbases $NUM_BASES \
            --limitGenomeGenerateRAM 16000000000
    """,
)


env = flyte.TaskEnvironment.from_task(
    "star",
    star_genome_generate.as_task(),
)
