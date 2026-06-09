"""STAR — spliced transcript aligner for RNA-seq.

Exposes:

- :data:`star_genome_generate` — build a STAR genome index from a FASTA +
  GTF.
- :data:`star_align` — align reads against a STAR index.

STAR emits a dozen prefix-named files whose exact names depend on the
options passed, so :data:`star_align` captures its whole output directory as
one :class:`Dir`; downstream steps pick out the genome BAM, transcriptome
BAM, ``SJ.out.tab``, logs, etc. by name. All STAR tasks share one
biocontainer (STAR + samtools + gawk).
"""



import flyte
from flyte.extras import shell
from flyte.io import Dir, File

# Pinned biocontainer URI (STAR + samtools + gawk).
STAR_IMAGE = "community.wave.seqera.io/library/htslib_samtools_star_gawk:ae438e9a604351a4"

INDEX_RESOURCES = flyte.Resources(cpu=4, memory="16Gi")
ALIGN_RESOURCES = flyte.Resources(cpu=4, memory="16Gi")


# For small genomes STAR needs a reduced --genomeSAindexNbases. Following
# upstream: if the caller already set it via `args`, use args as-is;
# otherwise compute it from the total sequence length (samtools writes the
# .fai next to its input, so symlink the FASTA into the writable work dir
# first). gawk is single-quoted so the shell leaves $-fields and {} blocks to
# gawk; none are `{x.y}`, so the task renderer leaves them alone too.
star_genome_generate = shell.create(
    name="star_genome_generate",
    image=STAR_IMAGE,
    resources=INDEX_RESOURCES,
    inputs={"fasta": File, "gtf": File, "args": str},
    defaults={"args": ""},
    outputs={"index": Dir},
    script=r"""
        ln -s {inputs.fasta} genome.fa
        ARGS={inputs.args}
        if [[ "$ARGS" == *"--genomeSAindexNbases"* ]]; then
            STAR \
                --runMode genomeGenerate \
                --genomeDir {outputs.index} \
                --genomeFastaFiles genome.fa \
                --sjdbGTFfile {inputs.gtf} \
                --runThreadN 4 \
                --limitGenomeGenerateRAM 16000000000 \
                $ARGS
        else
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
                --limitGenomeGenerateRAM 16000000000 \
                $ARGS
        fi
    """,
)


# `reads` is a 1-file (single-end) or 2-file (paired) list; the glob expands
# to "R1" or "R1 R2", which STAR reads as the mate(s). STAR-specific options
# (output type, wig, transcriptome quant, readFilesCommand for gzipped input,
# ...) ride in via `args`. Everything STAR writes lands in the output Dir.
star_align = shell.create(
    name="star_align",
    image=STAR_IMAGE,
    resources=ALIGN_RESOURCES,
    inputs={"reads": list[File], "index": Dir, "gtf": File, "args": str},
    defaults={"args": ""},
    outputs={"output": Dir},
    script=r"""
        STAR \
            --genomeDir {inputs.index} \
            --readFilesIn {inputs.reads} \
            --sjdbGTFfile {inputs.gtf} \
            --runThreadN 4 \
            --outFileNamePrefix {outputs.output}/ \
            {inputs.args}
    """,
)


env = flyte.TaskEnvironment.from_task(
    "star",
    star_genome_generate.as_task(),
    star_align.as_task(),
)
