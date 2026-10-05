"""salmon — transcript-level RNA-seq quantification.

Exposes:

- :data:`salmon_index` — build a salmon index from a transcript FASTA,
  optionally decoy-aware when a genome FASTA is supplied.
- :data:`salmon_quant_reads` — quantify FASTQ reads against an index
  (selective-alignment mode).
- :data:`salmon_quant_bam` — quantify a transcriptome-coordinate BAM against
  the transcript FASTA (alignment mode), as the star_salmon path uses.

All salmon tasks share one biocontainer.
"""



import flyte
from flyte.extras import shell
from flyte.io import Dir, File

# Pinned biocontainer URI. Upstream's main.nf names the bare
# `biocontainers/salmon:...` (Docker Hub, not anonymously pullable) — use the
# quay.io biocontainers mirror.
SALMON_IMAGE = "quay.io/biocontainers/salmon:1.10.3--h6dccd9a_2"

INDEX_RESOURCES = flyte.Resources(cpu=2, memory="12Gi")
QUANT_RESOURCES = flyte.Resources(cpu=2, memory="12Gi")


# When a genome FASTA is supplied, build a decoy-aware index: the genome
# sequence names become decoys and the transcript+genome concatenation
# ("gentrome") is indexed. Inputs are already decompressed upstream.
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


# salmon's results dir is nested (aux_info/, libParams/, logs/). We return it
# whole as a Dir — downstream (tximport) consumes the salmon results directory
# directly. (Requires Flyte platform >=2.0.23, which supports nested-directory
# blob uploads; flyteorg/flyte#7490.)
#
# Reads mode: `reads` is a 1-file (single-end) or 2-file (paired) list. An
# empty lib_type means auto-detect ('A').
salmon_quant_reads = shell.create(
    name="salmon_quant_reads",
    image=SALMON_IMAGE,
    resources=QUANT_RESOURCES,
    inputs={"reads": list[File], "index": Dir, "gtf": File, "lib_type": str},
    defaults={"lib_type": ""},
    outputs={"results": Dir},
    script=r"""
        LT={inputs.lib_type}
        [ -z "$LT" ] && LT="A"
        READS=({inputs.reads})
        if [ ${#READS[@]} -eq 1 ]; then
            RR="-r ${READS[0]}"
        else
            RR="-1 ${READS[0]} -2 ${READS[1]}"
        fi
        salmon quant \
            --geneMap {inputs.gtf} \
            --threads 2 \
            --libType=$LT \
            --index {inputs.index} \
            $RR \
            -o {outputs.results}
    """,
)


# Alignment mode (star_salmon): quantify a transcriptome-coordinate BAM
# against the transcript FASTA. No index — the BAM is already aligned.
salmon_quant_bam = shell.create(
    name="salmon_quant_bam",
    image=SALMON_IMAGE,
    resources=QUANT_RESOURCES,
    inputs={"bam": File, "transcript_fasta": File, "gtf": File, "lib_type": str},
    defaults={"lib_type": ""},
    outputs={"results": Dir},
    script=r"""
        LT={inputs.lib_type}
        [ -z "$LT" ] && LT="A"
        salmon quant \
            --geneMap {inputs.gtf} \
            --threads 2 \
            --libType=$LT \
            -t {inputs.transcript_fasta} \
            -a {inputs.bam} \
            -o {outputs.results}
    """,
)


env = flyte.TaskEnvironment.from_task(
    "salmon",
    salmon_index.as_task(),
    salmon_quant_reads.as_task(),
    salmon_quant_bam.as_task(),
)
