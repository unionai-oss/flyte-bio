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
# `list[File]` (0 or 1 item) rather than `File | None` until flyteorg/flyte#8118
# is deployed: copilot stages a set optional File as a bare path, not the
# per-input dir the shell glob expects, so it was silently ignored.
salmon_index = shell.create(
    name="salmon_index",
    image=SALMON_IMAGE,
    resources=INDEX_RESOURCES,
    inputs={"transcript_fasta": File, "genome_fasta": list[File]},
    defaults={"genome_fasta": []},
    outputs={"index": Dir},
    script=r"""
        shopt -s nullglob; GENOME=({inputs.genome_fasta}); shopt -u nullglob
        if [ ${#GENOME[@]} -gt 0 ]; then
            grep '^>' "${GENOME[0]}" | cut -d ' ' -f 1 | sed 's/>//g' > decoys.txt
            cat {inputs.transcript_fasta} "${GENOME[0]}" > gentrome.fa
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
# Reads mode: mates are separate inputs (`reads_2` empty for single-end) so
# R1/R2 order is explicit — a list[File] is staged under original basenames
# and globbed alphabetically, which can swap the mates. An empty lib_type
# means auto-detect ('A').
# `list[File]` (0 or 1 item) rather than `File | None` until flyteorg/flyte#8118
# is deployed: copilot stages a set optional File as a bare path, not the
# per-input dir the shell glob expects, so it was silently ignored.
salmon_quant_reads = shell.create(
    name="salmon_quant_reads",
    image=SALMON_IMAGE,
    resources=QUANT_RESOURCES,
    inputs={"reads_1": File, "reads_2": list[File], "index": Dir, "gtf": File, "lib_type": str},
    defaults={"reads_2": [], "lib_type": ""},
    outputs={"results": Dir},
    script=r"""
        LT={inputs.lib_type}
        [ -z "$LT" ] && LT="A"
        shopt -s nullglob; R2=({inputs.reads_2}); shopt -u nullglob
        if [ ${#R2[@]} -eq 0 ]; then
            RR="-r {inputs.reads_1}"
        else
            RR="-1 {inputs.reads_1} -2 ${R2[0]}"
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
