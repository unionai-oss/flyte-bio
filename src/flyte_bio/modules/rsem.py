"""rsem — RSEM reference preparation and expression estimation.

Exposes :func:`rsem_preparereference`, which builds an RSEM reference from a
genome FASTA and a GTF. With ``--star`` in ``args`` it first builds a STAR
genome index in the same directory, as upstream (``star_args`` then go to
``STAR --runMode genomeGenerate``).

The reference keeps upstream's layout: a Dir ``rsem`` holding the input FASTA
and ``genome.{chrlist,grp,idx.fa,n2g.idx.fa,seq,ti,transcripts.fa}`` (plus the
STAR index with ``--star``); the transcript FASTA ``genome.transcripts.fa`` is
also returned on its own.

:func:`rsem_calculateexpression` estimates one sample's gene and isoform
expression against that reference, from FASTQ reads (with ``--star`` in
``args``, aligned by STAR first) or from a transcriptome BAM (``--alignments``,
chosen as upstream when the reads end in ``.bam``). Outputs keep upstream's
names: ``<prefix>.genes.results``, ``<prefix>.isoforms.results``, the model
directory ``<prefix>.stat`` and, when written, ``<prefix>.log``,
``<prefix>.STAR.genome.bam``, ``<prefix>.genome.bam`` and
``<prefix>.transcript.bam``.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

RSEM_IMAGE = "community.wave.seqera.io/library/rsem_star:5acb4e8c03239c32"
# Upstream labels it process_high. Without --star, RSEM needs little; a human
# STAR index needs about 32 GiB, more than the default.
DEFAULT_RESOURCES = flyte.Resources(cpu=4, memory="16Gi")

# rsem-prepare-reference is Perl, whose die exits 255, which the cluster takes
# for an interruption and retries without end; any failure exits 1 instead. As
# upstream, STAR may use the task's memory less 100 MB.
rsem_preparereference_cmd = shell.create(
    name="rsem_preparereference",
    image=RSEM_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"fasta": File, "gtf": File, "args": str, "star_args": str},
    defaults={"args": "", "star_args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        FA=({{inputs.fasta}}); GTF=({{inputs.gtf}})
        ARGS={{inputs.args}}
        STAR_ARGS={{inputs.star_args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        mkdir rsem
        cp "${{FA[0]}}" rsem/
        FASTA=rsem/$(basename "${{FA[0]}}")

        # As upstream: `--star` selects STAR mode, and every option naming
        # --star (--star, --star-path ...) is then dropped from RSEM's options.
        STAR_MODE=false
        for A in $ARGS; do if [ "$A" = --star ]; then STAR_MODE=true; fi; done
        RSEM_ARGS=()
        for A in $ARGS; do
            if [ $STAR_MODE = true ] && [[ "$A" == *--star* ]]; then continue; fi
            RSEM_ARGS+=("$A")
        done
        if [ $STAR_MODE = true ]; then
            STAR \
                --runMode genomeGenerate \
                --genomeDir rsem/ \
                --genomeFastaFiles "$FASTA" \
                --sjdbGTFfile "${{GTF[0]}}" \
                --runThreadN $TASK_CPUS \
                --limitGenomeGenerateRAM $(( TASK_MEMORY_MB * 1048576 - 100000000 )) \
                $STAR_ARGS \
                || exit 1
        fi
        rsem-prepare-reference \
            --gtf "${{GTF[0]}}" \
            --num-threads $TASK_CPUS \
            "${{RSEM_ARGS[@]}}" \
            "$FASTA" \
            rsem/genome \
            || exit 1
        cp rsem/genome.transcripts.fa "$OUT"/
        mv rsem "$OUT"/rsem
    """,
)

# As upstream: `--paired-end` follows `paired` when given; unknown (""), it is
# detected from the BAM's flags or from the number of FASTQ files. Temporary
# files go to ./tmp/ unless `args` names a --temporary-folder.
rsem_calculateexpression_cmd = shell.create(
    name="rsem_calculateexpression",
    image=RSEM_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "reads_1": File,
        "reads_2": File | None,
        "index": Dir,
        "strandedness": str,
        "paired": str,
        "prefix": str,
        "args": str,
    },
    defaults={"strandedness": "", "paired": "", "args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        R1=({{inputs.reads_1}})
        shopt -s nullglob; R2=({{inputs.reads_2}}); shopt -u nullglob
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        STRANDEDNESS={{inputs.strandedness}}
        PAIRED={{inputs.paired}}
        OUT={{outputs.results}}
        INDEX=$(find -L {{inputs.index}}/ -name "*.grp" | sed 's/\.grp$//')
        W=$(mktemp -d); cd "$W"
        READS=("$(basename "${{R1[0]}}")"); ln -s "${{R1[0]}}" "${{READS[0]}}"
        if [ ${{#R2[@]}} -gt 0 ]; then READS+=("$(basename "${{R2[0]}}")"); ln -s "${{R2[0]}}" "${{READS[1]}}"; fi

        if [[ " $ARGS " != *--temporary-folder* ]]; then ARGS="$ARGS --temporary-folder ./tmp/"; fi
        STRAND=()
        if [ "$STRANDEDNESS" = forward ] || [ "$STRANDEDNESS" = reverse ]; then
            STRAND=(--strandedness "$STRANDEDNESS")
        fi
        IS_BAM=false
        if [[ "${{READS[0],,}}" == *.bam ]]; then IS_BAM=true; fi
        ALIGNMENTS=()
        if [ $IS_BAM = true ]; then ALIGNMENTS=(--alignments); fi
        PAIRED_END=()
        if [ "$PAIRED" = true ]; then
            PAIRED_END=(--paired-end)
        elif [ -z "$PAIRED" ]; then
            if [ $IS_BAM = true ]; then
                if samtools flagstat "${{READS[@]}}" | grep -q 'paired in sequencing'; then
                    PAIRED_END=(--paired-end)
                fi
            elif [ ${{#READS[@]}} -gt 1 ]; then
                PAIRED_END=(--paired-end)
            fi
        fi

        rsem-calculate-expression \
            --num-threads $TASK_CPUS \
            "${{ALIGNMENTS[@]}}" \
            "${{STRAND[@]}}" \
            "${{PAIRED_END[@]}}" \
            $ARGS \
            "${{READS[@]}}" \
            $INDEX \
            "$P" \
            || exit 1
        # The input reads are links here: move only the BAMs RSEM wrote.
        shopt -s nullglob
        mv "$P".genes.results "$P".isoforms.results "$P".stat ./*.log "$OUT"/
        for F in ./*.STAR.genome.bam "$P".genome.bam "$P".transcript.bam; do
            if [ -f "$F" ] && [ ! -L "$F" ]; then mv "$F" "$OUT"/; fi
        done
    """,
)

env = flyte.TaskEnvironment.from_task(
    "rsem", rsem_preparereference_cmd.as_task(), rsem_calculateexpression_cmd.as_task()
)


@dataclass
class RsemReference:
    index: Dir  # the reference directory, `rsem`
    transcript_fasta: File  # genome.transcripts.fa


async def rsem_preparereference(fasta: File, gtf: File, args: str = "", star_args: str = "") -> RsemReference:
    """An RSEM reference (and with ``--star`` in ``args``, a STAR index) for ``fasta`` and ``gtf``."""
    results = await rsem_preparereference_cmd(fasta=fasta, gtf=gtf, args=args, star_args=star_args)
    transcript_fasta = await results.get_file("genome.transcripts.fa")
    if transcript_fasta is None:
        raise FileNotFoundError("rsem-prepare-reference wrote no genome.transcripts.fa")
    return RsemReference(
        index=Dir.from_existing_remote(results.path.rstrip("/") + "/rsem"), transcript_fasta=transcript_fasta
    )


@dataclass
class RsemExpression:
    counts_gene: File  # <prefix>.genes.results
    counts_transcript: File  # <prefix>.isoforms.results
    stat: Dir  # <prefix>.stat: <prefix>.cnt, <prefix>.model, <prefix>.theta
    log: File | None  # <prefix>.log
    bam_star: File | None  # <prefix>.STAR.genome.bam (STAR mode)
    bam_genome: File | None  # <prefix>.genome.bam (with --output-genome-bam)
    bam_transcript: File | None  # <prefix>.transcript.bam
    results: Dir  # everything above


async def rsem_calculateexpression(
    reads_1: File,
    index: Dir,
    reads_2: File | None = None,
    strandedness: str = "",
    paired_end: bool | None = None,
    prefix: str = "sample",
    args: str = "",
) -> RsemExpression:
    """Estimate one sample's expression against an RSEM reference ``index``.

    ``reads_1`` is a FASTQ (``reads_2`` its mate) or a transcriptome BAM.
    ``strandedness`` is ``forward``, ``reverse`` or anything else for
    unstranded; ``paired_end`` None detects it as upstream.
    """
    results = await rsem_calculateexpression_cmd(
        reads_1=reads_1,
        reads_2=reads_2,
        index=index,
        strandedness=strandedness,
        paired="" if paired_end is None else str(paired_end).lower(),
        prefix=prefix,
        args=args,
    )

    async def pick(name: str) -> File:
        found = await results.get_file(name)
        if found is None:
            raise FileNotFoundError(f"rsem-calculate-expression wrote no {name}")
        return found

    return RsemExpression(
        counts_gene=await pick(f"{prefix}.genes.results"),
        counts_transcript=await pick(f"{prefix}.isoforms.results"),
        stat=Dir.from_existing_remote(results.path.rstrip("/") + f"/{prefix}.stat"),
        log=await results.get_file(f"{prefix}.log"),
        bam_star=await results.get_file(f"{prefix}.STAR.genome.bam"),
        bam_genome=await results.get_file(f"{prefix}.genome.bam"),
        bam_transcript=await results.get_file(f"{prefix}.transcript.bam"),
        results=results,
    )
