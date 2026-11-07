"""rsem — RSEM reference preparation.

Exposes :func:`rsem_preparereference`, which builds an RSEM reference from a
genome FASTA and a GTF. With ``--star`` in ``args`` it first builds a STAR
genome index in the same directory, as upstream (``star_args`` then go to
``STAR --runMode genomeGenerate``).

The reference keeps upstream's layout: a Dir ``rsem`` holding the input FASTA
and ``genome.{chrlist,grp,idx.fa,n2g.idx.fa,seq,ti,transcripts.fa}`` (plus the
STAR index with ``--star``); the transcript FASTA ``genome.transcripts.fa`` is
also returned on its own.
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

env = flyte.TaskEnvironment.from_task("rsem", rsem_preparereference_cmd.as_task())


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
