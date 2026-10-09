"""stringtie — transcript assembly and quantification.

Exposes :func:`stringtie`, which assembles transcripts from one BAM and
estimates their abundance, optionally guided by a reference annotation
(``-G``, which also writes coverage and Ballgown tables). Outputs keep
upstream's names: ``<prefix>.transcripts.gtf``,
``<prefix>.gene.abundance.txt``, ``<prefix>.coverage.gtf`` and the
``<prefix>.ballgown/`` directory.

StringTie records its command line (input names, ``-p`` threads) in the
transcripts GTF header. As upstream's staging does, the inputs are symlinked
into a scratch directory under their original basenames and passed as
relative paths, and the thread count is an explicit input.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

# Pinned biocontainer URI. Upstream's main.nf names the bare
# `biocontainers/stringtie:...` (Docker Hub) — use the quay.io mirror.
STRINGTIE_IMAGE = "quay.io/biocontainers/stringtie:2.2.3--h43eeafb_0"

DEFAULT_RESOURCES = flyte.Resources(cpu=4, memory="8Gi")

stringtie_cmd = shell.create(
    name="stringtie",
    image=STRINGTIE_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"bam": File, "gtf": File | None, "prefix": str, "strandedness": str, "threads": int, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        BAMS=({inputs.bam})
        shopt -s nullglob; GTFS=({inputs.gtf}); shopt -u nullglob
        P={inputs.prefix}
        ARGS={inputs.args}
        case {inputs.strandedness} in
            forward) STRAND=(--fr) ;;
            reverse) STRAND=(--rf) ;;
            *) STRAND=() ;;
        esac
        W=$(mktemp -d); cd "$W"
        B=$(basename "${BAMS[0]}"); ln -s "${BAMS[0]}" "$B"
        STAGED=("$B")
        REF=()
        if [ ${#GTFS[@]} -gt 0 ]; then
            G=$(basename "${GTFS[0]}"); ln -s "${GTFS[0]}" "$G"; STAGED+=("$G")
            REF=(-G "$G")
            EXTRA=(-C "$P.coverage.gtf" -b "$P.ballgown")
        else
            EXTRA=()
        fi
        stringtie \
            "$B" \
            "${STRAND[@]}" \
            "${REF[@]}" \
            -o "$P.transcripts.gtf" \
            -A "$P.gene.abundance.txt" \
            "${EXTRA[@]}" \
            -p {inputs.threads} \
            $ARGS
        for f in *; do
            keep=1; for s in "${STAGED[@]}"; do [ "$f" = "$s" ] && keep=0; done
            [ $keep = 1 ] && mv "$f" {outputs.results}/
        done
        true  # the loop's last test can be false (a staged input); that's not a failure
    """,
)


env = flyte.TaskEnvironment.from_task(
    "stringtie",
    stringtie_cmd.as_task(),
)


@dataclass
class StringTieResult:
    transcripts_gtf: File  # <prefix>.transcripts.gtf
    abundance: File  # <prefix>.gene.abundance.txt
    coverage_gtf: File | None  # <prefix>.coverage.gtf (with a reference GTF)
    results: Dir  # everything, incl. <prefix>.ballgown/ (with a reference GTF)


async def stringtie(
    bam: File,
    prefix: str,
    gtf: File | None = None,
    strandedness: str = "unstranded",
    threads: int = 4,
    args: str = "",
) -> StringTieResult:
    """Run StringTie on ``bam``; ``strandedness`` is forward / reverse / anything else (unstranded)."""
    results = await stringtie_cmd(
        bam=bam,
        gtf=gtf,
        prefix=prefix,
        strandedness=strandedness,
        threads=threads,
        args=args,
    )
    transcripts = await results.get_file(f"{prefix}.transcripts.gtf")
    abundance = await results.get_file(f"{prefix}.gene.abundance.txt")
    if transcripts is None or abundance is None:
        raise FileNotFoundError(f"stringtie wrote no {prefix}.transcripts.gtf / .gene.abundance.txt")
    coverage = await results.get_file(f"{prefix}.coverage.gtf") if gtf is not None else None
    return StringTieResult(transcripts, abundance, coverage, results)
