"""stringtie — transcript assembly and quantification.

Exposes :func:`stringtie`, which assembles transcripts from one BAM and
estimates their abundance, optionally guided by a reference annotation
(``-G``, which also writes coverage and Ballgown tables). Outputs keep
upstream's names: ``<prefix>.transcripts.gtf``,
``<prefix>.gene.abundance.txt``, ``<prefix>.coverage.gtf`` and the
``<prefix>.ballgown/`` directory.

It also exposes :func:`stringtie_merge`, which merges several samples'
transcripts GTFs into one non-redundant ``<prefix>.gtf``, optionally guided
by a reference annotation. Upstream runs it on a newer StringTie (3.0.3), so
it has its own environment, which ``env`` depends on.

StringTie records its command line (input names, ``-p`` threads) in the
transcripts GTF header. As upstream's staging does, the inputs are symlinked
into a scratch directory under their original basenames and passed as
relative paths, and the thread count is an explicit input.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

# Pinned biocontainer URI. Upstream's main.nf names the bare
# `biocontainers/stringtie:...` (Docker Hub) — use the quay.io mirror.
STRINGTIE_IMAGE = "quay.io/biocontainers/stringtie:2.2.3--h43eeafb_0"

DEFAULT_RESOURCES = flyte.Resources(cpu=4, memory="8Gi")

MERGE_IMAGE = "community.wave.seqera.io/library/stringtie:3.0.3--e8043d00caecd051"
# The merged GTF's header records `-p <threads>`; upstream's tests ran it with 4 CPUs.
MERGE_RESOURCES = flyte.Resources(cpu=4, memory="4Gi")

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

# As upstream, the GTFs are passed by their basenames (recorded in the merged GTF's header).
stringtie_merge_cmd = shell.create(
    name="stringtie_merge",
    image=MERGE_IMAGE,
    resources=MERGE_RESOURCES,
    inputs={"gtfs": list[File], "annotation_gtf": File | None, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        GTFS=({{inputs.gtfs}})
        shopt -s nullglob; ANNOTATION=({{inputs.annotation_gtf}}); shopt -u nullglob
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        W=$(mktemp -d); cd "$W"
        NAMES=()
        for f in "${{GTFS[@]}}"; do ln -s "$f" "$(basename "$f")"; NAMES+=("$(basename "$f")"); done
        REF=()
        if [ ${{#ANNOTATION[@]}} -gt 0 ]; then
            G=$(basename "${{ANNOTATION[0]}}"); ln -s "${{ANNOTATION[0]}}" "$G"; REF=(-G "$G")
        fi
        stringtie \
            --merge \
            "${{NAMES[@]}}" \
            "${{REF[@]}}" \
            -o "$P.gtf" \
            -p $TASK_CPUS \
            $ARGS
        mv "$P.gtf" {{outputs.results}}/
    """,
)


env = flyte.TaskEnvironment.from_task(
    "stringtie",
    stringtie_cmd.as_task(),
)
# A TaskEnvironment holds one image; `env` pulls in the merge env too.
merge_env = flyte.TaskEnvironment.from_task("stringtie_merge", stringtie_merge_cmd.as_task())
env.add_dependency(merge_env)


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


async def stringtie_merge(gtfs: list[File], prefix: str, annotation_gtf: File | None = None, args: str = "") -> File:
    """Merge the transcripts ``gtfs`` (in name order) into ``<prefix>.gtf``, guided by ``annotation_gtf``."""
    if not gtfs:
        raise ValueError("stringtie_merge needs at least one GTF")
    names = [f.name for f in gtfs]
    if len(set(names)) != len(names):
        raise ValueError(f"the GTFs to merge must have distinct names: {names}")
    results = await stringtie_merge_cmd(gtfs=gtfs, annotation_gtf=annotation_gtf, prefix=prefix, args=args)
    out = await results.get_file(f"{prefix}.gtf")
    if out is None:
        raise FileNotFoundError(f"stringtie --merge wrote no {prefix}.gtf")
    return out
