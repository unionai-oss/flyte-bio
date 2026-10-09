"""subread — read summarization with featureCounts.

Exposes :func:`featurecounts`, which counts one BAM's reads over the
features of an annotation and writes ``<prefix>.featureCounts.tsv`` and its
``.summary``.

featureCounts records its command line (input names, ``-T`` threads, ``-o``
output) in the counts header. As upstream's staging does, the inputs are
symlinked into a scratch directory under their original basenames and passed
as relative paths, and the thread count is an explicit input, so the header
matches upstream's for the same inputs.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

# Pinned biocontainer URI. Upstream's main.nf names the bare
# `biocontainers/subread:...` (Docker Hub) — use the quay.io mirror.
SUBREAD_IMAGE = "quay.io/biocontainers/subread:2.0.6--he4a0461_2"

DEFAULT_RESOURCES = flyte.Resources(cpu=4, memory="8Gi")

# Upstream maps the sample's strandedness to `-s` and adds `-p` for paired-end.
featurecounts_cmd = shell.create(
    name="featurecounts",
    image=SUBREAD_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "bam": File,
        "annotation": File,
        "prefix": str,
        "strandedness": str,
        "single_end": bool,
        "threads": int,
        "args": str,
    },
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        BAMS=({inputs.bam}); ANN=({inputs.annotation})
        P={inputs.prefix}
        ARGS={inputs.args}
        case {inputs.strandedness} in
            forward) S=1 ;;
            reverse) S=2 ;;
            *) S=0 ;;
        esac
        PE=()
        [ {inputs.single_end} = true ] || PE=(-p)
        W=$(mktemp -d); cd "$W"
        B=$(basename "${BAMS[0]}"); A=$(basename "${ANN[0]}")
        ln -s "${BAMS[0]}" "$B"; ln -s "${ANN[0]}" "$A"
        featureCounts \
            $ARGS \
            "${PE[@]}" \
            -T {inputs.threads} \
            -a "$A" \
            -s $S \
            -o "$P.featureCounts.tsv" \
            "$B"
        mv "$P.featureCounts.tsv" "$P.featureCounts.tsv.summary" {outputs.results}/
    """,
)


env = flyte.TaskEnvironment.from_task(
    "subread",
    featurecounts_cmd.as_task(),
)


@dataclass
class FeatureCountsResult:
    counts: File  # <prefix>.featureCounts.tsv
    summary: File  # <prefix>.featureCounts.tsv.summary


async def featurecounts(
    bam: File,
    annotation: File,
    prefix: str,
    strandedness: str = "unstranded",
    single_end: bool = False,
    threads: int = 4,
    args: str = "",
) -> FeatureCountsResult:
    """Run featureCounts on ``bam``; ``strandedness`` is forward / reverse / anything else (unstranded)."""
    results = await featurecounts_cmd(
        bam=bam,
        annotation=annotation,
        prefix=prefix,
        strandedness=strandedness,
        single_end=single_end,
        threads=threads,
        args=args,
    )
    counts = await results.get_file(f"{prefix}.featureCounts.tsv")
    summary = await results.get_file(f"{prefix}.featureCounts.tsv.summary")
    if counts is None or summary is None:
        raise FileNotFoundError(f"featureCounts wrote no {prefix}.featureCounts.tsv[.summary]")
    return FeatureCountsResult(counts=counts, summary=summary)
