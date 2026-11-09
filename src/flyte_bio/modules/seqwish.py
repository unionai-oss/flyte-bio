"""seqwish — induce a variation graph from pairwise alignments of sequences.

Exposes :func:`seqwish_induce`, which builds a GFA graph, ``<prefix>.gfa``,
from one or more PAF files of alignments between the sequences of a FASTA
(e.g. from :func:`flyte_bio.modules.wfmash.wfmash`).

As upstream, several PAFs are passed as one comma-separated
``--paf-alns``, and seqwish gets the task's CPUs as threads; options in
``args`` come last.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

SEQWISH_IMAGE = "quay.io/biocontainers/seqwish:0.7.11--h5ca1c30_1"

# Upstream labels seqwish process_medium.
SEQWISH_RESOURCES = flyte.Resources(cpu=4, memory="8Gi")

seqwish_induce_cmd = shell.create(
    name="seqwish_induce",
    image=SEQWISH_IMAGE,
    resources=SEQWISH_RESOURCES,
    inputs={"paf": list[File], "fasta": File, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        set -e
        PAF=({{inputs.paf}})
        FA=({{inputs.fasta}})
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        NAMES=()
        for f in "${{PAF[@]}}"; do ln -s "$f" "$(basename "$f")"; NAMES+=("$(basename "$f")"); done
        ln -s "${{FA[0]}}" "$(basename "${{FA[0]}}")"
        ALNS=$(IFS=,; echo "${{NAMES[*]}}")
        seqwish --threads $TASK_CPUS --paf-alns="$ALNS" --seqs="$(basename "${{FA[0]}}")" --gfa="$P.gfa" $ARGS || exit 1
        mv "$P.gfa" "$OUT/"
    """,
)

env = flyte.TaskEnvironment.from_task("seqwish", seqwish_induce_cmd.as_task())


@dataclass
class SeqwishResult:
    gfa: File  # <prefix>.gfa
    results: Dir


async def seqwish_induce(paf: list[File], fasta: File, prefix: str = "sample", args: str = "") -> SeqwishResult:
    """Induce a variation graph from the alignments in ``paf`` between ``fasta``'s sequences: ``<prefix>.gfa``."""
    if not paf:
        raise ValueError("seqwish_induce: pass at least one PAF")
    results = await seqwish_induce_cmd(paf=paf, fasta=fasta, prefix=prefix, args=args)
    found = await results.get_file(f"{prefix}.gfa")
    if found is None:
        raise FileNotFoundError(f"seqwish wrote no {prefix}.gfa")
    return SeqwishResult(gfa=found, results=results)
