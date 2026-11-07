"""evigene — EvidentialGene transcript reduction (``tr2aacds``).

Exposes :func:`evigene_tr2aacds`, which reduces an over-assembled
transcriptome (e.g. several assemblers' transcripts pooled) to a
non-redundant set: it translates the transcripts to coding sequences and
proteins, clusters them, and keeps the best main and alternative
transcripts per locus. Outputs keep upstream's names, two directories:

- ``okayset/``: the kept transcripts, ``<prefix>.okay.{mrna,cds,aa}``, the
  culled ``<prefix>.cull.{mrna,cds,aa}``, ``<prefix>.mainalt.tab``,
  ``<prefix>.pubids``, ``<prefix>.ann.txt``, ``<prefix>.genesum.txt``, …;
- ``dropset/``: the dropped ones, ``<prefix>.drop.{tr,cds,aa}``.

tr2aacds names its outputs after the FASTA's name up to its first dot; as
upstream, that part of each output name is replaced with ``prefix`` when they
differ. ``-NCPU`` follows the task's CPUs and ``-MAXMEM`` is 80% of its
memory in MiB, as upstream.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

EVIGENE_IMAGE = "quay.io/biocontainers/evigene:23.7.15--hdfd78af_1"
# Upstream labels it process_medium (6 CPUs, 36 GB); large transcriptomes may need more.
DEFAULT_RESOURCES = flyte.Resources(cpu=4, memory="16Gi")

evigene_tr2aacds_cmd = shell.create(
    name="evigene_tr2aacds",
    image=EVIGENE_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"fasta": File, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"dropset": Dir, "okayset": Dir},
    script=rf"""
        {TASK_LIMITS}
        FA=({{inputs.fasta}})
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        W=$(mktemp -d); cd "$W"
        NAME=$(basename "${{FA[0]}}")
        ln -s "${{FA[0]}}" "$NAME"
        SIMPLE_NAME=${{NAME%%.*}}

        # The container's entrypoint sets EVIGENEHOME and the tool paths; source it, as upstream.
        if [ -f "/usr/local/env-activate.sh" ]; then
            set +u  # Otherwise, errors out because of various unbound variables
            . "/usr/local/env-activate.sh"
        fi

        # tr2aacds dies with exit code 255, which would be retried as a system error.
        $EVIGENEHOME/scripts/prot/tr2aacds.pl \
            $ARGS \
            -NCPU=$TASK_CPUS \
            -MAXMEM=$(( TASK_MEMORY_MB * 8 / 10 )) \
            -cdnaseq "$NAME" || exit 1
        if [ ! -d dropset ] || [ ! -d okayset ]; then
            echo "tr2aacds wrote no dropset/ or okayset/" >&2
            exit 1
        fi

        if [ "$SIMPLE_NAME" != "$P" ]; then
            for D in dropset okayset; do
                find "$D" -type f -exec sh -c 'mv "$1" "$(echo $1 | sed s/'"$SIMPLE_NAME"'/'"$P"'/)"' sh {{}} \;
            done
        fi
        cp -r dropset/. {{outputs.dropset}}/
        cp -r okayset/. {{outputs.okayset}}/
    """,
)

env = flyte.TaskEnvironment.from_task("evigene", evigene_tr2aacds_cmd.as_task())


@dataclass
class EvigeneTr2aacdsResult:
    okayset: Dir  # okayset/: the kept transcripts, CDS and proteins, and tables
    dropset: Dir  # dropset/: <prefix>.drop.{tr,cds,aa}


async def evigene_tr2aacds(fasta: File, prefix: str = "sample", args: str = "") -> EvigeneTr2aacdsResult:
    """Reduce the transcripts in ``fasta`` to a non-redundant set, named ``<prefix>.*``."""
    dropset, okayset = await evigene_tr2aacds_cmd(fasta=fasta, prefix=prefix, args=args)
    return EvigeneTr2aacdsResult(okayset=okayset, dropset=dropset)
