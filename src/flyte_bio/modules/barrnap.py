"""barrnap — Barrnap, prediction of ribosomal RNA genes.

Exposes :func:`barrnap`, which predicts the rRNA genes (5S, 16S, 23S, …) of a
FASTA (plain or gzipped) for one kingdom: ``bac`` (the default), ``arc``,
``euk`` or ``mito``. The output keeps upstream's name,
``<prefix>_<kingdom>.gff``.

As upstream, a gzipped FASTA is decompressed first (Barrnap reads plain
FASTA), and ``--threads`` follows the task's CPUs.
"""

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

BARRNAP_IMAGE = "community.wave.seqera.io/library/barrnap:0.9--847dcfa21d3169b7"
# Upstream: process_single.
DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory="4Gi")

barrnap_cmd = shell.create(
    name="barrnap",
    image=BARRNAP_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"fasta": File, "kingdom": str, "prefix": str, "args": str},
    defaults={"kingdom": "bac", "args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        FA=({{inputs.fasta}})
        DB={{inputs.kingdom}}
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        W=$(mktemp -d); cd "$W"
        NAME=$(basename "${{FA[0]}}")
        if [[ "$NAME" == *.gz ]]; then
            INPUT=${{NAME%.gz}}
            gunzip -c "${{FA[0]}}" > "$INPUT"
        else
            INPUT=$NAME
            ln -s "${{FA[0]}}" "$INPUT"
        fi

        barrnap \
            $ARGS \
            --threads $TASK_CPUS \
            --kingdom "$DB" \
            "$INPUT" \
            > "${{P}}_$DB.gff" || exit 1
        mv "${{P}}_$DB.gff" {{outputs.results}}/
    """,
)

env = flyte.TaskEnvironment.from_task("barrnap", barrnap_cmd.as_task())


async def barrnap(fasta: File, kingdom: str = "bac", prefix: str = "sample", args: str = "") -> File:
    """The rRNA genes of ``fasta`` for ``kingdom``: ``<prefix>_<kingdom>.gff``."""
    kingdom = kingdom or "bac"
    results = await barrnap_cmd(fasta=fasta, kingdom=kingdom, prefix=prefix, args=args)
    name = f"{prefix}_{kingdom}.gff"
    gff = await results.get_file(name)
    if gff is None:
        raise FileNotFoundError(f"barrnap wrote no {name}")
    return gff
