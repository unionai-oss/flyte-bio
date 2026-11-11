"""tiara — deep-learning classification of assembly contigs by domain and organelle.

Exposes :func:`tiara`, which classifies each sequence of a FASTA (may be
gzipped) as archaea, bacteria, eukarya, organelle (mitochondrion, plastid) or
unknown. Outputs keep upstream's names:

- ``<prefix>.txt``: one classification per sequence;
- ``log_<prefix>.txt``: tiara's run log;
- with ``--tf``/``--to-fasta`` in ``args`` (e.g. ``--tf bac``): the sequences
  of each chosen class, renamed as upstream to ``<class>_<prefix>.fasta``
  (``.fasta.gz`` when ``args`` asks tiara to gzip them).

``--threads`` follows the task's CPUs.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

TIARA_IMAGE = "community.wave.seqera.io/library/tiara:1.0.3--e367cb0eaef39fa3"
# Upstream's process_medium.
DEFAULT_RESOURCES = flyte.Resources(cpu=4, memory="16Gi")

tiara_cmd = shell.create(
    name="tiara",
    image=TIARA_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"fasta": File, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        FA=({{inputs.fasta}})
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        FASTA=$(basename "${{FA[0]}}")
        ln -s "${{FA[0]}}" "$FASTA"

        tiara -i "$FASTA" -o "$P.txt" --threads $TASK_CPUS $ARGS || exit 1

        # As upstream: name the per-class FASTAs <class>_<prefix>.fasta[.gz].
        if echo "$ARGS" | grep -qE "tf|to-fasta"; then
            SUFFIX=fasta
            if echo "$ARGS" | grep -q "gz"; then SUFFIX=fasta.gz; fi
            find . -name "*_$FASTA*" | while read -r F; do
                F=$(basename "$F")
                mv "$F" "${{F%%_*}}_$P.$SUFFIX"
            done
        fi
        shopt -s nullglob
        mv "$P".txt* log_*.txt* ./*.fasta ./*.fasta.gz "$OUT"/
    """,
)

env = flyte.TaskEnvironment.from_task("tiara", tiara_cmd.as_task())


@dataclass
class TiaraResult:
    classifications: File  # <prefix>.txt
    log: File  # log_<prefix>.txt
    fasta: list[File]  # <class>_<prefix>.fasta[.gz], with --tf in args
    results: Dir  # everything tiara wrote


def file_name(file: File) -> str:
    return file.path.rstrip("/").rsplit("/", 1)[-1]


async def tiara(fasta: File, prefix: str = "sample", args: str = "") -> TiaraResult:
    """Classify the sequences of ``fasta`` by domain and organelle."""
    results = await tiara_cmd(fasta=fasta, prefix=prefix, args=args)
    files = {file_name(f): f for f in await results.list_files()}
    for name in (f"{prefix}.txt", f"log_{prefix}.txt"):
        if name not in files:
            raise FileNotFoundError(f"tiara wrote no {name}")
    return TiaraResult(
        classifications=files[f"{prefix}.txt"],
        log=files[f"log_{prefix}.txt"],
        fasta=[files[n] for n in sorted(files) if n.endswith((".fasta", ".fasta.gz"))],
        results=results,
    )
