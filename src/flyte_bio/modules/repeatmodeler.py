"""repeatmodeler — de novo repeat family discovery with RepeatModeler.

Exposes:

- :func:`repeatmodeler_builddatabase` — a RepeatModeler (BLAST) database of a
  genome FASTA (``BuildDatabase``): ``<prefix>.tar.gz`` of the files
  ``<prefix>.{nhr,nin,njs,nnd,nni,nog,nsq,translation}``. A genome's database
  is about the genome's size, so it travels as one archive.
- :func:`repeatmodeler` — the repeat families of such a database
  (``RepeatModeler``): their consensus sequences ``<prefix>.fa`` (a library
  for RepeatMasker), their alignments ``<prefix>.stk`` and the run's log
  ``<prefix>.log``. ``-threads`` follows the task's CPUs.

The image (RepeatModeler 2.0.5 with RepeatMasker 4.1.5) carries only
RepeatMasker's bundled libraries: Dfam's curated consensus sequences (CC0) and
RepeatMasker's repeat proteins; no RepBase (which is licensed). RepeatModeler
builds its families from the genome itself.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

REPEATMODELER_IMAGE = "community.wave.seqera.io/library/repeatmasker_repeatmodeler:b718c9eee4ccffcb"
# Upstream: BuildDatabase process_single, RepeatModeler process_medium. These fit
# test-sized genomes, and the 4 CPUs match upstream's tests; a real genome needs
# far more (RepeatModeler runs for days), set with repeatmodeler_cmd.override(resources=...).
BUILDDATABASE_RESOURCES = flyte.Resources(cpu=1, memory="6Gi")
REPEATMODELER_RESOURCES = flyte.Resources(cpu=4, memory="16Gi")

repeatmodeler_builddatabase_cmd = shell.create(
    name="repeatmodeler_builddatabase",
    image=REPEATMODELER_IMAGE,
    resources=BUILDDATABASE_RESOURCES,
    inputs={"fasta": File, "prefix": str},
    outputs={"results": Dir},
    script=r"""
        FA=({inputs.fasta})
        P={inputs.prefix}
        W=$(mktemp -d); cd "$W"
        NAME=$(basename "${FA[0]}")
        ln -s "${FA[0]}" "$NAME"
        BuildDatabase -name "$P" "$NAME" || exit 1
        FILES=()
        for F in "$P".*; do
            if [ "$F" != "$NAME" ]; then FILES+=("$F"); fi
        done
        tar -czf {outputs.results}/"$P.tar.gz" "${FILES[@]}" || exit 1
    """,
)

repeatmodeler_cmd = shell.create(
    name="repeatmodeler",
    image=REPEATMODELER_IMAGE,
    resources=REPEATMODELER_RESOURCES,
    inputs={"db": File, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        DB=({{inputs.db}})
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        tar -xzf "${{DB[0]}}" || exit 1
        # The database's name: its files' common base name.
        NSQ=(*.nsq)
        if [ ! -f "${{NSQ[0]}}" ]; then echo "no BLAST database (*.nsq) in the archive" >&2; exit 1; fi
        NAME=${{NSQ[0]%.nsq}}

        RepeatModeler -database "$NAME" $ARGS -threads $TASK_CPUS || exit 1

        mv "$NAME-families.fa" "$OUT/$P.fa"
        mv "$NAME-families.stk" "$OUT/$P.stk"
        mv "$NAME-rmod.log" "$OUT/$P.log"
    """,
)

env = flyte.TaskEnvironment.from_task(
    "repeatmodeler", repeatmodeler_builddatabase_cmd.as_task(), repeatmodeler_cmd.as_task()
)


async def repeatmodeler_builddatabase(fasta: File, prefix: str = "sample") -> File:
    """A RepeatModeler database of ``fasta``: ``<prefix>.tar.gz`` of ``<prefix>.*``."""
    results = await repeatmodeler_builddatabase_cmd(fasta=fasta, prefix=prefix)
    found = await results.get_file(f"{prefix}.tar.gz")
    if found is None:
        raise FileNotFoundError(f"BuildDatabase wrote no {prefix}.tar.gz")
    return found


@dataclass
class RepeatmodelerResult:
    fasta: File  # <prefix>.fa: the families' consensus sequences
    stk: File  # <prefix>.stk: the families' alignments (Stockholm)
    log: File  # <prefix>.log


async def repeatmodeler(db: File, prefix: str = "sample", args: str = "") -> RepeatmodelerResult:
    """The repeat families of the database archive ``db`` (from :func:`repeatmodeler_builddatabase`)."""
    results = await repeatmodeler_cmd(db=db, prefix=prefix, args=args)

    async def pick(name: str) -> File:
        found = await results.get_file(name)
        if found is None:
            raise FileNotFoundError(f"RepeatModeler wrote no {name}")
        return found

    return RepeatmodelerResult(
        fasta=await pick(f"{prefix}.fa"), stk=await pick(f"{prefix}.stk"), log=await pick(f"{prefix}.log")
    )
