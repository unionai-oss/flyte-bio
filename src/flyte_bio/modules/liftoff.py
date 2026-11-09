"""liftoff — Liftoff, lifting genome annotations over from one assembly to another.

Exposes :func:`liftoff`, which maps the genes of a reference annotation
(``-g``, GFF3 or GTF) or of a gffutils database (``-db``) from the reference
FASTA onto the target FASTA. Either FASTA may be gzipped. Further Liftoff
options (``-polish``, ``-copies``, …) go in ``args``.

The outputs keep upstream's names:

- ``<prefix>.gff3``: the lifted-over annotation;
- ``<prefix>.polished.gff3``: with ``-polish``, the annotation after
  re-aligning exons to restore their reading frames;
- ``<prefix>.unmapped.txt``: the features that could not be mapped.

As upstream, the FASTAs are read as ``target.fasta`` and
``reference.fasta``. The annotation is linked into the working directory,
where Liftoff writes its database. ``-p`` is the task's CPUs.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

LIFTOFF_IMAGE = "quay.io/biocontainers/liftoff:1.6.3--pyhdfd78af_0"
# Upstream's label asks for 12 CPUs and 72 GB; its tests got 4 CPUs.
DEFAULT_RESOURCES = flyte.Resources(cpu=4, memory="16Gi")

liftoff_cmd = shell.create(
    name="liftoff",
    image=LIFTOFF_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "target_fa": File,
        "ref_fa": File,
        "ref_annotation": File | None,
        "ref_db": File | None,
        "prefix": str,
        "args": str,
    },
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        TARGET=({{inputs.target_fa}}); REF=({{inputs.ref_fa}})
        shopt -s nullglob; ANNOTATION=({{inputs.ref_annotation}}); DB=({{inputs.ref_db}}); shopt -u nullglob
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"

        unpack() {{ if [[ "$1" == *.gz ]]; then zcat "$1" > "$2"; else ln -s "$1" "$2"; fi; }}
        unpack "${{TARGET[0]}}" target.fasta
        unpack "${{REF[0]}}" reference.fasta
        FEATURES=()
        if [ ${{#ANNOTATION[@]}} -gt 0 ]; then
            ln -s "${{ANNOTATION[0]}}" "$(basename "${{ANNOTATION[0]}}")"
            FEATURES+=(-g "$(basename "${{ANNOTATION[0]}}")")
        fi
        if [ ${{#DB[@]}} -gt 0 ]; then
            ln -s "${{DB[0]}}" "$(basename "${{DB[0]}}")"
            FEATURES+=(-db "$(basename "${{DB[0]}}")")
        fi

        liftoff \
            "${{FEATURES[@]}}" \
            -p $TASK_CPUS \
            -o "$P.gff3" \
            -u "$P.unmapped.txt" \
            $ARGS \
            target.fasta \
            reference.fasta

        mv "$P.gff3" "$P.unmapped.txt" "$OUT"/
        if [ -f "$P.gff3_polished" ]; then mv "$P.gff3_polished" "$OUT/$P.polished.gff3"; fi
        true
    """,
)

env = flyte.TaskEnvironment.from_task("liftoff", liftoff_cmd.as_task())


@dataclass
class LiftoffResult:
    gff3: File
    polished_gff3: File | None  # only with -polish
    unmapped_txt: File
    results: Dir  # everything above


async def liftoff(
    target_fa: File,
    ref_fa: File,
    ref_annotation: File | None = None,
    ref_db: File | None = None,
    prefix: str = "sample",
    args: str = "",
) -> LiftoffResult:
    """Lift ``ref_annotation`` (or ``ref_db``) over from ``ref_fa`` to ``target_fa``."""
    if ref_annotation is None and ref_db is None:
        raise ValueError("liftoff needs a reference annotation or database")
    results = await liftoff_cmd(
        target_fa=target_fa, ref_fa=ref_fa, ref_annotation=ref_annotation, ref_db=ref_db, prefix=prefix, args=args
    )
    gff3, unmapped = await results.get_file(f"{prefix}.gff3"), await results.get_file(f"{prefix}.unmapped.txt")
    if gff3 is None or unmapped is None:
        raise FileNotFoundError(f"liftoff wrote no {prefix}.gff3 or {prefix}.unmapped.txt")
    return LiftoffResult(
        gff3=gff3,
        polished_gff3=await results.get_file(f"{prefix}.polished.gff3"),
        unmapped_txt=unmapped,
        results=results,
    )
