"""megahit — MEGAHIT de novo assembly of short reads (metagenomes, single genomes).

Exposes :func:`megahit`, which assembles single-end reads (``-r``) or paired
reads (``-1``/``-2``). Several files per mate (a co-assembly) are joined with
commas in list order, as upstream.

The outputs keep upstream's names, gzipped with ``pigz --no-name``
(``args2`` goes to pigz):

- ``<prefix>.contigs.fa.gz``: the final contigs;
- ``intermediate_contigs/k<k>.{contigs,addi,local,final.contigs}.fa.gz``:
  the contigs of each k-mer iteration;
- ``<prefix>.log``: MEGAHIT's log.

``-t`` and pigz's ``-p`` are the task's CPUs. Unlike upstream, ``-m`` (the
memory MEGAHIT plans its graph for, by default 90% of the machine's) is 90%
of the task's memory; an ``-m`` in ``args`` overrides it.
"""

import os
from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

MEGAHIT_IMAGE = "community.wave.seqera.io/library/megahit_pigz:87a590163e594224"

# Upstream's label asks for 12 CPUs and 72 GB. These defaults fit small genomes
# and the test data; raise them with megahit_cmd.override(resources=...)
# for metagenomes. The 4 CPUs match upstream's tests.
DEFAULT_RESOURCES = flyte.Resources(cpu=4, memory="16Gi")

# A list[File] is staged under its basenames and globbed alphabetically, so the
# list order travels as "/"-joined names ("/" can't occur in a basename).
megahit_cmd = shell.create(
    name="megahit",
    image=MEGAHIT_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "reads_1": list[File],
        "names_1": str,
        "reads_2": list[File],
        "names_2": str,
        "prefix": str,
        "args": str,
        "args2": str,
    },
    defaults={"args": "", "args2": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        shopt -s nullglob
        R1=({{inputs.reads_1}}); R2=({{inputs.reads_2}})
        shopt -u nullglob
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        ARGS2={{inputs.args2}}
        OUT={{outputs.results}}

        # Comma-join a mate's files in the order given by its "/"-joined names.
        joined() {{
            local DIR=$1 NAMES=$2 LIST="" N
            IFS=/ read -ra N <<< "$NAMES"
            for n in "${{N[@]}}"; do LIST="$LIST${{LIST:+,}}$DIR/$n"; done
            echo "$LIST"
        }}
        if [ ${{#R2[@]}} -gt 0 ]; then
            READS=(-1 "$(joined "$(dirname "${{R1[0]}}")" {{inputs.names_1}})" \
                -2 "$(joined "$(dirname "${{R2[0]}}")" {{inputs.names_2}})")
        else
            READS=(-r "$(joined "$(dirname "${{R1[0]}}")" {{inputs.names_1}})")
        fi

        W=$(mktemp -d); cd "$W"
        megahit \
            -m $(( TASK_MEMORY_MB * 1048576 * 9 / 10 )) \
            $ARGS \
            -t $TASK_CPUS \
            "${{READS[@]}}" \
            --out-prefix "$P"

        pigz \
            --no-name \
            -p $TASK_CPUS \
            $ARGS2 \
            megahit_out/*.fa \
            megahit_out/intermediate_contigs/*.fa

        mv "megahit_out/$P.contigs.fa.gz" "megahit_out/$P.log" megahit_out/intermediate_contigs "$OUT"/
    """,
)

env = flyte.TaskEnvironment.from_task("megahit", megahit_cmd.as_task())


@dataclass
class MegahitResult:
    contigs: File  # <prefix>.contigs.fa.gz
    k_contigs: list[File]  # intermediate_contigs/k*.contigs.fa.gz (as upstream's glob, with k*.final.contigs.fa.gz)
    addi_contigs: list[File]  # intermediate_contigs/k*.addi.fa.gz
    local_contigs: list[File]  # intermediate_contigs/k*.local.fa.gz
    kfinal_contigs: list[File]  # intermediate_contigs/k*.final.contigs.fa.gz
    log: File
    results: Dir  # everything above


def file_name(file: File) -> str:
    return file.name or os.path.basename(file.path)


def names_of(files: list[File]) -> str:
    names = [file_name(f) for f in files]
    if len(set(names)) != len(names):
        raise ValueError(f"megahit needs reads with distinct names, got {names}")
    return "/".join(names)


async def megahit(
    reads_1: File | list[File],
    reads_2: File | list[File] | None = None,
    prefix: str = "sample",
    args: str = "",
    args2: str = "",
) -> MegahitResult:
    """Assemble single-end ``reads_1``, or ``reads_1``/``reads_2`` pairs (one or more files per mate)."""
    r1 = reads_1 if isinstance(reads_1, list) else [reads_1]
    r2 = [] if reads_2 is None else reads_2 if isinstance(reads_2, list) else [reads_2]
    if not r1:
        raise ValueError("megahit needs reads")
    if r2 and len(r1) != len(r2):
        raise ValueError(f"megahit needs as many R2 files as R1 files, got {len(r1)} and {len(r2)}")
    results = await megahit_cmd(
        reads_1=r1,
        names_1=names_of(r1),
        reads_2=r2,
        names_2=names_of(r2) if r2 else "",
        prefix=prefix,
        args=args,
        args2=args2,
    )

    files = sorted([f async for f in results.walk()], key=file_name)

    def named(name: str) -> File:
        found = [f for f in files if file_name(f) == name]
        if not found:
            raise FileNotFoundError(f"megahit wrote no {name}")
        return found[0]

    def intermediate(suffix: str) -> list[File]:
        return [
            f
            for f in files
            if "intermediate_contigs" in f.path and file_name(f).startswith("k") and file_name(f).endswith(suffix)
        ]

    return MegahitResult(
        contigs=named(f"{prefix}.contigs.fa.gz"),
        k_contigs=intermediate(".contigs.fa.gz"),
        addi_contigs=intermediate(".addi.fa.gz"),
        local_contigs=intermediate(".local.fa.gz"),
        kfinal_contigs=intermediate(".final.contigs.fa.gz"),
        log=named(f"{prefix}.log"),
        results=results,
    )
