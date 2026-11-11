"""bwameth — bisulfite read alignment with bwa-meth.

Exposes:

- :func:`bwameth_index` — the bwa-meth index of a FASTA (``bwameth.py index``,
  or ``index-mem2`` for bwa-mem2): a Dir of ``<fasta name>.bwameth.c2t`` (the
  in-silico converted genome) and its BWA index files
  (``.c2t.{amb,ann,bwt,pac,sa}``, or ``.c2t.{0123,amb,ann,bwt.2bit.64,pac}``
  for bwa-mem2). As upstream, the index doesn't keep the FASTA.
- :func:`bwameth_align` — align one sample's reads (``bwameth.py``), piped to
  ``samtools view -bhS``: ``<prefix>.bam``.

As upstream, the index files and the FASTA are linked side by side
(``<index>/<fasta name>``, the reference bwa-meth expects), both tools get the
task's CPUs as threads, and bwa-meth's check that the index is newer than the
FASTA is skipped. ``args`` go to ``bwameth.py`` (e.g. ``-R`` for a read
group), ``args2`` to ``samtools view``.
"""

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

BWAMETH_IMAGE = "community.wave.seqera.io/library/bwameth:0.2.10--a589797f382631bc"
BWAMETH_SAMTOOLS_IMAGE = "community.wave.seqera.io/library/bwameth_samtools:1ff0dd4d9a699bdc"

# Upstream labels both process_high. These defaults fit small genomes; raise
# them with <cmd>.override(resources=...) for larger ones.
DEFAULT_RESOURCES = flyte.Resources(cpu=4, memory="16Gi")

bwameth_index_cmd = shell.create(
    name="bwameth_index",
    image=BWAMETH_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"fasta": File, "use_mem2": bool},
    defaults={"use_mem2": False},
    outputs={"index": Dir},
    script=rf"""
        {TASK_LIMITS}
        FA=({{inputs.fasta}})
        OUT={{outputs.index}}
        COMMAND=index
        if [ {{inputs.use_mem2}} = true ]; then COMMAND=index-mem2; fi
        NAME=$(basename "${{FA[0]}}")
        ln -s "${{FA[0]}}" "$OUT/$NAME"
        cd "$OUT"
        bwameth.py $COMMAND "$NAME"
        rm "$NAME"
    """,
)

bwameth_align_cmd = shell.create(
    name="bwameth_align",
    image=BWAMETH_SAMTOOLS_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "reads_1": File,
        "reads_2": File | None,
        "fasta": File,
        "index": Dir,
        "prefix": str,
        "args": str,
        "args2": str,
    },
    defaults={"args": "", "args2": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        set -o pipefail
        R1=({{inputs.reads_1}})
        shopt -s nullglob; R2=({{inputs.reads_2}}); shopt -u nullglob
        FA=({{inputs.fasta}})
        INDEX_DIR={{inputs.index}}
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        ARGS2={{inputs.args2}}
        W=$(mktemp -d); cd "$W"
        mkdir index
        for F in "$INDEX_DIR"/*; do ln -s "$F" index/; done
        NAME=$(basename "${{FA[0]}}")
        [ -e "index/$NAME" ] || ln -s "${{FA[0]}}" "index/$NAME"
        export BWA_METH_SKIP_TIME_CHECKS=1
        bwameth.py $ARGS -t $TASK_CPUS --reference "index/$NAME" "${{R1[@]}}" "${{R2[@]}}" \
            | samtools view $ARGS2 -@ $TASK_CPUS -bhS -o {{outputs.results}}/"$P.bam" -
    """,
)

# One image per environment: alignment adds samtools, so `env` has two.
bwameth_index_env = flyte.TaskEnvironment.from_task("bwameth_index", bwameth_index_cmd.as_task())
bwameth_align_env = flyte.TaskEnvironment.from_task("bwameth_align", bwameth_align_cmd.as_task())
env = flyte.TaskEnvironment(name="bwameth", depends_on=[bwameth_index_env, bwameth_align_env])


async def bwameth_index(fasta: File, use_mem2: bool = False) -> Dir:
    """The bwa-meth index of ``fasta`` (bwa, or bwa-mem2 with ``use_mem2``), without the FASTA."""
    return await bwameth_index_cmd(fasta=fasta, use_mem2=use_mem2)


async def bwameth_align(
    reads: list[File], fasta: File, index: Dir, prefix: str, args: str = "", args2: str = ""
) -> File:
    """Align one sample's reads (one file, or an R1/R2 pair) to a :func:`bwameth_index` index: ``<prefix>.bam``."""
    if len(reads) not in (1, 2):
        raise ValueError(f"bwameth_align takes one or two read files, got {len(reads)}")
    results = await bwameth_align_cmd(
        reads_1=reads[0],
        reads_2=reads[1] if len(reads) == 2 else None,
        fasta=fasta,
        index=index,
        prefix=prefix,
        args=args,
        args2=args2,
    )
    bam = await results.get_file(f"{prefix}.bam")
    if bam is None:
        raise FileNotFoundError(f"bwameth_align wrote no {prefix}.bam")
    return bam
