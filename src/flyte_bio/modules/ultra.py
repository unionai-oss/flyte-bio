"""ultra — uLTRA: splice alignment of long transcript reads, guided by a gene annotation.

Exposes:

- :func:`ultra_index` — ``uLTRA index`` of a genome FASTA and its GTF
  (sorted by position, e.g. with :func:`flyte_bio.modules.gnusort.gnu_sort`
  and ``-k1,1 -k4,4n``): a Dir of ``database.db`` (the annotation's
  gffutils database) and the ``*.pickle`` files uLTRA reads at alignment.
- :func:`ultra_align` — ``uLTRA align`` of one sample's reads (uncompressed
  FASTA/FASTQ) against the genome with that index, sorted by
  ``samtools sort`` into ``<prefix>.bam``.

As upstream, align gets the task's CPUs as uLTRA's ``--t`` and samtools'
``--threads``, the index files are linked into the working directory uLTRA
reads its ``--index`` from, and the genome and reads are linked under their
own names. ``args`` go to uLTRA; align's ``args2`` go to ``samtools sort``.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

# A multi-tool container of uLTRA 0.1, gffutils and samtools 1.17.
ULTRA_IMAGE = (
    "quay.io/biocontainers/mulled-v2-4b749ef583d6de806ddbf51c2d235ac8c14763c6"
    ":c2c0cd48e7ed1cf3f365b421c7389d04e6bfa812-0"
)
# Upstream labels index process_low and align process_medium.
INDEX_RESOURCES = flyte.Resources(cpu=1, memory="8Gi")
ALIGN_RESOURCES = flyte.Resources(cpu=4, memory="16Gi")

ultra_index_cmd = shell.create(
    name="ultra_index",
    image=ULTRA_IMAGE,
    resources=INDEX_RESOURCES,
    inputs={"fasta": File, "gtf": File, "args": str},
    defaults={"args": ""},
    outputs={"index": Dir},
    script=r"""
        FA=({inputs.fasta}); GTF=({inputs.gtf})
        ARGS={inputs.args}
        OUT={outputs.index}
        W=$(mktemp -d); cd "$W"
        FA_NAME=$(basename "${FA[0]}"); ln -s "${FA[0]}" "$FA_NAME"
        GTF_NAME=$(basename "${GTF[0]}"); ln -s "${GTF[0]}" "$GTF_NAME"
        mkdir index
        uLTRA index $ARGS "$FA_NAME" "$GTF_NAME" index/ || exit 1
        mv index/*.db index/*.pickle "$OUT"/
    """,
)

ultra_align_cmd = shell.create(
    name="ultra_align",
    image=ULTRA_IMAGE,
    resources=ALIGN_RESOURCES,
    inputs={"reads": File, "genome": File, "index": Dir, "prefix": str, "args": str, "args2": str},
    defaults={"args": "", "args2": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        READS=({{inputs.reads}}); GENOME=({{inputs.genome}})
        INDEX={{inputs.index}}
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        ARGS2={{inputs.args2}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        for f in "$INDEX"/*; do ln -s "$f" "$(basename "$f")"; done
        READS_NAME=$(basename "${{READS[0]}}"); ln -s "${{READS[0]}}" "$READS_NAME"
        GENOME_NAME=$(basename "${{GENOME[0]}}"); ln -s "${{GENOME[0]}}" "$GENOME_NAME"
        uLTRA align --t $TASK_CPUS --prefix "$P" --index ./ $ARGS "$GENOME_NAME" "$READS_NAME" ./ || exit 1
        samtools sort --threads $TASK_CPUS -o "$OUT/$P.bam" -O BAM $ARGS2 "$P.sam" || exit 1
    """,
)

env = flyte.TaskEnvironment.from_task("ultra", ultra_index_cmd.as_task(), ultra_align_cmd.as_task())


@dataclass
class UltraIndexResult:
    database: File  # database.db: the annotation's gffutils database
    index: Dir  # database.db and the *.pickle files


async def ultra_index(fasta: File, gtf: File, args: str = "") -> UltraIndexResult:
    """uLTRA's index of genome ``fasta`` and its position-sorted annotation ``gtf``."""
    index = await ultra_index_cmd(fasta=fasta, gtf=gtf, args=args)
    database = await index.get_file("database.db")
    if database is None:
        raise FileNotFoundError("uLTRA index wrote no database.db")
    return UltraIndexResult(database=database, index=index)


async def ultra_align(
    reads: File, genome: File, index: Dir, prefix: str = "sample", args: str = "", args2: str = ""
) -> File:
    """``reads`` (FASTA/FASTQ) splice-aligned to ``genome`` with uLTRA ``index``, sorted: ``<prefix>.bam``."""
    results = await ultra_align_cmd(reads=reads, genome=genome, index=index, prefix=prefix, args=args, args2=args2)
    bam = await results.get_file(f"{prefix}.bam")
    if bam is None:
        raise FileNotFoundError(f"uLTRA align wrote no {prefix}.bam")
    return bam
