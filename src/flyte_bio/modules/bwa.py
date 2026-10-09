"""bwa — BWA short-read alignment.

Exposes:

- :func:`bwa_index` — build a BWA index of a FASTA. The index files
  (``<name>.amb``, ``.ann``, ``.bwt``, ``.pac``, ``.sa``, named after the FASTA
  without its extension, as upstream) land in one output Dir.
- :func:`bwa_mem` — align reads (one or two mates) with BWA-MEM, piped into
  ``samtools sort`` (or ``view``) to ``<prefix>.bam``.
"""

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

BWA_IMAGE = "community.wave.seqera.io/library/bwa_htslib_samtools:83b50ff84ead50d0"

# Indexing needs ~5.4x the genome size in memory (BWA's own figure): ~17 GB
# for a human genome.
INDEX_RESOURCES = flyte.Resources(cpu=1, memory="24Gi")
BWA_MEM_CPUS = 4
MEM_RESOURCES = flyte.Resources(cpu=BWA_MEM_CPUS, memory="16Gi")

bwa_index_cmd = shell.create(
    name="bwa_index",
    image=BWA_IMAGE,
    resources=INDEX_RESOURCES,
    inputs={"fasta": File, "args": str},
    defaults={"args": ""},
    outputs={"index": Dir},
    script=r"""
        FA=({inputs.fasta})
        NAME=$(basename "${FA[0]}")
        ARGS={inputs.args}
        bwa index $ARGS -p {outputs.index}/"${NAME%.*}" "${FA[0]}"
    """,
)


# As upstream: the index is found by its .amb file, BWA-MEM's SAM is piped
# into samtools (sort, or view without `sort`), and threads follow the task's
# CPUs. The read group is its own input so it reaches BWA as one argument
# (``@RG\tID:...``, BWA expands the \t).
bwa_mem_cmd = shell.create(
    name="bwa_mem",
    image=BWA_IMAGE,
    resources=MEM_RESOURCES,
    inputs={
        "reads_1": File,
        "reads_2": File | None,
        "index": Dir,
        "prefix": str,
        "read_group": str,
        "sort": bool,
        "args": str,
        "args2": str,
    },
    defaults={"read_group": "", "sort": True, "args": "", "args2": ""},
    outputs={"results": Dir},
    script=rf"""
        shopt -s nullglob; R2=({{inputs.reads_2}}); shopt -u nullglob
        R1=({{inputs.reads_1}})
        ARGS={{inputs.args}}
        ARGS2={{inputs.args2}}
        RG=({{inputs.read_group}})
        RGARG=(); [ -n "${{RG[0]}}" ] && RGARG=(-R "${{RG[0]}}")
        INDEX=$(find -L {{inputs.index}} -name "*.amb" | sed 's/\.amb$//')
        SAMTOOLS=view; [ {{inputs.sort}} = true ] && SAMTOOLS=sort
        bwa mem $ARGS "${{RGARG[@]}}" -t {BWA_MEM_CPUS} "$INDEX" "${{R1[0]}}" "${{R2[@]}}" \
            | samtools $SAMTOOLS $ARGS2 --threads {BWA_MEM_CPUS} -o {{outputs.results}}/{{inputs.prefix}}.bam -
    """,
)


env = flyte.TaskEnvironment.from_task(
    "bwa",
    bwa_index_cmd.as_task(),
    bwa_mem_cmd.as_task(),
)


async def bwa_index(fasta: File, args: str = "") -> Dir:
    """BWA index of ``fasta``: a Dir of ``<fasta name>.{amb,ann,bwt,pac,sa}``."""
    return await bwa_index_cmd(fasta=fasta, args=args)


async def bwa_mem(
    reads_1: File,
    reads_2: File | None,
    index: Dir,
    prefix: str,
    read_group: str = "",
    sort: bool = True,
    args: str = "",
    args2: str = "",
) -> File:
    """Align reads with BWA-MEM to ``<prefix>.bam`` (coordinate-sorted unless ``sort=False``).

    ``read_group`` is a BWA ``-R`` value such as ``@RG\\tID:x\\tSM:y``;
    ``args`` go to BWA, ``args2`` to samtools.
    """
    results = await bwa_mem_cmd(
        reads_1=reads_1,
        reads_2=reads_2,
        index=index,
        prefix=prefix,
        read_group=read_group,
        sort=sort,
        args=args,
        args2=args2,
    )
    bam = await results.get_file(f"{prefix}.bam")
    if bam is None:
        raise FileNotFoundError(f"bwa mem wrote no {prefix}.bam")
    return bam
