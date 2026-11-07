"""bwa — BWA short-read alignment.

Exposes:

- :func:`bwa_index` — build a BWA index of a FASTA. The index files
  (``<name>.amb``, ``.ann``, ``.bwt``, ``.pac``, ``.sa``, named after the FASTA
  without its extension, as upstream) land in one output Dir.
- :func:`bwa_mem` — align reads (one or two mates) with BWA-MEM, piped into
  ``samtools sort`` (or ``view``) to ``<prefix>.bam``.
- :func:`bwa_aln` — find suffix-array coordinates of reads with BWA-backtrack
  (``bwa aln``): ``<prefix>.sai`` for single-end reads, ``<prefix>.1.sai`` and
  ``<prefix>.2.sai`` for pairs.
- :func:`bwa_samse` / :func:`bwa_sampe` — turn single-end / paired-end
  ``bwa aln`` coordinates into alignments, coordinate-sorted by
  ``samtools sort`` to ``<prefix>.bam``.
"""

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

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
        {TASK_LIMITS}
        shopt -s nullglob; R2=({{inputs.reads_2}}); shopt -u nullglob
        R1=({{inputs.reads_1}})
        ARGS={{inputs.args}}
        ARGS2={{inputs.args2}}
        RG=({{inputs.read_group}})
        RGARG=(); [ -n "${{RG[0]}}" ] && RGARG=(-R "${{RG[0]}}")
        INDEX=$(find -L {{inputs.index}} -name "*.amb" | sed 's/\.amb$//')
        SAMTOOLS=view; [ {{inputs.sort}} = true ] && SAMTOOLS=sort
        bwa mem $ARGS "${{RGARG[@]}}" -t $TASK_CPUS "$INDEX" "${{R1[0]}}" "${{R2[@]}}" \
            | samtools $SAMTOOLS $ARGS2 --threads $TASK_CPUS -o {{outputs.results}}/{{inputs.prefix}}.bam -
    """,
)


# As upstream: one ``bwa aln`` per mate, threads follow the task's CPUs.
bwa_aln_cmd = shell.create(
    name="bwa_aln",
    image=BWA_IMAGE,
    resources=MEM_RESOURCES,
    inputs={"reads_1": File, "reads_2": File | None, "index": Dir, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        shopt -s nullglob; R2=({{inputs.reads_2}}); shopt -u nullglob
        R1=({{inputs.reads_1}})
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        INDEX=$(find -L {{inputs.index}} -name "*.amb" | sed 's/\.amb$//')
        if [ ${{#R2[@]}} -eq 0 ]; then
            bwa aln $ARGS -t $TASK_CPUS -f "$OUT/$P.sai" "$INDEX" "${{R1[0]}}"
        else
            bwa aln $ARGS -t $TASK_CPUS -f "$OUT/$P.1.sai" "$INDEX" "${{R1[0]}}"
            bwa aln $ARGS -t $TASK_CPUS -f "$OUT/$P.2.sai" "$INDEX" "${{R2[0]}}"
        fi
    """,
)

# As upstream: samse/sampe's SAM is coordinate-sorted by samtools (samse with
# one thread fewer than the task's CPUs, as upstream). The read group is its
# own input so it reaches BWA as one argument (``-r '@RG\tID:...'``).
bwa_samse_cmd = shell.create(
    name="bwa_samse",
    image=BWA_IMAGE,
    resources=MEM_RESOURCES,
    inputs={"reads": File, "sai": File, "index": Dir, "prefix": str, "read_group": str, "args": str},
    defaults={"read_group": "", "args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        READS=({{inputs.reads}}); SAI=({{inputs.sai}})
        ARGS={{inputs.args}}
        RG=({{inputs.read_group}})
        RGARG=(); [ -n "${{RG[0]}}" ] && RGARG=(-r "${{RG[0]}}")
        INDEX=$(find -L {{inputs.index}} -name "*.amb" | sed 's/\.amb$//')
        bwa samse $ARGS "${{RGARG[@]}}" "$INDEX" "${{SAI[0]}}" "${{READS[0]}}" \
            | samtools sort -@ $(( TASK_CPUS - 1 )) -O bam - > {{outputs.results}}/{{inputs.prefix}}.bam
    """,
)

bwa_sampe_cmd = shell.create(
    name="bwa_sampe",
    image=BWA_IMAGE,
    resources=MEM_RESOURCES,
    inputs={
        "reads_1": File,
        "reads_2": File,
        "sai_1": File,
        "sai_2": File,
        "index": Dir,
        "prefix": str,
        "read_group": str,
        "args": str,
    },
    defaults={"read_group": "", "args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        R1=({{inputs.reads_1}}); R2=({{inputs.reads_2}})
        S1=({{inputs.sai_1}}); S2=({{inputs.sai_2}})
        ARGS={{inputs.args}}
        RG=({{inputs.read_group}})
        RGARG=(); [ -n "${{RG[0]}}" ] && RGARG=(-r "${{RG[0]}}")
        INDEX=$(find -L {{inputs.index}} -name "*.amb" | sed 's/\.amb$//')
        bwa sampe $ARGS "${{RGARG[@]}}" "$INDEX" "${{S1[0]}}" "${{S2[0]}}" "${{R1[0]}}" "${{R2[0]}}" \
            | samtools sort -@ $TASK_CPUS -O bam - > {{outputs.results}}/{{inputs.prefix}}.bam
    """,
)


env = flyte.TaskEnvironment.from_task(
    "bwa",
    bwa_index_cmd.as_task(),
    bwa_mem_cmd.as_task(),
    bwa_aln_cmd.as_task(),
    bwa_samse_cmd.as_task(),
    bwa_sampe_cmd.as_task(),
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


async def bam_named(results: Dir, prefix: str, tool: str) -> File:
    bam = await results.get_file(f"{prefix}.bam")
    if bam is None:
        raise FileNotFoundError(f"{tool} wrote no {prefix}.bam")
    return bam


async def bwa_aln(reads_1: File, reads_2: File | None, index: Dir, prefix: str, args: str = "") -> list[File]:
    """``bwa aln`` coordinates of the reads: ``[<prefix>.sai]``, or ``[<prefix>.1.sai, <prefix>.2.sai]`` for pairs."""
    results = await bwa_aln_cmd(reads_1=reads_1, reads_2=reads_2, index=index, prefix=prefix, args=args)
    names = [f"{prefix}.sai"] if reads_2 is None else [f"{prefix}.1.sai", f"{prefix}.2.sai"]
    sais = []
    for name in names:
        sai = await results.get_file(name)
        if sai is None:
            raise FileNotFoundError(f"bwa aln wrote no {name}")
        sais.append(sai)
    return sais


async def bwa_samse(reads: File, sai: File, index: Dir, prefix: str, read_group: str = "", args: str = "") -> File:
    """Single-end alignments from ``bwa aln``'s ``sai``, coordinate-sorted to ``<prefix>.bam``.

    ``read_group`` is a BWA ``-r`` value such as ``@RG\\tID:x\\tSM:y``.
    """
    results = await bwa_samse_cmd(reads=reads, sai=sai, index=index, prefix=prefix, read_group=read_group, args=args)
    return await bam_named(results, prefix, "bwa samse")


async def bwa_sampe(
    reads_1: File,
    reads_2: File,
    sai: list[File],
    index: Dir,
    prefix: str,
    read_group: str = "",
    args: str = "",
) -> File:
    """Paired-end alignments from :func:`bwa_aln`'s two ``sai`` files, coordinate-sorted to ``<prefix>.bam``.

    ``read_group`` is a BWA ``-r`` value such as ``@RG\\tID:x\\tSM:y``.
    """
    if len(sai) != 2:
        raise ValueError(f"bwa sampe needs one sai per mate, got {len(sai)}")
    results = await bwa_sampe_cmd(
        reads_1=reads_1,
        reads_2=reads_2,
        sai_1=sai[0],
        sai_2=sai[1],
        index=index,
        prefix=prefix,
        read_group=read_group,
        args=args,
    )
    return await bam_named(results, prefix, "bwa sampe")
