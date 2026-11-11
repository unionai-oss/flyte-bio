"""bwamem2 — BWA-MEM2 short-read alignment.

Exposes:

- :func:`bwamem2_index` — build a BWA-MEM2 index of a FASTA (plain or
  gzipped). The index files (``<name>.0123``, ``.amb``, ``.ann``,
  ``.bwt.2bit.64``, ``.pac``) are named after the whole FASTA file name, as
  upstream (``genome.fasta.0123``), and land in one output Dir.
- :func:`bwamem2_mem` — align reads (one or two mates) with BWA-MEM2, piped
  into ``samtools sort`` (or ``view``) to ``<prefix>.bam``. As upstream, an
  ``--output-fmt``/``-O`` in ``args2`` picks the format and extension
  (``sam``, ``cram``); CRAM needs the reference ``fasta``.
"""

import re

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

BWAMEM2_IMAGE = "community.wave.seqera.io/library/bwa-mem2_htslib_samtools:db98f81f55b64113"

# Indexing is single-threaded and needs ~28 bytes of memory per reference
# base (upstream sizes it at 280 MB per 10 MB of FASTA, x4 when gzipped).
# The default covers genomes up to ~1 Gbp; a human genome needs ~96 GiB:
# bwamem2_index_cmd.override(resources=flyte.Resources(cpu=1, memory="96Gi")).
INDEX_RESOURCES = flyte.Resources(cpu=1, memory="32Gi")
# Upstream labels this process_high; the human index takes ~16 GB in memory.
MEM_RESOURCES = flyte.Resources(cpu=8, memory="24Gi")

bwamem2_index_cmd = shell.create(
    name="bwamem2_index",
    image=BWAMEM2_IMAGE,
    resources=INDEX_RESOURCES,
    inputs={"fasta": File, "args": str},
    defaults={"args": ""},
    outputs={"index": Dir},
    script=r"""
        FA=({inputs.fasta})
        NAME=$(basename "${FA[0]}")
        ARGS={inputs.args}
        bwa-mem2 index $ARGS -p {outputs.index}/"$NAME" "${FA[0]}"
    """,
)


# As upstream: the index is found by its .amb file, BWA-MEM2's SAM is piped
# into samtools (sort, or view without `sort`), and threads follow the task's
# CPUs. The read group is its own input so it reaches BWA-MEM2 as one argument
# (``@RG\tID:...``, BWA-MEM2 expands the \t).
bwamem2_mem_cmd = shell.create(
    name="bwamem2_mem",
    image=BWAMEM2_IMAGE,
    resources=MEM_RESOURCES,
    inputs={
        "reads_1": File,
        "reads_2": File | None,
        "index": Dir,
        "fasta": File | None,
        "prefix": str,
        "extension": str,
        "read_group": str,
        "sort": bool,
        "args": str,
        "args2": str,
    },
    defaults={"extension": "bam", "read_group": "", "sort": True, "args": "", "args2": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        shopt -s nullglob; R2=({{inputs.reads_2}}); FA=({{inputs.fasta}}); shopt -u nullglob
        R1=({{inputs.reads_1}})
        ARGS={{inputs.args}}
        ARGS2={{inputs.args2}}
        EXT={{inputs.extension}}
        RG=({{inputs.read_group}})
        RGARG=(); [ -n "${{RG[0]}}" ] && RGARG=(-R "${{RG[0]}}")
        REF=(); [ "$EXT" = cram ] && REF=(--reference "${{FA[0]}}")
        INDEX=$(find -L {{inputs.index}} -name "*.amb" | sed 's/\.amb$//')
        SAMTOOLS=view; [ {{inputs.sort}} = true ] && SAMTOOLS=sort
        bwa-mem2 mem $ARGS "${{RGARG[@]}}" -t $TASK_CPUS "$INDEX" "${{R1[0]}}" "${{R2[@]}}" \
            | samtools $SAMTOOLS $ARGS2 -@ $TASK_CPUS "${{REF[@]}}" -o {{outputs.results}}/{{inputs.prefix}}.$EXT -
    """,
)


env = flyte.TaskEnvironment.from_task(
    "bwamem2",
    bwamem2_index_cmd.as_task(),
    bwamem2_mem_cmd.as_task(),
)


async def bwamem2_index(fasta: File, args: str = "") -> Dir:
    """BWA-MEM2 index of ``fasta``: a Dir of ``<fasta file name>.{0123,amb,ann,bwt.2bit.64,pac}``."""
    return await bwamem2_index_cmd(fasta=fasta, args=args)


def output_extension(args2: str) -> str:
    """The alignment's extension, as upstream: the ``--output-fmt``/``-O`` value in ``args2``, else ``bam``."""
    match = re.search(r"(--output-fmt|-O)+\s+(\S+)", args2)
    return match.group(2).lower() if match else "bam"


async def bwamem2_mem(
    reads_1: File,
    reads_2: File | None,
    index: Dir,
    prefix: str,
    fasta: File | None = None,
    read_group: str = "",
    sort: bool = True,
    args: str = "",
    args2: str = "",
) -> File:
    """Align reads with BWA-MEM2 to ``<prefix>.bam`` (coordinate-sorted unless ``sort=False``).

    ``read_group`` is a BWA-MEM2 ``-R`` value such as ``@RG\\tID:x\\tSM:y``;
    ``args`` go to BWA-MEM2, ``args2`` to samtools. ``fasta`` is only needed
    (and required) for CRAM output.
    """
    extension = output_extension(args2)
    if extension == "cram" and fasta is None:
        raise ValueError("bwamem2_mem: a fasta reference is required for CRAM output")
    results = await bwamem2_mem_cmd(
        reads_1=reads_1,
        reads_2=reads_2,
        index=index,
        fasta=fasta,
        prefix=prefix,
        extension=extension,
        read_group=read_group,
        sort=sort,
        args=args,
        args2=args2,
    )
    alignment = await results.get_file(f"{prefix}.{extension}")
    if alignment is None:
        raise FileNotFoundError(f"bwa-mem2 mem wrote no {prefix}.{extension}")
    return alignment
