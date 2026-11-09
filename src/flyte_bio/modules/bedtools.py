"""bedtools — genome arithmetic CLI tools.

Wraps tools from the `bedtools suite <https://bedtools.readthedocs.io/>`_
shipped via the ``quay.io/biocontainers/bedtools`` image.

Currently exposed:

- :data:`bedtools_intersect` — report overlaps between two feature files.
- :data:`bedtools_sort` — sort BED/GFF/VCF features by chromosome and start.
- :data:`bedtools_merge` — combine overlapping or nearby features into one.
- :func:`bedtools_genomecov` — genome-wide coverage of a BAM (e.g. bedGraph).
- :func:`bedtools_getfasta` — the FASTA sequences of a BED's intervals.
- :func:`bedtools_maskfasta` — a FASTA with a BED's intervals masked.
- :func:`bedtools_bamtobed` — a BAM's alignments as BED intervals.
- :func:`bedtools_complement` — the intervals of a genome a BED doesn't cover.
- :func:`bedtools_makewindows` — fixed-size windows over a BED's intervals or
  over a genome (a ``.fai`` or chromosome sizes file).

The module-level :data:`env` is a single :class:`flyte.TaskEnvironment`
containing every bedtools task. All commands share the same biocontainer
image, so they live in the same env. Pipelines depend on this one env to
gain access to every bedtools subcommand at once::

    from flyteplugins.bio.bedtools import (
        bedtools_intersect,
        bedtools_sort,
        bedtools_merge,
        env as bedtools_env,
    )

    env = flyte.TaskEnvironment(name="my_pipeline", depends_on=[bedtools_env])
"""



import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

# Pinned biocontainer URI. Update when bumping bedtools version.
BEDTOOLS_IMAGE = "quay.io/biocontainers/bedtools:2.31.1--hf5e1c6e_0"

# Most bedtools subcommands are single-threaded; memory is the dimension
# that varies with genome / interval-set size. Override per call when a
# specific run needs more.
DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory="6Gi")


bedtools_intersect = shell.create(
    name="bedtools_intersect",
    image=BEDTOOLS_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "a": File,
        "b": list[File],
        # Output mode (mutually exclusive in practice)
        "wa": bool,
        "wb": bool,
        "loj": bool,
        "wo": bool,
        "wao": bool,
        "u": bool,
        "c": bool,
        "C": bool,
        "v": bool,
        # BAM output / BAM-to-BED conversion
        "ubam": bool,
        "bed_output": bool,  # mapped to -bed
        # Strand
        "s": bool,
        "S": bool,
        # Overlap fraction
        "f": float | None,
        "F": float | None,
        "r": bool,
        "e": bool,
        # Split BAM/BED12 entries
        "split": bool,
        # Sorted-input optimisation
        "sorted": bool,
        "g": File | None,
        "nonamecheck": bool,
        # Multi-database controls
        "names": str | None,  # space-separated aliases, one per -b file
        "filenames": bool,
        "sortout": bool,
        # Output control
        "header": bool,
        "nobuf": bool,
        "iobuf": str | None,  # e.g. "128M", "1G"
    },
    flag_aliases={
        "bed_output": "-bed",
    },
    defaults={
        "wa": False, "wb": False, "loj": False, "wo": False, "wao": False,
        "u": False, "c": False, "C": False, "v": False,
        "ubam": False, "bed_output": False,
        "s": False, "S": False,
        "r": False, "e": False,
        "split": False,
        "sorted": False, "nonamecheck": False,
        "filenames": False, "sortout": False,
        "header": False, "nobuf": False,
    },
    outputs={
        "bed": File,
    },
    script=r"""
        bedtools intersect \
            {flags.wa} {flags.wb} {flags.loj} {flags.wo} {flags.wao} \
            {flags.u} {flags.c} {flags.C} {flags.v} \
            {flags.ubam} {flags.bed_output} \
            {flags.s} {flags.S} \
            {flags.f} {flags.F} {flags.r} {flags.e} \
            {flags.split} \
            {flags.sorted} {flags.g} {flags.nonamecheck} \
            {flags.names} {flags.filenames} {flags.sortout} \
            {flags.header} {flags.nobuf} {flags.iobuf} \
            -a {inputs.a} \
            -b {inputs.b} \
            > {outputs.bed}
    """,
)

bedtools_sort = shell.create(
    name="bedtools_sort",
    image=BEDTOOLS_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "i": File,
        # Sort orderings — `chrThenSizeA` and friends are bool flags.
        "header": bool,
        # Genome file: enforces a chromosome sort order across files.
        "g": File | None,
    },
    defaults={
        "header": False,
    },
    outputs={
        "sorted": File,
    },
    script=r"""
        bedtools sort \
            {flags.header} {flags.g} \
            -i {inputs.i} \
            > {outputs.sorted}
    """,
)


bedtools_merge = shell.create(
    name="bedtools_merge",
    image=BEDTOOLS_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "i": File,  # Must be sorted.
        "s": bool,  # Require same strand.
        "S": str | None,  # Force one strand: "+" or "-".
        "d": int | None,  # Max distance between features to merge.
        "header": bool,
    },
    defaults={
        "s": False,
        "header": False,
    },
    outputs={
        "merged": File,
    },
    script=r"""
        bedtools merge \
            {flags.s} {flags.S} {flags.d} {flags.header} \
            -i {inputs.i} \
            > {outputs.merged}
    """,
)

# genomecov pipes through GNU `sort`, so (as upstream) it runs on the
# bedtools + coreutils image. Output is `<prefix>.<extension>` in a Dir.
GENOMECOV_IMAGE = "community.wave.seqera.io/library/bedtools_coreutils:a623c13f66d5262b"
GENOMECOV_CPUS = 2
GENOMECOV_MEMORY_GB = 6

bedtools_genomecov_cmd = shell.create(
    name="bedtools_genomecov",
    image=GENOMECOV_IMAGE,
    resources=flyte.Resources(cpu=GENOMECOV_CPUS, memory=f"{GENOMECOV_MEMORY_GB}Gi"),
    inputs={"bam": File, "prefix": str, "extension": str, "sort": bool, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        BAM=({{inputs.bam}})
        ARGS={{inputs.args}}
        OUT={{outputs.results}}/{{inputs.prefix}}.{{inputs.extension}}
        if [ {{inputs.sort}} = true ]; then
            bedtools genomecov -ibam "${{BAM[0]}}" $ARGS \
                | LC_ALL=C sort --parallel=$TASK_CPUS --buffer-size=$(( TASK_MEMORY_MB / 2 ))M -k1,1 -k2,2n \
                > "$OUT"
        else
            bedtools genomecov -ibam "${{BAM[0]}}" $ARGS > "$OUT"
        fi
    """,
)


env = flyte.TaskEnvironment.from_task(
    "bedtools",
    bedtools_intersect.as_task(),
    bedtools_sort.as_task(),
    bedtools_merge.as_task(),
)

# genomecov's image differs, and a TaskEnvironment holds one image.
genomecov_env = flyte.TaskEnvironment.from_task("bedtools_genomecov", bedtools_genomecov_cmd.as_task())

# getfasta and maskfasta run on upstream's newer bedtools image. Both write
# `<prefix>.fa`; getfasta indexes the FASTA beside itself, so (as upstream
# stages it) the FASTA is linked into a scratch directory first.
FASTA_IMAGE = "community.wave.seqera.io/library/bedtools:2.31.1--897920cd303d6413"


def bedtools_fasta_task(subcommand: str):
    return shell.create(
        name=f"bedtools_{subcommand}",
        image=FASTA_IMAGE,
        resources=DEFAULT_RESOURCES,
        inputs={"bed": File, "fasta": File, "prefix": str, "args": str},
        defaults={"args": ""},
        outputs={"results": Dir},
        script=rf"""
        BED=({{inputs.bed}}); FA=({{inputs.fasta}})
        ARGS={{inputs.args}}
        OUT={{outputs.results}}/{{inputs.prefix}}.fa
        W=$(mktemp -d); cd "$W"
        FA_LINK=$(basename "${{FA[0]}}"); ln -s "${{FA[0]}}" "$FA_LINK"
        bedtools {subcommand} $ARGS -fi "$FA_LINK" -bed "${{BED[0]}}" -fo "$OUT"
    """,
    )


bedtools_getfasta_cmd = bedtools_fasta_task("getfasta")
bedtools_maskfasta_cmd = bedtools_fasta_task("maskfasta")

# bamtobed and complement run on the same image; each writes `<prefix>.bed`.
bedtools_bamtobed_cmd = shell.create(
    name="bedtools_bamtobed",
    image=FASTA_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"bam": File, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        BAM=({inputs.bam})
        ARGS={inputs.args}
        bedtools bamtobed $ARGS -i "${BAM[0]}" > {outputs.results}/{inputs.prefix}.bed
    """,
)

bedtools_complement_cmd = shell.create(
    name="bedtools_complement",
    image=FASTA_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"bed": File, "sizes": File, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        BED=({inputs.bed}); SIZES=({inputs.sizes})
        ARGS={inputs.args}
        bedtools complement -i "${BED[0]}" -g "${SIZES[0]}" $ARGS > {outputs.results}/{inputs.prefix}.bed
    """,
)

# As upstream, a .bed or .tab input is a BED of intervals (-b); anything else
# is a genome file such as a .fai (-g).
bedtools_makewindows_cmd = shell.create(
    name="bedtools_makewindows",
    image=FASTA_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"regions": File, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        REGIONS=({inputs.regions})
        ARGS={inputs.args}
        case "${REGIONS[0]}" in
            *.bed | *.tab) INPUT=(-b "${REGIONS[0]}") ;;
            *) INPUT=(-g "${REGIONS[0]}") ;;
        esac
        bedtools makewindows "${INPUT[@]}" $ARGS > {outputs.results}/{inputs.prefix}.bed || exit 1
    """,
)

# `env` stays the one to depend on: it pulls in the newer image's env too.
fasta_env = flyte.TaskEnvironment.from_task(
    "bedtools_fasta",
    bedtools_getfasta_cmd.as_task(),
    bedtools_maskfasta_cmd.as_task(),
    bedtools_bamtobed_cmd.as_task(),
    bedtools_complement_cmd.as_task(),
    bedtools_makewindows_cmd.as_task(),
)
env.add_dependency(fasta_env)


async def bedtools_genomecov(
    bam: File, prefix: str, extension: str = "bedGraph", sort: bool = True, args: str = ""
) -> File:
    """Coverage of ``bam`` as ``<prefix>.<extension>``, sorted by chrom/start when ``sort`` (as upstream)."""
    results = await bedtools_genomecov_cmd(bam=bam, prefix=prefix, extension=extension, sort=sort, args=args)
    out = await results.get_file(f"{prefix}.{extension}")
    if out is None:
        raise FileNotFoundError(f"bedtools genomecov wrote no {prefix}.{extension}")
    return out


async def bedtools_getfasta(bed: File, fasta: File, prefix: str, args: str = "") -> File:
    """The sequences of ``fasta`` under ``bed``'s intervals, as ``<prefix>.fa``."""
    if fasta.name == f"{prefix}.fa":
        raise ValueError("input and output names are the same; choose another prefix")
    results = await bedtools_getfasta_cmd(bed=bed, fasta=fasta, prefix=prefix, args=args)
    out = await results.get_file(f"{prefix}.fa")
    if out is None:
        raise FileNotFoundError(f"bedtools getfasta wrote no {prefix}.fa")
    return out


async def bedtools_maskfasta(bed: File, fasta: File, prefix: str, args: str = "") -> File:
    """``fasta`` with ``bed``'s intervals masked (soft-masked with ``args="-soft"``), as ``<prefix>.fa``."""
    results = await bedtools_maskfasta_cmd(bed=bed, fasta=fasta, prefix=prefix, args=args)
    out = await results.get_file(f"{prefix}.fa")
    if out is None:
        raise FileNotFoundError(f"bedtools maskfasta wrote no {prefix}.fa")
    return out


async def bedtools_bamtobed(bam: File, prefix: str, args: str = "") -> File:
    """The alignments of ``bam`` as BED intervals, ``<prefix>.bed``."""
    results = await bedtools_bamtobed_cmd(bam=bam, prefix=prefix, args=args)
    out = await results.get_file(f"{prefix}.bed")
    if out is None:
        raise FileNotFoundError(f"bedtools bamtobed wrote no {prefix}.bed")
    return out


async def bedtools_complement(bed: File, sizes: File, prefix: str, args: str = "") -> File:
    """The intervals of the genome in ``sizes`` (chromosome, length) that ``bed`` doesn't cover, ``<prefix>.bed``."""
    if bed.name == f"{prefix}.bed":
        raise ValueError("input and output names are the same; choose another prefix")
    results = await bedtools_complement_cmd(bed=bed, sizes=sizes, prefix=prefix, args=args)
    out = await results.get_file(f"{prefix}.bed")
    if out is None:
        raise FileNotFoundError(f"bedtools complement wrote no {prefix}.bed")
    return out


async def bedtools_makewindows(regions: File, prefix: str, args: str = "") -> File:
    """Windows over ``regions``, as ``<prefix>.bed``; set their size with ``-w`` (or ``-n``) in ``args``.

    ``regions`` is a BED of intervals (``.bed`` or ``.tab``) or a genome file
    of chromosome lengths, such as a ``.fai``.
    """
    if regions.name == f"{prefix}.bed":
        raise ValueError("Input and output names are the same, set prefix to disambiguate!")
    results = await bedtools_makewindows_cmd(regions=regions, prefix=prefix, args=args)
    out = await results.get_file(f"{prefix}.bed")
    if out is None:
        raise FileNotFoundError(f"bedtools makewindows wrote no {prefix}.bed")
    return out
