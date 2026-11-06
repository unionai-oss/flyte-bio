"""bcftools — VCF/BCF utilities.

Exposes:

- :func:`bcftools_stats` — ``bcftools stats`` of a VCF, written to
  ``<prefix>.bcftools_stats.txt``. The VCF (and a reference FASTA, for indel
  context) are linked into the working directory under their own names, so the
  command line bcftools records matches upstream's.
- :func:`bcftools_sort` — sort a VCF/BCF.
- :func:`bcftools_norm` — normalize indels against a reference and split or
  join multiallelic sites.
- :func:`bcftools_view` — subset and convert a VCF/BCF, optionally by regions,
  targets or samples files.
- :func:`bcftools_filter` — filter a VCF/BCF.
- :func:`bcftools_annotate` — add, remove or rename annotations and contigs.
- :func:`bcftools_merge` — merge several samples' VCFs into one.
- :func:`bcftools_mpileup` — call variants from a BAM: ``bcftools mpileup``
  piped into ``bcftools call``, the sample renamed, then ``bcftools stats`` of
  the calls.

All but :func:`bcftools_mpileup` (whose ``args3`` play the part) take the output
format from ``--output-type``/``-O`` in ``args``, as upstream, and write
``<prefix>.vcf``, ``.vcf.gz``, ``.bcf`` or ``.bcf.gz``; with ``--write-index``
in the args the index lands beside it as ``<prefix>.<extension>.csi`` (or
``.tbi`` with ``--write-index=tbi``). Inputs are linked into a working
directory under their own names, indexes next to their VCFs, and the tools
run there, so the command lines bcftools records match upstream's. Tools that
take ``--threads`` use the task's CPUs; :func:`bcftools_sort` gets 90% of the
task's memory as ``--max-mem``.

bcftools stats keeps the bcftools 1.23.1 image the variant-calling pipeline
uses; the other tools use upstream's bcftools 1.25 image.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

BCFTOOLS_IMAGE = "community.wave.seqera.io/library/bcftools_htslib:1.23.1--9f08ec665533d64a"
BCFTOOLS_TOOLS_IMAGE = "community.wave.seqera.io/library/bcftools_htslib:1.25--130aa7603bf3bf5c"

DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory="2Gi")
# Upstream labels these tools process_medium (6 CPUs, 36 GB); test-sized VCFs
# need far less. Raise with task.override(resources=...) for whole genomes.
TOOLS_RESOURCES = flyte.Resources(cpu=2, memory="4Gi")

bcftools_stats_cmd = shell.create(
    name="bcftools_stats",
    image=BCFTOOLS_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"vcf": File, "tbi": File | None, "fasta": File | None, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        V=({inputs.vcf})
        shopt -s nullglob; T=({inputs.tbi}); FA=({inputs.fasta}); shopt -u nullglob
        ARGS={inputs.args}
        W=$(mktemp -d); cd "$W"
        NAME=$(basename "${V[0]}"); ln -s "${V[0]}" "$NAME"
        if [ ${#T[@]} -gt 0 ]; then ln -s "${T[0]}" "$NAME.tbi"; fi
        REF=()
        if [ ${#FA[@]} -gt 0 ]; then
            FA_NAME=$(basename "${FA[0]}"); ln -s "${FA[0]}" "$FA_NAME"; REF=(--fasta-ref "$FA_NAME")
        fi
        bcftools stats $ARGS "${REF[@]}" "$NAME" > {outputs.results}/{inputs.prefix}.bcftools_stats.txt
    """,
)

# Linking helpers, for a working directory: `link FILE` links FILE under its
# own name and prints that name; `link_index FILE INDEX` links INDEX beside
# FILE as FILE.csi or FILE.tbi.
LINK_HELPERS = r"""
        link() { ln -s "$1" "$(basename "$1")"; basename "$1"; }
        link_index() { case "$2" in *.csi) ln -s "$2" "$1.csi" ;; *) ln -s "$2" "$1.tbi" ;; esac; }
"""

# Starts a script in a fresh working directory with the VCF (array V) linked
# under its own name, NAME, and its index (array I, may be empty) beside it.
LINK_VCF = (
    LINK_HELPERS
    + r"""
        W=$(mktemp -d); cd "$W"
        NAME=$(link "${V[0]}")
        if [ ${#I[@]} -gt 0 ]; then link_index "$NAME" "${I[0]}"; fi
"""
)

# Moves the output, and the index --write-index wrote beside it, into OUT.
COLLECT_OUTPUT = r"""
        mv "$P.$EXT" "$OUT"/
        for SUFFIX in csi tbi; do
            if [ -f "$P.$EXT.$SUFFIX" ]; then mv "$P.$EXT.$SUFFIX" "$OUT"/; fi
        done
"""

# Every tool but mpileup writes <prefix>.<extension> (from the wrapper) into
# one results Dir.
VCF_IO = {"prefix": str, "extension": str, "args": str}

bcftools_sort_cmd = shell.create(
    name="bcftools_sort",
    image=BCFTOOLS_TOOLS_IMAGE,
    resources=TOOLS_RESOURCES,
    inputs={"vcf": File, **VCF_IO},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        V=({{inputs.vcf}}); I=()
        P={{inputs.prefix}}; EXT={{inputs.extension}}; OUT={{outputs.results}}
        ARGS={{inputs.args}}
        {LINK_VCF}
        bcftools sort --output "$P.$EXT" --temp-dir . --max-mem "$(( TASK_MEMORY_MB * 9 / 10 ))M" $ARGS "$NAME"
        {COLLECT_OUTPUT}
    """,
)

bcftools_norm_cmd = shell.create(
    name="bcftools_norm",
    image=BCFTOOLS_TOOLS_IMAGE,
    resources=TOOLS_RESOURCES,
    inputs={"vcf": File, "index": File | None, "fasta": File, **VCF_IO},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        V=({{inputs.vcf}}); FA=({{inputs.fasta}})
        shopt -s nullglob; I=({{inputs.index}}); shopt -u nullglob
        P={{inputs.prefix}}; EXT={{inputs.extension}}; OUT={{outputs.results}}
        ARGS={{inputs.args}}
        {LINK_VCF}
        FA_NAME=$(link "${{FA[0]}}")
        bcftools norm --fasta-ref "$FA_NAME" --output "$P.$EXT" $ARGS --threads $TASK_CPUS "$NAME"
        {COLLECT_OUTPUT}
    """,
)

bcftools_view_cmd = shell.create(
    name="bcftools_view",
    image=BCFTOOLS_TOOLS_IMAGE,
    resources=TOOLS_RESOURCES,
    inputs={
        "vcf": File,
        "index": File | None,
        "regions": File | None,
        "targets": File | None,
        "samples": File | None,
        "include": str,
        "exclude": str,
        **VCF_IO,
    },
    defaults={"include": "", "exclude": "", "args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        V=({{inputs.vcf}})
        shopt -s nullglob
        I=({{inputs.index}}); R=({{inputs.regions}}); T=({{inputs.targets}}); S=({{inputs.samples}})
        shopt -u nullglob
        P={{inputs.prefix}}; EXT={{inputs.extension}}; OUT={{outputs.results}}
        set -f
        ARGS={{inputs.args}}
        {LINK_VCF}
        SUBSET=()
        if [ ${{#R[@]}} -gt 0 ]; then SUBSET+=(--regions-file "$(link "${{R[0]}}")"); fi
        if [ ${{#T[@]}} -gt 0 ]; then SUBSET+=(--targets-file "$(link "${{T[0]}}")"); fi
        if [ ${{#S[@]}} -gt 0 ]; then SUBSET+=(--samples-file "$(link "${{S[0]}}")"); fi
        EXPR=()
        INCLUDE={{inputs.include}}; EXCLUDE={{inputs.exclude}}
        if [ -n "$INCLUDE" ]; then EXPR+=(--include "$INCLUDE"); fi
        if [ -n "$EXCLUDE" ]; then EXPR+=(--exclude "$EXCLUDE"); fi
        bcftools view --output "$P.$EXT" "${{SUBSET[@]}}" "${{EXPR[@]}}" $ARGS --threads $TASK_CPUS "$NAME"
        {COLLECT_OUTPUT}
    """,
)

bcftools_filter_cmd = shell.create(
    name="bcftools_filter",
    image=BCFTOOLS_TOOLS_IMAGE,
    resources=TOOLS_RESOURCES,
    inputs={"vcf": File, "index": File | None, "include": str, "exclude": str, **VCF_IO},
    defaults={"include": "", "exclude": "", "args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        V=({{inputs.vcf}})
        shopt -s nullglob; I=({{inputs.index}}); shopt -u nullglob
        P={{inputs.prefix}}; EXT={{inputs.extension}}; OUT={{outputs.results}}
        set -f
        ARGS={{inputs.args}}
        {LINK_VCF}
        EXPR=()
        INCLUDE={{inputs.include}}; EXCLUDE={{inputs.exclude}}
        if [ -n "$INCLUDE" ]; then EXPR+=(--include "$INCLUDE"); fi
        if [ -n "$EXCLUDE" ]; then EXPR+=(--exclude "$EXCLUDE"); fi
        bcftools filter --output "$P.$EXT" --threads $TASK_CPUS "${{EXPR[@]}}" $ARGS "$NAME"
        {COLLECT_OUTPUT}
    """,
)

# As upstream, a VCF without an index is indexed first (annotate needs one).
bcftools_annotate_cmd = shell.create(
    name="bcftools_annotate",
    image=BCFTOOLS_TOOLS_IMAGE,
    resources=TOOLS_RESOURCES,
    inputs={
        "vcf": File,
        "index": File | None,
        "annotations": File | None,
        "annotations_index": File | None,
        "columns": File | None,
        "header_lines": File | None,
        "rename_chrs": File | None,
        **VCF_IO,
    },
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        V=({{inputs.vcf}})
        shopt -s nullglob
        I=({{inputs.index}}); A=({{inputs.annotations}}); AI=({{inputs.annotations_index}})
        C=({{inputs.columns}}); H=({{inputs.header_lines}}); RC=({{inputs.rename_chrs}})
        shopt -u nullglob
        P={{inputs.prefix}}; EXT={{inputs.extension}}; OUT={{outputs.results}}
        set -f
        ARGS={{inputs.args}}
        {LINK_VCF}
        if [ ${{#I[@]}} -eq 0 ]; then bcftools index "$NAME"; fi
        EXTRA=()
        if [ ${{#A[@]}} -gt 0 ]; then
            A_NAME=$(link "${{A[0]}}"); EXTRA+=(--annotations "$A_NAME")
            if [ ${{#AI[@]}} -gt 0 ]; then link_index "$A_NAME" "${{AI[0]}}"; fi
        fi
        if [ ${{#C[@]}} -gt 0 ]; then EXTRA+=(--columns-file "$(link "${{C[0]}}")"); fi
        if [ ${{#H[@]}} -gt 0 ]; then EXTRA+=(--header-lines "$(link "${{H[0]}}")"); fi
        if [ ${{#RC[@]}} -gt 0 ]; then EXTRA+=(--rename-chrs "$(link "${{RC[0]}}")"); fi
        bcftools annotate $ARGS "${{EXTRA[@]}}" --output "$P.$EXT" --threads $TASK_CPUS "$NAME"
        {COLLECT_OUTPUT}
    """,
)

# The VCFs (in name order, as upstream) and their indexes are linked side by
# side; so are the FASTA and its .fai, which args such as --gvcf name.
bcftools_merge_cmd = shell.create(
    name="bcftools_merge",
    image=BCFTOOLS_TOOLS_IMAGE,
    resources=TOOLS_RESOURCES,
    inputs={
        "vcfs": list[File],
        "indexes": list[File],
        "bed": File | None,
        "fasta": File | None,
        "fai": File | None,
        **VCF_IO,
    },
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        shopt -s nullglob
        VCFS=({{inputs.vcfs}}); IDX=({{inputs.indexes}})
        B=({{inputs.bed}}); FA=({{inputs.fasta}}); FAI=({{inputs.fai}})
        shopt -u nullglob
        P={{inputs.prefix}}; EXT={{inputs.extension}}; OUT={{outputs.results}}
        set -f
        ARGS={{inputs.args}}
        {LINK_HELPERS}
        W=$(mktemp -d); cd "$W"
        NAMES=()
        for F in "${{VCFS[@]}}"; do NAMES+=("$(link "$F")"); done
        for F in "${{IDX[@]}}" "${{FA[@]}}" "${{FAI[@]}}"; do link "$F" > /dev/null; done
        REGIONS=()
        if [ ${{#B[@]}} -gt 0 ]; then REGIONS=(--regions-file "$(link "${{B[0]}}")"); fi
        bcftools merge $ARGS "${{REGIONS[@]}}" --threads $TASK_CPUS --output "$P.$EXT" "${{NAMES[@]}}"
        {COLLECT_OUTPUT}
    """,
)

# As upstream: mpileup | (tee <prefix>.mpileup) | call | reheader to the sample
# name | view, then bcftools stats of the calls. Interval files are linked under
# mpileup_intervals/ and call_intervals/, where upstream stages them.
bcftools_mpileup_cmd = shell.create(
    name="bcftools_mpileup",
    image=BCFTOOLS_TOOLS_IMAGE,
    resources=TOOLS_RESOURCES,
    inputs={
        "bam": File,
        "intervals_mpileup": File | None,
        "intervals_call": File | None,
        "fasta": File,
        "fai": File | None,
        "sample": str,
        "save_mpileup": bool,
        "prefix": str,
        "extension": str,
        "args": str,
        "args2": str,
        "args3": str,
    },
    defaults={"save_mpileup": False, "args": "", "args2": "", "args3": ""},
    outputs={"results": Dir},
    script=rf"""
        set -o pipefail
        BAM=({{inputs.bam}}); FA=({{inputs.fasta}})
        shopt -s nullglob
        IM=({{inputs.intervals_mpileup}}); IC=({{inputs.intervals_call}}); FAI=({{inputs.fai}})
        shopt -u nullglob
        P={{inputs.prefix}}; EXT={{inputs.extension}}; OUT={{outputs.results}}
        set -f
        ARGS={{inputs.args}}; ARGS2={{inputs.args2}}; ARGS3={{inputs.args3}}
        W=$(mktemp -d); cd "$W"
        BAM_NAME=$(basename "${{BAM[0]}}"); ln -s "${{BAM[0]}}" "$BAM_NAME"
        FA_NAME=$(basename "${{FA[0]}}"); ln -s "${{FA[0]}}" "$FA_NAME"
        if [ ${{#FAI[@]}} -gt 0 ]; then ln -s "${{FAI[0]}}" "$FA_NAME.fai"; fi
        MPILEUP_TARGETS=(); CALL_TARGETS=()
        if [ ${{#IM[@]}} -gt 0 ]; then
            mkdir mpileup_intervals; ln -s "${{IM[0]}}" "mpileup_intervals/$(basename "${{IM[0]}}")"
            MPILEUP_TARGETS=(-T "mpileup_intervals/$(basename "${{IM[0]}}")")
        fi
        if [ ${{#IC[@]}} -gt 0 ]; then
            mkdir call_intervals; ln -s "${{IC[0]}}" "call_intervals/$(basename "${{IC[0]}}")"
            CALL_TARGETS=(-T "call_intervals/$(basename "${{IC[0]}}")")
        fi
        echo {{inputs.sample}} > sample_name.list
        SAVE=(cat)
        if [ {{inputs.save_mpileup}} = true ]; then SAVE=(tee "$P.mpileup"); fi
        bcftools mpileup --fasta-ref "$FA_NAME" $ARGS "$BAM_NAME" "${{MPILEUP_TARGETS[@]}}" \
            | "${{SAVE[@]}}" \
            | bcftools call --output-type v $ARGS2 "${{CALL_TARGETS[@]}}" \
            | bcftools reheader --samples sample_name.list \
            | bcftools view --output-file "$P.$EXT" $ARGS3
        if [ {{inputs.save_mpileup}} = true ]; then bgzip "$P.mpileup"; mv "$P.mpileup.gz" "$OUT"/; fi
        bcftools stats "$P.$EXT" > "$OUT/$P.bcftools_stats.txt"
        {COLLECT_OUTPUT}
    """,
)

bcftools_stats_env = flyte.TaskEnvironment.from_task("bcftools_stats", bcftools_stats_cmd.as_task())
bcftools_tools_env = flyte.TaskEnvironment.from_task(
    "bcftools_tools",
    bcftools_sort_cmd.as_task(),
    bcftools_norm_cmd.as_task(),
    bcftools_view_cmd.as_task(),
    bcftools_filter_cmd.as_task(),
    bcftools_annotate_cmd.as_task(),
    bcftools_merge_cmd.as_task(),
    bcftools_mpileup_cmd.as_task(),
)
env = flyte.TaskEnvironment(name="bcftools", depends_on=[bcftools_stats_env, bcftools_tools_env])


async def bcftools_stats(
    vcf: File, prefix: str, tbi: File | None = None, fasta: File | None = None, args: str = ""
) -> File:
    """``<prefix>.bcftools_stats.txt``: bcftools' statistics of ``vcf``."""
    results = await bcftools_stats_cmd(vcf=vcf, tbi=tbi, fasta=fasta, prefix=prefix, args=args)
    stats = await results.get_file(f"{prefix}.bcftools_stats.txt")
    if stats is None:
        raise FileNotFoundError(f"bcftools stats wrote no {prefix}.bcftools_stats.txt")
    return stats


def output_extension(args: str, default: str = "vcf") -> str:
    """The output's extension, read from ``--output-type``/``-O`` in ``args`` as upstream does."""
    for options, extension in (
        (("--output-type b", "-Ob"), "bcf.gz"),
        (("--output-type u", "-Ou"), "bcf"),
        (("--output-type z", "-Oz"), "vcf.gz"),
        (("--output-type v", "-Ov"), "vcf"),
    ):
        if any(option in args for option in options):
            return extension
    return default


def check_output_name(vcf: File, prefix: str, extension: str) -> None:
    """Refuse to write the output under the input's own name, as upstream does."""
    if vcf.name == f"{prefix}.{extension}":
        raise ValueError("Input and output names are the same, set prefix to disambiguate!")


@dataclass
class BcftoolsResult:
    vcf: File  # <prefix>.vcf, .vcf.gz, .bcf or .bcf.gz
    index: File | None  # <prefix>.<extension>.csi or .tbi, with --write-index in the args
    results: Dir  # everything the task wrote


async def collect(results: Dir, prefix: str, extension: str, tool: str) -> BcftoolsResult:
    vcf = await results.get_file(f"{prefix}.{extension}")
    if vcf is None:
        raise FileNotFoundError(f"bcftools {tool} wrote no {prefix}.{extension}")
    index = None
    for suffix in ("csi", "tbi"):
        index = index or await results.get_file(f"{prefix}.{extension}.{suffix}")
    return BcftoolsResult(vcf=vcf, index=index, results=results)


async def bcftools_sort(vcf: File, prefix: str, args: str = "--output-type z") -> BcftoolsResult:
    """Sort ``vcf`` into ``<prefix>.<extension>`` (by default bgzipped VCF, as upstream)."""
    extension = output_extension(args)
    check_output_name(vcf, prefix, extension)
    results = await bcftools_sort_cmd(vcf=vcf, prefix=prefix, extension=extension, args=args)
    return await collect(results, prefix, extension, "sort")


async def bcftools_norm(
    vcf: File, fasta: File, prefix: str, index: File | None = None, args: str = "--output-type z"
) -> BcftoolsResult:
    """Normalize ``vcf`` against ``fasta`` (plain or bgzipped) into ``<prefix>.<extension>``.

    The output is a bgzipped VCF unless ``args`` set ``--output-type``.
    """
    extension = output_extension(args, default="vcf.gz")
    check_output_name(vcf, prefix, extension)
    results = await bcftools_norm_cmd(vcf=vcf, index=index, fasta=fasta, prefix=prefix, extension=extension, args=args)
    return await collect(results, prefix, extension, "norm")


async def bcftools_view(
    vcf: File,
    prefix: str,
    index: File | None = None,
    regions: File | None = None,
    targets: File | None = None,
    samples: File | None = None,
    include: str = "",
    exclude: str = "",
    args: str = "",
) -> BcftoolsResult:
    """``bcftools view`` of ``vcf`` into ``<prefix>.<extension>``.

    ``regions``, ``targets`` and ``samples`` are passed as ``--regions-file``,
    ``--targets-file`` and ``--samples-file``; ``regions`` needs ``index``.
    ``include`` and ``exclude`` are filter expressions (``--include``/``--exclude``),
    passed whole, so unlike ``args`` they can contain spaces: ``"QUAL>30 && DP>10"``.
    """
    extension = output_extension(args)
    check_output_name(vcf, prefix, extension)
    results = await bcftools_view_cmd(
        vcf=vcf,
        index=index,
        regions=regions,
        targets=targets,
        samples=samples,
        include=include,
        exclude=exclude,
        prefix=prefix,
        extension=extension,
        args=args,
    )
    return await collect(results, prefix, extension, "view")


async def bcftools_filter(
    vcf: File, prefix: str, index: File | None = None, include: str = "", exclude: str = "", args: str = ""
) -> BcftoolsResult:
    """Filter ``vcf`` into ``<prefix>.<extension>``; region options (``-r``) need ``index``.

    ``include`` and ``exclude`` are filter expressions (``--include``/``--exclude``),
    passed whole, so unlike ``args`` they can contain spaces: ``"QUAL>30 && DP>10"``.
    """
    extension = output_extension(args)
    check_output_name(vcf, prefix, extension)
    results = await bcftools_filter_cmd(
        vcf=vcf, index=index, include=include, exclude=exclude, prefix=prefix, extension=extension, args=args
    )
    return await collect(results, prefix, extension, "filter")


async def bcftools_annotate(
    vcf: File,
    prefix: str,
    index: File | None = None,
    annotations: File | None = None,
    annotations_index: File | None = None,
    columns: File | None = None,
    header_lines: File | None = None,
    rename_chrs: File | None = None,
    args: str = "",
) -> BcftoolsResult:
    """Annotate ``vcf`` into ``<prefix>.<extension>``.

    ``annotations`` (with its index), ``columns``, ``header_lines`` and
    ``rename_chrs`` are passed as ``--annotations``, ``--columns-file``,
    ``--header-lines`` and ``--rename-chrs``. Without ``index``, ``vcf`` is
    indexed first.
    """
    extension = output_extension(args)
    check_output_name(vcf, prefix, extension)
    results = await bcftools_annotate_cmd(
        vcf=vcf,
        index=index,
        annotations=annotations,
        annotations_index=annotations_index,
        columns=columns,
        header_lines=header_lines,
        rename_chrs=rename_chrs,
        prefix=prefix,
        extension=extension,
        args=args,
    )
    return await collect(results, prefix, extension, "annotate")


async def bcftools_merge(
    vcfs: list[File],
    indexes: list[File],
    prefix: str,
    bed: File | None = None,
    fasta: File | None = None,
    fai: File | None = None,
    args: str = "",
) -> BcftoolsResult:
    """Merge ``vcfs`` (with their ``.tbi``/``.csi`` ``indexes``) into ``<prefix>.<extension>``.

    ``bed`` restricts the merge (``--regions-file``). ``fasta`` and ``fai`` are
    linked into the working directory for args that name them, such as
    ``--gvcf genome.fasta``.
    """
    extension = output_extension(args)
    results = await bcftools_merge_cmd(
        vcfs=vcfs, indexes=indexes, bed=bed, fasta=fasta, fai=fai, prefix=prefix, extension=extension, args=args
    )
    return await collect(results, prefix, extension, "merge")


@dataclass
class BcftoolsMpileupResult:
    vcf: File  # <prefix>.vcf, .vcf.gz, .bcf or .bcf.gz
    index: File | None  # <prefix>.<extension>.csi or .tbi, with --write-index in args3
    stats: File  # <prefix>.bcftools_stats.txt
    mpileup: File | None  # <prefix>.mpileup.gz, with save_mpileup
    results: Dir  # everything the task wrote


async def bcftools_mpileup(
    bam: File,
    fasta: File,
    prefix: str,
    sample: str | None = None,
    intervals_mpileup: File | None = None,
    intervals_call: File | None = None,
    fai: File | None = None,
    save_mpileup: bool = False,
    args: str = "",
    args2: str = "",
    args3: str = "",
) -> BcftoolsMpileupResult:
    """Call variants in ``bam`` against ``fasta``; ``args`` go to mpileup, ``args2`` to call, ``args3`` to view.

    The calls' sample is renamed ``sample`` (default ``prefix``). The
    ``intervals_*`` files restrict mpileup and call (``-T``). ``args3``'s
    ``--output-type`` sets the output format (default VCF).
    """
    extension = output_extension(args3)
    results = await bcftools_mpileup_cmd(
        bam=bam,
        intervals_mpileup=intervals_mpileup,
        intervals_call=intervals_call,
        fasta=fasta,
        fai=fai,
        sample=sample or prefix,
        save_mpileup=save_mpileup,
        prefix=prefix,
        extension=extension,
        args=args,
        args2=args2,
        args3=args3,
    )
    calls = await collect(results, prefix, extension, "mpileup")
    stats = await results.get_file(f"{prefix}.bcftools_stats.txt")
    if stats is None:
        raise FileNotFoundError(f"bcftools mpileup wrote no {prefix}.bcftools_stats.txt")
    mpileup = await results.get_file(f"{prefix}.mpileup.gz") if save_mpileup else None
    if save_mpileup and mpileup is None:
        raise FileNotFoundError(f"bcftools mpileup wrote no {prefix}.mpileup.gz")
    return BcftoolsMpileupResult(vcf=calls.vcf, index=calls.index, stats=stats, mpileup=mpileup, results=results)
