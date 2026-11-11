"""survivor — SURVIVOR, a toolkit for structural variant VCFs.

Exposes:

- :func:`survivor_filter` — keep a VCF's SVs by size, allele frequency and
  read support (``RE``), dropping those in a BED's regions:
  ``<prefix>.vcf``.
- :func:`survivor_merge` — merge (plain) VCFs from several callers or samples
  into ``<prefix>.vcf``, matching SVs whose breakpoints lie within a
  distance.
- :func:`survivor_stats` — summary statistics of a VCF's SVs (by type and
  size): ``<prefix>.stats``. A gzipped VCF is decompressed first, as
  upstream.

As upstream, the numeric thresholds are positional SURVIVOR arguments, ``-1``
disabling a size or read-support limit.
"""

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

# Upstream's main.nf names the bare biocontainer (Docker Hub); this is its quay.io mirror.
SURVIVOR_IMAGE = "quay.io/biocontainers/survivor:1.0.7--h9a82719_1"
# Upstream's process_single and process_low, scaled to the tool: SURVIVOR is single-threaded.
DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory="2Gi")

survivor_filter_cmd = shell.create(
    name="survivor_filter",
    image=SURVIVOR_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "vcf": File,
        "bed": File | None,
        "min_sv_size": int,
        "max_sv_size": int,
        "min_allele_freq": float,
        "min_read_support": int,
        "prefix": str,
    },
    outputs={"results": Dir},
    script=r"""
        VCF=({inputs.vcf})
        shopt -s nullglob; BED=({inputs.bed}); shopt -u nullglob
        BED_FILE=NA
        if [ ${#BED[@]} -gt 0 ]; then BED_FILE="${BED[0]}"; fi
        SURVIVOR filter "${VCF[0]}" "$BED_FILE" {inputs.min_sv_size} {inputs.max_sv_size} \
            {inputs.min_allele_freq} {inputs.min_read_support} {outputs.results}/{inputs.prefix}.vcf || exit 1
    """,
)

# As upstream, the VCFs are listed in name order (`ls *.vcf`) in the file SURVIVOR reads them from.
survivor_merge_cmd = shell.create(
    name="survivor_merge",
    image=SURVIVOR_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "vcfs": list[File],
        "max_distance_breakpoints": float,
        "min_supporting_callers": int,
        "account_for_type": int,
        "account_for_sv_strands": int,
        "estimate_distance_by_sv_size": int,
        "min_sv_size": int,
        "prefix": str,
    },
    outputs={"results": Dir},
    script=r"""
        VCFS=({inputs.vcfs})
        W=$(mktemp -d); cd "$W"
        mkdir vcfs
        for f in "${VCFS[@]}"; do ln -s "$f" "vcfs/$(basename "$f")"; done
        (cd vcfs && ls *.vcf) | sed "s|^|$W/vcfs/|" > vcf_list.txt
        SURVIVOR merge vcf_list.txt {inputs.max_distance_breakpoints} {inputs.min_supporting_callers} \
            {inputs.account_for_type} {inputs.account_for_sv_strands} {inputs.estimate_distance_by_sv_size} \
            {inputs.min_sv_size} {outputs.results}/{inputs.prefix}.vcf || exit 1
    """,
)

survivor_stats_cmd = shell.create(
    name="survivor_stats",
    image=SURVIVOR_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"vcf": File, "min_sv_size": int, "max_sv_size": int, "min_read_support": int, "prefix": str},
    outputs={"results": Dir},
    script=r"""
        VCF=({inputs.vcf})
        W=$(mktemp -d); cd "$W"
        NAME=$(basename "${VCF[0]}")
        if [[ "$NAME" == *.gz ]]; then
            NAME="${NAME%.gz}"
            gzip -c -d "${VCF[0]}" > "$NAME"
        else
            ln -s "${VCF[0]}" "$NAME"
        fi
        SURVIVOR stats "$NAME" {inputs.min_sv_size} {inputs.max_sv_size} {inputs.min_read_support} \
            {outputs.results}/{inputs.prefix}.stats || exit 1
    """,
)

env = flyte.TaskEnvironment.from_task(
    "survivor", survivor_filter_cmd.as_task(), survivor_merge_cmd.as_task(), survivor_stats_cmd.as_task()
)


async def pick(results: Dir, name: str, tool: str) -> File:
    found = await results.get_file(name)
    if found is None:
        raise FileNotFoundError(f"SURVIVOR {tool} wrote no {name}")
    return found


async def survivor_filter(
    vcf: File,
    prefix: str,
    min_sv_size: int,
    max_sv_size: int,
    min_allele_freq: float,
    min_read_support: int,
    bed: File | None = None,
) -> File:
    """``vcf``'s SVs within the size, allele-frequency and read-support limits, outside ``bed``: ``<prefix>.vcf``."""
    if vcf.name == f"{prefix}.vcf":
        raise ValueError("input and output names are the same; choose another prefix")
    results = await survivor_filter_cmd(
        vcf=vcf,
        bed=bed,
        min_sv_size=min_sv_size,
        max_sv_size=max_sv_size,
        min_allele_freq=min_allele_freq,
        min_read_support=min_read_support,
        prefix=prefix,
    )
    return await pick(results, f"{prefix}.vcf", "filter")


async def survivor_merge(
    vcfs: list[File],
    prefix: str,
    max_distance_breakpoints: float,
    min_supporting_callers: int,
    account_for_type: int,
    account_for_sv_strands: int,
    estimate_distance_by_sv_size: int,
    min_sv_size: int,
) -> File:
    """Merge the (uncompressed, distinctly named) ``vcfs`` into ``<prefix>.vcf``.

    ``max_distance_breakpoints`` is in bases, or a fraction of the SV's size
    when below 1; the ``account_for_*`` and ``estimate_*`` flags are 1 or 0.
    """
    if any(v.name.endswith(".gz") for v in vcfs):
        raise ValueError("SURVIVOR merge needs uncompressed VCFs")
    results = await survivor_merge_cmd(
        vcfs=vcfs,
        max_distance_breakpoints=max_distance_breakpoints,
        min_supporting_callers=min_supporting_callers,
        account_for_type=account_for_type,
        account_for_sv_strands=account_for_sv_strands,
        estimate_distance_by_sv_size=estimate_distance_by_sv_size,
        min_sv_size=min_sv_size,
        prefix=prefix,
    )
    return await pick(results, f"{prefix}.vcf", "merge")


async def survivor_stats(
    vcf: File, prefix: str, min_sv_size: int = -1, max_sv_size: int = -1, min_read_support: int = -1
) -> File:
    """Summary statistics of ``vcf``'s SVs (plain or gzipped): ``<prefix>.stats``."""
    results = await survivor_stats_cmd(
        vcf=vcf, min_sv_size=min_sv_size, max_sv_size=max_sv_size, min_read_support=min_read_support, prefix=prefix
    )
    return await pick(results, f"{prefix}.stats", "stats")
