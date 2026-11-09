"""glimpse2 — GLIMPSE2 genotype imputation and phasing of low-coverage data.

Exposes the steps of GLIMPSE 2.0.1:

- :func:`glimpse2_chunk` — split a region of a reference panel into
  imputation chunks with the ``recursive`` or ``sequential`` model,
  ``<prefix>.txt``.
- :func:`glimpse2_phase` — impute and phase samples in one chunk, from their
  reads (one or more BAM/CRAM files, or a BAM list) or their genotype
  likelihoods (one VCF/BCF), into ``<prefix>.<suffix>``; with reads it also
  writes ``<prefix>_stats_coverage.txt.gz``. Upstream's default prefix is
  ``<id>_<output region>`` (``:`` turned into ``_``), or ``<id>_<reference
  panel name>`` without regions.
- :func:`glimpse2_ligate` — ligate the phased chunks (each with its index)
  into ``<prefix>.<suffix>``. As upstream, the chunks are listed in name order.
- :func:`glimpse2_concordance` — the concordance of imputed genotypes with
  truth genotypes, by allele frequency bins (``<prefix>.error.cal.txt.gz``,
  ``.error.grp.txt.gz``, ``.error.spl.txt.gz``, ``.rsquare.grp.txt.gz``,
  ``.rsquare.spl.txt.gz`` and, with ``--out-r2-per-site``,
  ``<prefix>_r2_sites.txt.gz``).

Inputs are linked into a working directory under their own names, indexes
beside their files. ``--thread``/``--threads`` follows the task's CPUs.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

GLIMPSE2_IMAGE = "quay.io/biocontainers/glimpse-bio:2.0.1--h46b9e50_1"
# Upstream's process_single for chunk and phase (phase's output depends on
# its threads, and upstream's tests run it on 1 CPU), process_low for ligate
# and concordance.
SINGLE_RESOURCES = flyte.Resources(cpu=1, memory="6Gi")
LOW_RESOURCES = flyte.Resources(cpu=2, memory="4Gi")

# `link FILE` links FILE into the working directory under its own name and
# prints that name; `link_index FILE INDEX` links INDEX beside FILE under its
# own suffix (.csi, .tbi, .bai, .crai).
LINK_HELPERS = r"""
        link() { ln -s "$1" "$(basename "$1")"; basename "$1"; }
        link_index() { ln -s "$2" "$1.${2##*.}"; }
        W=$(mktemp -d); cd "$W"
"""

glimpse2_chunk_cmd = shell.create(
    name="glimpse2_chunk",
    image=GLIMPSE2_IMAGE,
    resources=SINGLE_RESOURCES,
    inputs={
        "vcf": File,
        "index": File,
        "region": str,
        "genetic_map": File | None,
        "model": str,
        "prefix": str,
        "args": str,
    },
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        V=({{inputs.vcf}}); I=({{inputs.index}})
        shopt -s nullglob; M=({{inputs.genetic_map}}); shopt -u nullglob
        REGION={{inputs.region}}; MODEL={{inputs.model}}; P={{inputs.prefix}}
        ARGS={{inputs.args}}
        {LINK_HELPERS}
        NAME=$(link "${{V[0]}}"); link_index "$NAME" "${{I[0]}}"
        MAP=()
        if [ ${{#M[@]}} -gt 0 ]; then MAP=(--map "$(link "${{M[0]}}")"); fi
        GLIMPSE2_chunk $ARGS "${{MAP[@]}}" --"$MODEL" --input "$NAME" --region "$REGION" \
            --threads $TASK_CPUS --output "$P.txt" || exit 1
        mv "$P.txt" {{outputs.results}}/
    """,
)

# As upstream: a BAM list wins; several BAM/CRAMs are listed (all_bam.txt) in
# name order; else one BAM/CRAM (--bam-file) or one VCF/BCF of genotype
# likelihoods (--input-gl).
glimpse2_phase_cmd = shell.create(
    name="glimpse2_phase",
    image=GLIMPSE2_IMAGE,
    resources=SINGLE_RESOURCES,
    inputs={
        "inputs": list[File],
        "indexes": list[File],
        "input_type": str,
        "bamlist": File | None,
        "samples_file": File | None,
        "input_region": str,
        "output_region": str,
        "reference": File,
        "reference_index": File,
        "genetic_map": File | None,
        "fasta": File | None,
        "fasta_index": File | None,
        "prefix": str,
        "suffix": str,
        "args": str,
    },
    defaults={"input_region": "", "output_region": "", "args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        REF=({{inputs.reference}}); REF_I=({{inputs.reference_index}})
        shopt -s nullglob
        IN=({{inputs.inputs}}); IDX=({{inputs.indexes}}); BAMLIST=({{inputs.bamlist}}); S=({{inputs.samples_file}})
        M=({{inputs.genetic_map}}); FA=({{inputs.fasta}}); FAI=({{inputs.fasta_index}})
        shopt -u nullglob
        INPUT_TYPE={{inputs.input_type}}
        INPUT_REGION={{inputs.input_region}}; OUTPUT_REGION={{inputs.output_region}}
        P={{inputs.prefix}}; SUFFIX={{inputs.suffix}}
        ARGS={{inputs.args}}
        {LINK_HELPERS}
        NAMES=()
        for F in "${{IN[@]}}" ; do NAMES+=("$(link "$F")"); done
        for F in "${{IDX[@]}}"; do link "$F" > /dev/null; done
        REF_NAME=$(link "${{REF[0]}}"); link_index "$REF_NAME" "${{REF_I[0]}}"
        if [ ${{#BAMLIST[@]}} -gt 0 ]; then
            INPUT=(--bam-list "$(link "${{BAMLIST[0]}}")")
        elif [ ${{#NAMES[@]}} -gt 1 ]; then
            ls -1 | grep '\.cram$\|\.bam$' | sort > all_bam.txt
            INPUT=(--bam-list all_bam.txt)
        elif [ "$INPUT_TYPE" = bam ]; then
            INPUT=(--bam-file "${{NAMES[0]}}")
        else
            INPUT=(--input-gl "${{NAMES[0]}}")
        fi
        EXTRA=()
        if [ ${{#M[@]}} -gt 0 ]; then EXTRA+=(--map "$(link "${{M[0]}}")"); fi
        if [ ${{#FA[@]}} -gt 0 ]; then
            FA_NAME=$(link "${{FA[0]}}"); EXTRA+=(--fasta "$FA_NAME")
            if [ ${{#FAI[@]}} -gt 0 ]; then ln -s "${{FAI[0]}}" "$FA_NAME.fai"; fi
        fi
        if [ ${{#S[@]}} -gt 0 ]; then EXTRA+=(--samples-file "$(link "${{S[0]}}")"); fi
        if [ -n "$INPUT_REGION" ]; then EXTRA+=(--input-region "$INPUT_REGION"); fi
        if [ -n "$OUTPUT_REGION" ]; then EXTRA+=(--output-region "$OUTPUT_REGION"); fi
        GLIMPSE2_phase $ARGS "${{INPUT[@]}}" --reference "$REF_NAME" "${{EXTRA[@]}}" \
            --thread $TASK_CPUS --output "$P.$SUFFIX" || exit 1
        mv "$P.$SUFFIX" {{outputs.results}}/
        if [ -f "$P"_stats_coverage.txt.gz ]; then mv "$P"_stats_coverage.txt.gz {{outputs.results}}/; fi
    """,
)

# As upstream, the chunks are listed (all_files.txt) in name order.
glimpse2_ligate_cmd = shell.create(
    name="glimpse2_ligate",
    image=GLIMPSE2_IMAGE,
    resources=LOW_RESOURCES,
    inputs={"vcfs": list[File], "indexes": list[File], "prefix": str, "suffix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        shopt -s nullglob; VCFS=({{inputs.vcfs}}); IDX=({{inputs.indexes}}); shopt -u nullglob
        P={{inputs.prefix}}; SUFFIX={{inputs.suffix}}
        ARGS={{inputs.args}}
        {LINK_HELPERS}
        for F in "${{VCFS[@]}}" "${{IDX[@]}}"; do link "$F" > /dev/null; done
        for F in "${{VCFS[@]}}"; do basename "$F"; done | LC_ALL=C sort > all_files.txt
        GLIMPSE2_ligate $ARGS --input all_files.txt --thread $TASK_CPUS --output "$P.$SUFFIX" || exit 1
        mv "$P.$SUFFIX" {{outputs.results}}/
    """,
)

# As upstream: one input line per region, "<region> <freq> <truth> <estimate>".
# bins, ac_bins and allele_counts are whitespace-separated lists, passed as
# separate arguments.
glimpse2_concordance_cmd = shell.create(
    name="glimpse2_concordance",
    image=GLIMPSE2_IMAGE,
    resources=LOW_RESOURCES,
    inputs={
        "estimate": File,
        "estimate_index": File,
        "truth": File,
        "truth_index": File,
        "freq": File,
        "freq_index": File,
        "samples": File | None,
        "groups": File | None,
        "regions": str,
        "bins": str,
        "ac_bins": str,
        "allele_counts": str,
        "min_val_gl": str,
        "min_val_dp": str,
        "prefix": str,
        "args": str,
    },
    defaults={"bins": "", "ac_bins": "", "allele_counts": "", "min_val_gl": "", "min_val_dp": "", "args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        E=({{inputs.estimate}}); EI=({{inputs.estimate_index}}); T=({{inputs.truth}}); TI=({{inputs.truth_index}})
        F=({{inputs.freq}}); FI=({{inputs.freq_index}})
        shopt -s nullglob; S=({{inputs.samples}}); G=({{inputs.groups}}); shopt -u nullglob
        REGIONS={{inputs.regions}}
        BINS={{inputs.bins}}; AC_BINS={{inputs.ac_bins}}; ALLELE_COUNTS={{inputs.allele_counts}}
        MIN_VAL_GL={{inputs.min_val_gl}}; MIN_VAL_DP={{inputs.min_val_dp}}
        P={{inputs.prefix}}
        set -f
        ARGS={{inputs.args}}
        {LINK_HELPERS}
        ESTIMATE=$(link "${{E[0]}}"); link_index "$ESTIMATE" "${{EI[0]}}"
        TRUTH=$(link "${{T[0]}}"); link_index "$TRUTH" "${{TI[0]}}"
        FREQ=$(link "${{F[0]}}"); link_index "$FREQ" "${{FI[0]}}"
        EXTRA=()
        if [ ${{#S[@]}} -gt 0 ]; then EXTRA+=(--samples "$(link "${{S[0]}}")"); fi
        if [ ${{#G[@]}} -gt 0 ]; then EXTRA+=(--groups "$(link "${{G[0]}}")"); fi
        if [ -n "$BINS" ]; then EXTRA+=(--bins $BINS); fi
        if [ -n "$AC_BINS" ]; then EXTRA+=(--ac-bins $AC_BINS); fi
        if [ -n "$ALLELE_COUNTS" ]; then EXTRA+=(--allele-counts $ALLELE_COUNTS); fi
        if [ -n "$MIN_VAL_GL" ]; then EXTRA+=(--min-val-gl "$MIN_VAL_GL"); fi
        if [ -n "$MIN_VAL_DP" ]; then EXTRA+=(--min-val-dp "$MIN_VAL_DP"); fi
        printf '%s' "$REGIONS" > regions.txt
        sed "s/\$/ $FREQ $TRUTH $ESTIMATE/" regions.txt > input.txt
        GLIMPSE2_concordance $ARGS "${{EXTRA[@]}}" --input input.txt --thread $TASK_CPUS --output "$P" || exit 1
        set +f
        mv "$P".*.txt.gz {{outputs.results}}/
        if [ -f "$P"_r2_sites.txt.gz ]; then mv "$P"_r2_sites.txt.gz {{outputs.results}}/; fi
    """,
)

env = flyte.TaskEnvironment.from_task(
    "glimpse2",
    glimpse2_chunk_cmd.as_task(),
    glimpse2_phase_cmd.as_task(),
    glimpse2_ligate_cmd.as_task(),
    glimpse2_concordance_cmd.as_task(),
)


async def pick(results: Dir, name: str, tool: str) -> File:
    found = await results.get_file(name)
    if found is None:
        raise FileNotFoundError(f"{tool} wrote no {name}")
    return found


async def glimpse2_chunk(
    vcf: File,
    index: File,
    region: str,
    prefix: str,
    model: str = "recursive",
    genetic_map: File | None = None,
    args: str = "",
) -> File:
    """``<prefix>.txt``: the imputation chunks of ``region`` of the panel ``vcf``.

    ``model`` is ``recursive`` or ``sequential`` (``--recursive``,
    ``--sequential``); ``genetic_map`` is passed as ``--map``.
    """
    results = await glimpse2_chunk_cmd(
        vcf=vcf, index=index, region=region, genetic_map=genetic_map, model=model, prefix=prefix, args=args
    )
    return await pick(results, f"{prefix}.txt", "GLIMPSE2_chunk")


def glimpse2_input_type(name: str) -> str:
    """``bam`` for reads, ``gl`` for genotype likelihoods, else the extension (as upstream)."""
    if name.endswith(("cram", "bam")):
        return "bam"
    if name.endswith(("vcf", "bcf", "vcf.gz")):
        return "gl"
    return name.rsplit(".", 1)[-1]


@dataclass
class Glimpse2PhaseResult:
    vcf: File  # <prefix>.<suffix>
    stats_coverage: File | None  # <prefix>_stats_coverage.txt.gz, with reads as input
    results: Dir


async def glimpse2_phase(
    inputs: list[File],
    indexes: list[File],
    reference: File,
    reference_index: File,
    prefix: str,
    input_region: str = "",
    output_region: str = "",
    bamlist: File | None = None,
    samples_file: File | None = None,
    genetic_map: File | None = None,
    fasta: File | None = None,
    fasta_index: File | None = None,
    suffix: str = "vcf.gz",
    args: str = "",
) -> Glimpse2PhaseResult:
    """Impute and phase ``inputs`` against the panel ``reference``: ``<prefix>.<suffix>``.

    ``inputs`` are BAM/CRAM files (with their ``indexes``) or one VCF/BCF of
    genotype likelihoods. ``bamlist`` lists the BAMs by name, optionally with
    a sample name each. ``fasta`` (with its ``fasta_index``) decodes CRAMs.
    ``input_region`` and ``output_region`` come from :func:`glimpse2_chunk`.
    """
    if suffix not in ("vcf", "vcf.gz", "bcf", "bgen"):
        raise ValueError(f"Output suffix must be one of vcf, vcf.gz, bcf, bgen. Found: {suffix}")
    types = list(dict.fromkeys(glimpse2_input_type(f.name) for f in inputs))
    if len(types) != 1 or types[0] not in ("bam", "gl"):
        raise ValueError(
            "Input files must be of the same type and either .bam/.cram or .vcf/.vcf.gz/.bcf format. "
            f"Found: [{', '.join(types)}]"
        )
    if types[0] == "gl" and len(inputs) > 1:
        raise ValueError("Only one input .vcf/.vcf.gz/.bcf file can be provided")
    results = await glimpse2_phase_cmd(
        inputs=inputs,
        indexes=indexes,
        input_type=types[0],
        bamlist=bamlist,
        samples_file=samples_file,
        input_region=input_region,
        output_region=output_region,
        reference=reference,
        reference_index=reference_index,
        genetic_map=genetic_map,
        fasta=fasta,
        fasta_index=fasta_index,
        prefix=prefix,
        suffix=suffix,
        args=args,
    )
    vcf = await pick(results, f"{prefix}.{suffix}", "GLIMPSE2_phase")
    stats = await results.get_file(f"{prefix}_stats_coverage.txt.gz")
    return Glimpse2PhaseResult(vcf=vcf, stats_coverage=stats, results=results)


async def glimpse2_ligate(
    vcfs: list[File], indexes: list[File], prefix: str, suffix: str = "vcf.gz", args: str = ""
) -> File:
    """Ligate phased chunks ``vcfs`` (with their ``indexes``, named after them) into ``<prefix>.<suffix>``."""
    names = [f.name for f in vcfs]
    if len(set(names)) != len(names):
        raise ValueError(f"GLIMPSE2_ligate inputs need distinct names, got {names}")
    results = await glimpse2_ligate_cmd(vcfs=vcfs, indexes=indexes, prefix=prefix, suffix=suffix, args=args)
    return await pick(results, f"{prefix}.{suffix}", "GLIMPSE2_ligate")


@dataclass
class Glimpse2ConcordanceResult:
    errors_cal: File  # <prefix>.error.cal.txt.gz
    errors_grp: File  # <prefix>.error.grp.txt.gz
    errors_spl: File  # <prefix>.error.spl.txt.gz
    rsquare_grp: File  # <prefix>.rsquare.grp.txt.gz
    rsquare_spl: File  # <prefix>.rsquare.spl.txt.gz
    rsquare_per_site: File | None  # <prefix>_r2_sites.txt.gz, with --out-r2-per-site
    results: Dir


async def glimpse2_concordance(
    estimate: File,
    estimate_index: File,
    truth: File,
    truth_index: File,
    freq: File,
    freq_index: File,
    regions: list[str],
    prefix: str,
    samples: File | None = None,
    groups: File | None = None,
    bins: str = "",
    ac_bins: str = "",
    allele_counts: str = "",
    min_val_gl: str = "",
    min_val_dp: str = "",
    args: str = "",
) -> Glimpse2ConcordanceResult:
    """The concordance of the imputed ``estimate`` with ``truth`` in ``regions``, by ``freq``'s allele frequencies.

    Exactly one of ``groups`` (a file), ``bins``, ``ac_bins`` and
    ``allele_counts`` (whitespace-separated values, e.g. ``"0 0.01 0.05"``)
    sets the frequency bins, as upstream.
    """
    if sum(bool(x) for x in (groups, bins, ac_bins, allele_counts)) != 1:
        raise ValueError("One and only one argument should be selected between groups, bins, ac_bins, allele_counts")
    if "--gt-val" in args and (min_val_gl or min_val_dp):
        raise ValueError("If --gt-val is set, --min-val-gl nor --min-val-dp must be set")
    results = await glimpse2_concordance_cmd(
        estimate=estimate,
        estimate_index=estimate_index,
        truth=truth,
        truth_index=truth_index,
        freq=freq,
        freq_index=freq_index,
        samples=samples,
        groups=groups,
        regions="\n".join(regions),
        bins=bins,
        ac_bins=ac_bins,
        allele_counts=allele_counts,
        min_val_gl=min_val_gl,
        min_val_dp=min_val_dp,
        prefix=prefix,
        args=args,
    )
    tool = "GLIMPSE2_concordance"
    return Glimpse2ConcordanceResult(
        errors_cal=await pick(results, f"{prefix}.error.cal.txt.gz", tool),
        errors_grp=await pick(results, f"{prefix}.error.grp.txt.gz", tool),
        errors_spl=await pick(results, f"{prefix}.error.spl.txt.gz", tool),
        rsquare_grp=await pick(results, f"{prefix}.rsquare.grp.txt.gz", tool),
        rsquare_spl=await pick(results, f"{prefix}.rsquare.spl.txt.gz", tool),
        rsquare_per_site=await results.get_file(f"{prefix}_r2_sites.txt.gz"),
        results=results,
    )
