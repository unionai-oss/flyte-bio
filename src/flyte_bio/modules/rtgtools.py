"""rtgtools — RTG Tools.

Exposes :func:`rtgtools_format`, which converts sequences into RTG's SDF
format (``rtg format``): a FASTA reference (as ``rtg vcfeval`` takes), FASTQ
reads or SAM/BAM reads, as ``args``' ``--format`` says (FASTA by default).
With ``single_end`` (the default, as upstream) ``input1`` is formatted on its
own; otherwise ``input1`` and ``input2`` are the left and right reads.
``sam_rg`` is a file holding the read group header line to record
(``--sam-rg``).

The SDF is upstream's ``<prefix>.sdf`` directory.

It also benchmarks variant calls against a truth set:

- :func:`rtgtools_vcfeval` — small variants, compared haplotype-aware against
  a reference SDF (``rtg vcfeval``), optionally within ``regions_bed``
  (``--bed-regions``) and ``truth_bed`` (``--evaluation-regions``).
- :func:`rtgtools_bndeval` — breakends (``rtg bndeval``), optionally within
  ``regions_bed``.
- :func:`rtgtools_cnveval` — copy-number variants within evaluation regions
  (``rtg cnveval``).
- :func:`rtgtools_svdecompose` — split structural variants into breakends
  (``rtg svdecompose``) for ``bndeval``: ``<prefix>.vcf.gz`` and its
  ``.tbi``.

As upstream, each evaluation's outputs are renamed ``<prefix>.<name>``
(``<prefix>.tp.vcf.gz``, ``<prefix>.summary.txt`` …), and VCFs given without
an index are indexed first (``rtg index``). RTG records its command line in
its VCF and ROC headers, so every input is linked into the working directory
and named relative to it, as upstream; ``vcfeval`` and ``cnveval`` record
their ``--threads`` (the task's CPUs) too. :func:`vcfeval_result` and
:func:`cnveval_result` read the outputs of a task run with other resources
(``rtgtools_vcfeval_cmd.override(resources=...)``).

RTG's memory (``RTG_MEM``) is 80% of the task's memory (upstream gives it all
of it, which leaves the JVM no headroom in a container).
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import JAVA_HEAP_MB, TASK_LIMITS

RTGTOOLS_IMAGE = "community.wave.seqera.io/library/rtg-tools:3.13--3465421f1b0be0ce"
# Upstream's process_low.
DEFAULT_RESOURCES = flyte.Resources(cpu=2, memory="4Gi")
# Upstream's process_medium, scaled down: the evaluations are light on test data.
EVAL_RESOURCES = flyte.Resources(cpu=4, memory="8Gi")

# Links the query and truth VCFs (and indexes, when given) into the working directory, indexing any VCF without one,
# as upstream. Run after TASK_LIMITS and `cd "$W"`, with QUERY, QUERY_TBI, TRUTH and TRUTH_TBI set.
LINK_AND_INDEX_VCFS = r"""
        for f in "${QUERY[0]}" "${QUERY_TBI[@]}" "${TRUTH[0]}" "${TRUTH_TBI[@]}"; do
            ln -s "$f" "$(basename "$f")"
        done
        QUERY_VCF=$(basename "${QUERY[0]}")
        TRUTH_VCF=$(basename "${TRUTH[0]}")
        if [ ${#QUERY_TBI[@]} -eq 0 ]; then rtg index "$QUERY_VCF" || exit 1; fi
        if [ ${#TRUTH_TBI[@]} -eq 0 ]; then rtg index "$TRUTH_VCF" || exit 1; fi
""".strip()

# As upstream: the evaluation's files are renamed <prefix>.<name>, but for its done and progress markers. Run with P,
# OUT and EVAL_DIR (the evaluation's output directory) set.
RENAME_OUTPUTS = r"""
        rm -f "$EVAL_DIR/done" "$EVAL_DIR/progress"
        for f in "$EVAL_DIR"/*; do mv "$f" "$OUT/$P.${f##*/}"; done
""".strip()

rtgtools_format_cmd = shell.create(
    name="rtgtools_format",
    image=RTGTOOLS_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "input1": File,
        "input2": File | None,
        "sam_rg": File | None,
        "single_end": bool,
        "prefix": str,
        "args": str,
    },
    defaults={"single_end": True, "args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        IN1=({{inputs.input1}})
        shopt -s nullglob; IN2=({{inputs.input2}}); RG=({{inputs.sam_rg}}); shopt -u nullglob
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"

        if [ {{inputs.single_end}} = true ]; then
            INPUT=("${{IN1[0]}}")
        else
            if [ ${{#IN2[@]}} -eq 0 ]; then
                echo "paired-end formatting needs input2" >&2
                exit 1
            fi
            INPUT=(--left "${{IN1[0]}}" --right "${{IN2[0]}}")
        fi
        RG_OPTION=()
        if [ ${{#RG[@]}} -gt 0 ]; then RG_OPTION=(--sam-rg "${{RG[0]}}"); fi

        rtg RTG_MEM={JAVA_HEAP_MB}M format \
            $ARGS \
            "${{RG_OPTION[@]}}" \
            --output "$P.sdf" \
            "${{INPUT[@]}}"
        mv "$P.sdf" "$OUT"/
    """,
)

rtgtools_vcfeval_cmd = shell.create(
    name="rtgtools_vcfeval",
    image=RTGTOOLS_IMAGE,
    resources=EVAL_RESOURCES,
    inputs={
        "query_vcf": File,
        "query_tbi": File | None,
        "truth_vcf": File,
        "truth_tbi": File | None,
        "truth_bed": File | None,
        "regions_bed": File | None,
        "sdf": Dir,
        "sdf_name": str,
        "prefix": str,
        "args": str,
    },
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        QUERY=({{inputs.query_vcf}}); TRUTH=({{inputs.truth_vcf}})
        shopt -s nullglob
        QUERY_TBI=({{inputs.query_tbi}}); TRUTH_TBI=({{inputs.truth_tbi}})
        TRUTH_BED=({{inputs.truth_bed}}); REGIONS=({{inputs.regions_bed}})
        shopt -u nullglob
        SDF={{inputs.sdf_name}}
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        ln -s {{inputs.sdf}} "$SDF"
        {LINK_AND_INDEX_VCFS}
        OPTIONAL=()
        if [ ${{#REGIONS[@]}} -gt 0 ]; then
            ln -s "${{REGIONS[0]}}" "$(basename "${{REGIONS[0]}}")"
            OPTIONAL+=(--bed-regions="$(basename "${{REGIONS[0]}}")")
        fi
        if [ ${{#TRUTH_BED[@]}} -gt 0 ]; then
            ln -s "${{TRUTH_BED[0]}}" "$(basename "${{TRUTH_BED[0]}}")"
            OPTIONAL+=(--evaluation-regions="$(basename "${{TRUTH_BED[0]}}")")
        fi

        rtg RTG_MEM={JAVA_HEAP_MB}M vcfeval \
            $ARGS \
            --baseline="$TRUTH_VCF" \
            "${{OPTIONAL[@]}}" \
            --calls="$QUERY_VCF" \
            --output=output \
            --template="$SDF" \
            --threads=$TASK_CPUS || exit 1
        EVAL_DIR=output
        {RENAME_OUTPUTS}
    """,
)

rtgtools_bndeval_cmd = shell.create(
    name="rtgtools_bndeval",
    image=RTGTOOLS_IMAGE,
    resources=EVAL_RESOURCES,
    inputs={
        "query_vcf": File,
        "query_tbi": File | None,
        "truth_vcf": File,
        "truth_tbi": File | None,
        "regions_bed": File | None,
        "prefix": str,
        "args": str,
    },
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        QUERY=({{inputs.query_vcf}}); TRUTH=({{inputs.truth_vcf}})
        shopt -s nullglob
        QUERY_TBI=({{inputs.query_tbi}}); TRUTH_TBI=({{inputs.truth_tbi}}); REGIONS=({{inputs.regions_bed}})
        shopt -u nullglob
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        {LINK_AND_INDEX_VCFS}
        OPTIONAL=()
        if [ ${{#REGIONS[@]}} -gt 0 ]; then
            ln -s "${{REGIONS[0]}}" "$(basename "${{REGIONS[0]}}")"
            OPTIONAL+=(--bed-regions="$(basename "${{REGIONS[0]}}")")
        fi

        rtg RTG_MEM={JAVA_HEAP_MB}M bndeval \
            $ARGS \
            --baseline="$TRUTH_VCF" \
            "${{OPTIONAL[@]}}" \
            --calls="$QUERY_VCF" \
            --output="$P" || exit 1
        EVAL_DIR="$P"
        {RENAME_OUTPUTS}
    """,
)

rtgtools_cnveval_cmd = shell.create(
    name="rtgtools_cnveval",
    image=RTGTOOLS_IMAGE,
    resources=EVAL_RESOURCES,
    inputs={
        "query_vcf": File,
        "query_tbi": File | None,
        "truth_vcf": File,
        "truth_tbi": File | None,
        "evaluation_regions_bed": File,
        "prefix": str,
        "args": str,
    },
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        QUERY=({{inputs.query_vcf}}); TRUTH=({{inputs.truth_vcf}}); REGIONS=({{inputs.evaluation_regions_bed}})
        shopt -s nullglob; QUERY_TBI=({{inputs.query_tbi}}); TRUTH_TBI=({{inputs.truth_tbi}}); shopt -u nullglob
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        {LINK_AND_INDEX_VCFS}
        ln -s "${{REGIONS[0]}}" "$(basename "${{REGIONS[0]}}")"

        rtg RTG_MEM={JAVA_HEAP_MB}M cnveval \
            $ARGS \
            --baseline="$TRUTH_VCF" \
            --calls="$QUERY_VCF" \
            --evaluation-regions="$(basename "${{REGIONS[0]}}")" \
            --output="$P" \
            --threads=$TASK_CPUS || exit 1
        EVAL_DIR="$P"
        {RENAME_OUTPUTS}
    """,
)

rtgtools_svdecompose_cmd = shell.create(
    name="rtgtools_svdecompose",
    image=RTGTOOLS_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"vcf": File, "tbi": File | None, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        VCF=({{inputs.vcf}})
        shopt -s nullglob; TBI=({{inputs.tbi}}); shopt -u nullglob
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        for f in "${{VCF[0]}}" "${{TBI[@]}}"; do ln -s "$f" "$(basename "$f")"; done
        INPUT=$(basename "${{VCF[0]}}")
        if [ ${{#TBI[@]}} -eq 0 ]; then rtg index "$INPUT" || exit 1; fi

        rtg RTG_MEM={JAVA_HEAP_MB}M svdecompose \
            $ARGS \
            --input="$INPUT" \
            --output="$P.vcf.gz" || exit 1
        rtg index "$P.vcf.gz" || exit 1
        mv "$P.vcf.gz" "$P.vcf.gz.tbi" "$OUT"/
    """,
)

env = flyte.TaskEnvironment.from_task(
    "rtgtools",
    rtgtools_format_cmd.as_task(),
    rtgtools_vcfeval_cmd.as_task(),
    rtgtools_bndeval_cmd.as_task(),
    rtgtools_cnveval_cmd.as_task(),
    rtgtools_svdecompose_cmd.as_task(),
)


async def rtgtools_format(
    input1: File,
    input2: File | None = None,
    sam_rg: File | None = None,
    single_end: bool = True,
    prefix: str = "sample",
    args: str = "",
) -> Dir:
    """``input1`` (and, paired, ``input2``) in RTG's SDF format: the ``<prefix>.sdf`` Dir."""
    results = await rtgtools_format_cmd(
        input1=input1, input2=input2, sam_rg=sam_rg, single_end=single_end, prefix=prefix, args=args
    )
    return Dir.from_existing_remote(results.path.rstrip("/") + f"/{prefix}.sdf")


async def pick(results: Dir, name: str, tool: str) -> File:
    found = await results.get_file(name)
    if found is None:
        raise FileNotFoundError(f"rtg {tool} wrote no {name}")
    return found


@dataclass
class VcfevalResult:
    tp_vcf: File  # <prefix>.tp.vcf.gz: the calls that match the truth (+ tp_tbi)
    tp_tbi: File
    fn_vcf: File  # <prefix>.fn.vcf.gz: the truth's variants the calls miss (+ fn_tbi)
    fn_tbi: File
    fp_vcf: File  # <prefix>.fp.vcf.gz: the calls not in the truth (+ fp_tbi)
    fp_tbi: File
    baseline_vcf: File  # <prefix>.tp-baseline.vcf.gz: the truth's variants the calls match (+ baseline_tbi)
    baseline_tbi: File
    snp_roc: File  # <prefix>.snp_roc.tsv.gz
    non_snp_roc: File  # <prefix>.non_snp_roc.tsv.gz
    weighted_roc: File  # <prefix>.weighted_roc.tsv.gz
    summary: File  # <prefix>.summary.txt
    phasing: File  # <prefix>.phasing.txt
    results: Dir  # everything, with <prefix>.vcfeval.log


async def vcfeval_result(results: Dir, prefix: str) -> VcfevalResult:
    """The outputs ``rtgtools_vcfeval_cmd`` wrote to ``results``."""
    names = {
        "tp_vcf": "tp.vcf.gz",
        "fn_vcf": "fn.vcf.gz",
        "fp_vcf": "fp.vcf.gz",
        "baseline_vcf": "tp-baseline.vcf.gz",
        "snp_roc": "snp_roc.tsv.gz",
        "non_snp_roc": "non_snp_roc.tsv.gz",
        "weighted_roc": "weighted_roc.tsv.gz",
        "summary": "summary.txt",
        "phasing": "phasing.txt",
    }
    files = {key: await pick(results, f"{prefix}.{name}", "vcfeval") for key, name in names.items()}
    for vcf in ("tp", "fn", "fp", "baseline"):
        files[f"{vcf}_tbi"] = await pick(results, f"{prefix}.{names[vcf + '_vcf']}.tbi", "vcfeval")
    return VcfevalResult(**files, results=results)


async def rtgtools_vcfeval(
    query_vcf: File,
    truth_vcf: File,
    sdf: Dir,
    prefix: str,
    query_tbi: File | None = None,
    truth_tbi: File | None = None,
    truth_bed: File | None = None,
    regions_bed: File | None = None,
    sdf_name: str = "",
    args: str = "",
) -> VcfevalResult:
    """Compare ``query_vcf``'s small variants with ``truth_vcf``'s against the reference ``sdf``.

    ``regions_bed`` restricts both to its regions (``--bed-regions``);
    ``truth_bed`` holds the truth's confident regions (``--evaluation-regions``).
    The SDF is linked into the working directory as ``sdf_name`` (by default
    its own name), which vcfeval records in its outputs' command lines.
    """
    sdf_name = sdf_name or sdf.path.rstrip("/").rsplit("/", 1)[-1]
    results = await rtgtools_vcfeval_cmd(
        query_vcf=query_vcf,
        query_tbi=query_tbi,
        truth_vcf=truth_vcf,
        truth_tbi=truth_tbi,
        truth_bed=truth_bed,
        regions_bed=regions_bed,
        sdf=sdf,
        sdf_name=sdf_name,
        prefix=prefix,
        args=args,
    )
    return await vcfeval_result(results, prefix)


@dataclass
class BndevalResult:
    tp_vcf: File  # <prefix>.tp.vcf.gz (+ tp_tbi)
    tp_tbi: File
    fn_vcf: File  # <prefix>.fn.vcf.gz (+ fn_tbi)
    fn_tbi: File
    fp_vcf: File  # <prefix>.fp.vcf.gz (+ fp_tbi)
    fp_tbi: File
    baseline_vcf: File  # <prefix>.tp-baseline.vcf.gz (+ baseline_tbi)
    baseline_tbi: File
    weighted_roc: File  # <prefix>.weighted_roc.tsv.gz
    summary: File  # <prefix>.summary.txt
    results: Dir  # everything, with <prefix>.bndeval.log


async def rtgtools_bndeval(
    query_vcf: File,
    truth_vcf: File,
    prefix: str,
    query_tbi: File | None = None,
    truth_tbi: File | None = None,
    regions_bed: File | None = None,
    args: str = "",
) -> BndevalResult:
    """Compare ``query_vcf``'s breakends with ``truth_vcf``'s, within ``regions_bed`` when given."""
    results = await rtgtools_bndeval_cmd(
        query_vcf=query_vcf,
        query_tbi=query_tbi,
        truth_vcf=truth_vcf,
        truth_tbi=truth_tbi,
        regions_bed=regions_bed,
        prefix=prefix,
        args=args,
    )
    files = {}
    for key, name in (("tp", "tp"), ("fn", "fn"), ("fp", "fp"), ("baseline", "tp-baseline")):
        files[f"{key}_vcf"] = await pick(results, f"{prefix}.{name}.vcf.gz", "bndeval")
        files[f"{key}_tbi"] = await pick(results, f"{prefix}.{name}.vcf.gz.tbi", "bndeval")
    return BndevalResult(
        **files,
        weighted_roc=await pick(results, f"{prefix}.weighted_roc.tsv.gz", "bndeval"),
        summary=await pick(results, f"{prefix}.summary.txt", "bndeval"),
        results=results,
    )


@dataclass
class CnvevalResult:
    baseline_bed: File  # <prefix>.baseline.bed.gz: the truth's regions, annotated
    calls_bed: File  # <prefix>.calls.bed.gz: the calls' regions, annotated
    weighted_roc: File  # <prefix>.weighted_roc.tsv.gz
    summary: File  # <prefix>.summary.txt
    results: Dir  # everything, with <prefix>.cnveval.log


async def cnveval_result(results: Dir, prefix: str) -> CnvevalResult:
    """The outputs ``rtgtools_cnveval_cmd`` wrote to ``results``."""
    return CnvevalResult(
        baseline_bed=await pick(results, f"{prefix}.baseline.bed.gz", "cnveval"),
        calls_bed=await pick(results, f"{prefix}.calls.bed.gz", "cnveval"),
        weighted_roc=await pick(results, f"{prefix}.weighted_roc.tsv.gz", "cnveval"),
        summary=await pick(results, f"{prefix}.summary.txt", "cnveval"),
        results=results,
    )


async def rtgtools_cnveval(
    query_vcf: File,
    truth_vcf: File,
    evaluation_regions_bed: File,
    prefix: str,
    query_tbi: File | None = None,
    truth_tbi: File | None = None,
    args: str = "",
) -> CnvevalResult:
    """Compare ``query_vcf``'s copy-number variants with ``truth_vcf``'s over ``evaluation_regions_bed``."""
    results = await rtgtools_cnveval_cmd(
        query_vcf=query_vcf,
        query_tbi=query_tbi,
        truth_vcf=truth_vcf,
        truth_tbi=truth_tbi,
        evaluation_regions_bed=evaluation_regions_bed,
        prefix=prefix,
        args=args,
    )
    return await cnveval_result(results, prefix)


@dataclass
class SvdecomposeResult:
    vcf: File  # <prefix>.vcf.gz
    tbi: File  # <prefix>.vcf.gz.tbi


async def rtgtools_svdecompose(vcf: File, prefix: str, tbi: File | None = None, args: str = "") -> SvdecomposeResult:
    """``vcf``'s structural variants split into breakends: ``<prefix>.vcf.gz`` and its ``.tbi``."""
    if vcf.name == f"{prefix}.vcf.gz":
        raise ValueError("input and output names are the same; choose another prefix")
    results = await rtgtools_svdecompose_cmd(vcf=vcf, tbi=tbi, prefix=prefix, args=args)
    return SvdecomposeResult(
        vcf=await pick(results, f"{prefix}.vcf.gz", "svdecompose"),
        tbi=await pick(results, f"{prefix}.vcf.gz.tbi", "svdecompose"),
    )
