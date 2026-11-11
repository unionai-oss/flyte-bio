"""happy — hap.py, Illumina's small variant benchmarking, and its reports.

Exposes:

- :func:`happy_happy` — compare a query VCF's small variants with a truth
  VCF's (``hap.py``) against a reference, optionally within confident regions
  (``-R``) and targets (``-T``), with a false-positive BED
  (``--false-positives``) or stratified by the regions a TSV names
  (``--stratification``; the TSV's BED paths are relative to the working
  directory, where the ``stratification_beds`` are linked). Outputs keep
  hap.py's names: ``<prefix>.summary.csv``, ``<prefix>.extended.csv``, the
  ROC tables ``<prefix>.roc.all.csv.gz`` and
  ``<prefix>.roc.Locations.{SNP,INDEL}[.PASS].csv.gz``,
  ``<prefix>.runinfo.json``, ``<prefix>.metrics.json.gz`` and the annotated
  ``<prefix>.vcf.gz`` with its ``.tbi``. hap.py gets the task's CPUs as
  threads.
- :func:`happy_sompy` — compare somatic calls with a truth VCF (``som.py``),
  optionally within regions and targets, with a false-positive BED,
  ambiguous-region BEDs (``--ambiguous``) and BAMs for features
  (``--bam``): ``<prefix>.stats.csv``, ``<prefix>.metrics.json`` and, when
  ``args`` asks for a feature table, ``<prefix>.features.csv``.
- :func:`happy_report` — an HTML report of hap.py ROC tables
  (``rep.py``), each labelled by its method: ``<prefix>.html``.

hap.py and som.py run in upstream's hap.py image, the report in its own. The
reference is linked beside its ``.fai`` and every input is named relative to
the working directory, as upstream stages them.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

HAPPY_IMAGE = "quay.io/biocontainers/hap.py:0.3.15--py27hcb73b3d_0"
HAPPY_REPORT_IMAGE = "community.wave.seqera.io/library/pip_ga4gh-happy-report:72e6bf058d448a4c"
# Upstream's process_medium, scaled down: the comparisons are light on test data.
DEFAULT_RESOURCES = flyte.Resources(cpu=4, memory="8Gi")
# Upstream's process_single.
REPORT_RESOURCES = flyte.Resources(cpu=1, memory="2Gi")

# Links every staged input into the working directory, under its own name. Run after `cd "$W"`, with LINKED set to
# the inputs' paths.
LINK_INPUTS = r"""
        for f in "${LINKED[@]}"; do ln -s "$f" "$(basename "$f")"; done
""".strip()

happy_happy_cmd = shell.create(
    name="happy_happy",
    image=HAPPY_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "query_vcf": File,
        "truth_vcf": File,
        "regions_bed": File | None,
        "targets_bed": File | None,
        "fasta": File,
        "fai": File,
        "false_positives_bed": File | None,
        "stratification_tsv": File | None,
        "stratification_beds": list[File],
        "prefix": str,
        "args": str,
    },
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        QUERY=({{inputs.query_vcf}}); TRUTH=({{inputs.truth_vcf}}); FA=({{inputs.fasta}}); FAI=({{inputs.fai}})
        shopt -s nullglob
        REGIONS=({{inputs.regions_bed}}); TARGETS=({{inputs.targets_bed}})
        FALSE_POSITIVES=({{inputs.false_positives_bed}})
        STRAT_TSV=({{inputs.stratification_tsv}}); STRAT_BEDS=({{inputs.stratification_beds}})
        shopt -u nullglob
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        W=$(mktemp -d); cd "$W"
        LINKED=("${{QUERY[0]}}" "${{TRUTH[0]}}" "${{FA[0]}}" "${{FAI[0]}}" "${{REGIONS[@]}}" "${{TARGETS[@]}}"
            "${{FALSE_POSITIVES[@]}}" "${{STRAT_TSV[@]}}" "${{STRAT_BEDS[@]}}")
        {LINK_INPUTS}
        OPTIONAL=()
        if [ ${{#REGIONS[@]}} -gt 0 ]; then OPTIONAL+=(-R "$(basename "${{REGIONS[0]}}")"); fi
        if [ ${{#TARGETS[@]}} -gt 0 ]; then OPTIONAL+=(-T "$(basename "${{TARGETS[0]}}")"); fi
        if [ ${{#FALSE_POSITIVES[@]}} -gt 0 ]; then
            OPTIONAL+=(--false-positives "$(basename "${{FALSE_POSITIVES[0]}}")")
        fi
        if [ ${{#STRAT_TSV[@]}} -gt 0 ] && [ ${{#STRAT_BEDS[@]}} -gt 0 ]; then
            OPTIONAL+=(--stratification "$(basename "${{STRAT_TSV[0]}}")")
        fi
        hap.py \
            "$(basename "${{TRUTH[0]}}")" \
            "$(basename "${{QUERY[0]}}")" \
            $ARGS \
            --reference "$(basename "${{FA[0]}}")" \
            --threads $TASK_CPUS \
            "${{OPTIONAL[@]}}" \
            -o "$P" || exit 1
        mv "$P".* {{outputs.results}}/
    """,
)

happy_sompy_cmd = shell.create(
    name="happy_sompy",
    image=HAPPY_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "query_vcf": File,
        "truth_vcf": File,
        "regions_bed": File | None,
        "targets_bed": File | None,
        "fasta": File,
        "fai": File,
        "false_positives_bed": File | None,
        "ambiguous_beds": list[File],
        "bams": list[File],
        "prefix": str,
        "args": str,
    },
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        QUERY=({{inputs.query_vcf}}); TRUTH=({{inputs.truth_vcf}}); FA=({{inputs.fasta}}); FAI=({{inputs.fai}})
        shopt -s nullglob
        REGIONS=({{inputs.regions_bed}}); TARGETS=({{inputs.targets_bed}})
        FALSE_POSITIVES=({{inputs.false_positives_bed}})
        AMBIGUOUS=({{inputs.ambiguous_beds}}); BAMS=({{inputs.bams}})
        shopt -u nullglob
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        W=$(mktemp -d); cd "$W"
        LINKED=("${{QUERY[0]}}" "${{TRUTH[0]}}" "${{FA[0]}}" "${{FAI[0]}}" "${{REGIONS[@]}}" "${{TARGETS[@]}}"
            "${{FALSE_POSITIVES[@]}}" "${{AMBIGUOUS[@]}}" "${{BAMS[@]}}")
        {LINK_INPUTS}
        OPTIONAL=()
        if [ ${{#REGIONS[@]}} -gt 0 ]; then OPTIONAL+=(-R "$(basename "${{REGIONS[0]}}")"); fi
        if [ ${{#TARGETS[@]}} -gt 0 ]; then OPTIONAL+=(-T "$(basename "${{TARGETS[0]}}")"); fi
        if [ ${{#FALSE_POSITIVES[@]}} -gt 0 ]; then
            OPTIONAL+=(--false-positives "$(basename "${{FALSE_POSITIVES[0]}}")")
        fi
        # As upstream, the files follow a single option each (`--ambiguous a.bed b.bed`, `--bam a.bam b.bam`).
        if [ ${{#AMBIGUOUS[@]}} -gt 0 ]; then
            OPTIONAL+=(--ambiguous); for f in "${{AMBIGUOUS[@]}}"; do OPTIONAL+=("$(basename "$f")"); done
        fi
        if [ ${{#BAMS[@]}} -gt 0 ]; then
            OPTIONAL+=(--bam); for f in "${{BAMS[@]}}"; do OPTIONAL+=("$(basename "$f")"); done
        fi
        som.py \
            "$(basename "${{TRUTH[0]}}")" \
            "$(basename "${{QUERY[0]}}")" \
            $ARGS \
            --reference "$(basename "${{FA[0]}}")" \
            "${{OPTIONAL[@]}}" \
            -o "$P" || exit 1
        shopt -s nullglob
        mv "$P".features.csv "$P".metrics.json "$P".stats.csv {{outputs.results}}/
    """,
)

# As upstream, each ROC table is passed as `<label>:<file>`: `labelled` lists them, in order, by the tables' names
# (staged list inputs arrive in name order).
happy_report_cmd = shell.create(
    name="happy_report",
    image=HAPPY_REPORT_IMAGE,
    resources=REPORT_RESOURCES,
    inputs={"csvs": list[File], "labelled": str, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        CSVS=({inputs.csvs})
        LABELLED={inputs.labelled}
        P={inputs.prefix}
        ARGS={inputs.args}
        W=$(mktemp -d); cd "$W"
        for f in "${CSVS[@]}"; do ln -s "$f" "$(basename "$f")"; done
        rep.py $ARGS $LABELLED -o "$P.html" || exit 1
        mv "$P.html" {outputs.results}/
    """,
)

happy_env = flyte.TaskEnvironment.from_task("happy_happy", happy_happy_cmd.as_task(), happy_sompy_cmd.as_task())
happy_report_env = flyte.TaskEnvironment.from_task("happy_report", happy_report_cmd.as_task())
env = flyte.TaskEnvironment(name="happy", depends_on=[happy_env, happy_report_env])


async def pick(results: Dir, name: str, tool: str) -> File:
    found = await results.get_file(name)
    if found is None:
        raise FileNotFoundError(f"{tool} wrote no {name}")
    return found


@dataclass
class HappyResult:
    summary_csv: File  # <prefix>.summary.csv
    extended_csv: File  # <prefix>.extended.csv
    roc_all_csv: File  # <prefix>.roc.all.csv.gz
    roc_indel_locations_csv: File  # <prefix>.roc.Locations.INDEL.csv.gz
    roc_indel_locations_pass_csv: File  # <prefix>.roc.Locations.INDEL.PASS.csv.gz
    roc_snp_locations_csv: File  # <prefix>.roc.Locations.SNP.csv.gz
    roc_snp_locations_pass_csv: File  # <prefix>.roc.Locations.SNP.PASS.csv.gz
    runinfo: File  # <prefix>.runinfo.json
    metrics_json: File  # <prefix>.metrics.json.gz
    vcf: File | None  # <prefix>.vcf.gz: the truth's and the query's variants, annotated
    tbi: File | None  # <prefix>.vcf.gz.tbi
    results: Dir


async def happy_happy(
    query_vcf: File,
    truth_vcf: File,
    fasta: File,
    fai: File,
    prefix: str,
    regions_bed: File | None = None,
    targets_bed: File | None = None,
    false_positives_bed: File | None = None,
    stratification_tsv: File | None = None,
    stratification_beds: list[File] | None = None,
    args: str = "",
) -> HappyResult:
    """Compare ``query_vcf``'s small variants with ``truth_vcf``'s against ``fasta`` with hap.py.

    ``stratification_tsv`` lists ``<name>\\t<bed>`` strata, naming the
    ``stratification_beds`` by file name; both are needed to stratify.
    """
    results = await happy_happy_cmd(
        query_vcf=query_vcf,
        truth_vcf=truth_vcf,
        regions_bed=regions_bed,
        targets_bed=targets_bed,
        fasta=fasta,
        fai=fai,
        false_positives_bed=false_positives_bed,
        stratification_tsv=stratification_tsv,
        stratification_beds=stratification_beds or [],
        prefix=prefix,
        args=args,
    )
    names = {
        "summary_csv": "summary.csv",
        "extended_csv": "extended.csv",
        "roc_all_csv": "roc.all.csv.gz",
        "roc_indel_locations_csv": "roc.Locations.INDEL.csv.gz",
        "roc_indel_locations_pass_csv": "roc.Locations.INDEL.PASS.csv.gz",
        "roc_snp_locations_csv": "roc.Locations.SNP.csv.gz",
        "roc_snp_locations_pass_csv": "roc.Locations.SNP.PASS.csv.gz",
        "runinfo": "runinfo.json",
        "metrics_json": "metrics.json.gz",
    }
    files = {key: await pick(results, f"{prefix}.{name}", "hap.py") for key, name in names.items()}
    return HappyResult(
        **files,
        vcf=await results.get_file(f"{prefix}.vcf.gz"),
        tbi=await results.get_file(f"{prefix}.vcf.gz.tbi"),
        results=results,
    )


@dataclass
class SompyResult:
    stats: File  # <prefix>.stats.csv
    metrics: File  # <prefix>.metrics.json
    features: File | None  # <prefix>.features.csv, when args asks for a feature table


async def happy_sompy(
    query_vcf: File,
    truth_vcf: File,
    fasta: File,
    fai: File,
    prefix: str,
    regions_bed: File | None = None,
    targets_bed: File | None = None,
    false_positives_bed: File | None = None,
    ambiguous_beds: list[File] | None = None,
    bams: list[File] | None = None,
    args: str = "",
) -> SompyResult:
    """Compare ``query_vcf``'s somatic calls with ``truth_vcf``'s against ``fasta`` with som.py."""
    results = await happy_sompy_cmd(
        query_vcf=query_vcf,
        truth_vcf=truth_vcf,
        regions_bed=regions_bed,
        targets_bed=targets_bed,
        fasta=fasta,
        fai=fai,
        false_positives_bed=false_positives_bed,
        ambiguous_beds=ambiguous_beds or [],
        bams=bams or [],
        prefix=prefix,
        args=args,
    )
    return SompyResult(
        stats=await pick(results, f"{prefix}.stats.csv", "som.py"),
        metrics=await pick(results, f"{prefix}.metrics.json", "som.py"),
        features=await results.get_file(f"{prefix}.features.csv"),
    )


async def happy_report(csvs: list[File], labels: list[str], prefix: str, args: str = "") -> File:
    """An HTML report of hap.py's ROC tables ``csvs``: ``<prefix>.html``.

    Each table is labelled by its method, as upstream builds it from a
    method and an optional comparison method (``gatk4_vcfeval``); the tables
    must have distinct names.
    """
    if len(csvs) != len(labels):
        raise ValueError(f"{len(csvs)} ROC tables but {len(labels)} labels")
    names = [csv.name for csv in csvs]
    if len(set(names)) != len(names):
        raise ValueError(f"the ROC tables' names must be distinct: {names}")
    if any(" " in label or not label for label in labels):
        raise ValueError(f"labels must be non-empty and without spaces: {labels}")
    labelled = " ".join(f"{label}:{name}" for label, name in zip(labels, names, strict=True))
    results = await happy_report_cmd(csvs=csvs, labelled=labelled, prefix=prefix, args=args)
    return await pick(results, f"{prefix}.html", "rep.py")
