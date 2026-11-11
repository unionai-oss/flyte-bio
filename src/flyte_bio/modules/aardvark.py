"""aardvark — Aardvark, small variant benchmarking by comparing haplotype sequences.

Exposes :func:`aardvark_compare`, which compares a query VCF's variants with a
truth VCF's against a reference within a BED's regions, optionally
stratified by the regions a TSV names (``--stratification``; the TSV's BED
paths are relative to the working directory, where the ``stratification_beds``
are linked). As upstream, the outputs (with the debug outputs, written to the
same directory) are renamed ``<prefix>.<name>``:

- ``<prefix>.summary.tsv``: precision and recall by variant type (and
  stratification);
- ``<prefix>.truth.vcf.gz`` / ``<prefix>.query.vcf.gz``: the variants,
  labelled with their match status;
- ``<prefix>.region_sequences.tsv.gz`` / ``<prefix>.region_summary.tsv.gz``:
  the compared haplotype sequences, and the matches per region;
- ``<prefix>.cli_settings.json``: the run's settings.

The inputs are linked into the working directory and named relative to it,
since the settings record them; they record the thread count (the task's
CPUs) too. :func:`aardvark_compare_result` reads the outputs of a task run with
other resources (``aardvark_compare_cmd.override(resources=...)``).
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

AARDVARK_IMAGE = "community.wave.seqera.io/library/aardvark:0.10.5--50402c072fe171c0"
# Upstream's process_low.
DEFAULT_RESOURCES = flyte.Resources(cpu=2, memory="8Gi")

aardvark_compare_cmd = shell.create(
    name="aardvark_compare",
    image=AARDVARK_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "query_vcf": File,
        "query_tbi": File,
        "truth_vcf": File,
        "truth_tbi": File,
        "regions_bed": File,
        "fasta": File,
        "stratification_tsv": File | None,
        "stratification_beds": list[File],
        "prefix": str,
        "args": str,
    },
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        QUERY=({{inputs.query_vcf}}); QUERY_TBI=({{inputs.query_tbi}})
        TRUTH=({{inputs.truth_vcf}}); TRUTH_TBI=({{inputs.truth_tbi}})
        REGIONS=({{inputs.regions_bed}}); FA=({{inputs.fasta}})
        shopt -s nullglob; STRAT_TSV=({{inputs.stratification_tsv}}); STRAT_BEDS=({{inputs.stratification_beds}})
        shopt -u nullglob
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        for f in "${{QUERY[0]}}" "${{QUERY_TBI[0]}}" "${{TRUTH[0]}}" "${{TRUTH_TBI[0]}}" "${{REGIONS[0]}}" \
            "${{FA[0]}}" "${{STRAT_TSV[@]}}" "${{STRAT_BEDS[@]}}"; do
            ln -s "$f" "$(basename "$f")"
        done
        STRATIFICATION=()
        if [ ${{#STRAT_TSV[@]}} -gt 0 ] && [ ${{#STRAT_BEDS[@]}} -gt 0 ]; then
            STRATIFICATION=(--stratification "$(basename "${{STRAT_TSV[0]}}")")
        fi
        # aardvark is a Rust tool: a panic's exit code could be retried forever, so any failure exits 1.
        aardvark compare \
            --threads $TASK_CPUS \
            --reference "$(basename "${{FA[0]}}")" \
            --truth-vcf "$(basename "${{TRUTH[0]}}")" \
            --query-vcf "$(basename "${{QUERY[0]}}")" \
            --regions "$(basename "${{REGIONS[0]}}")" \
            "${{STRATIFICATION[@]}}" \
            --output-dir "$P" \
            --output-debug "$P" \
            $ARGS || exit 1
        for f in "$P"/*; do mv "$f" "$OUT/$P.${{f##*/}}"; done
    """,
)

env = flyte.TaskEnvironment.from_task("aardvark", aardvark_compare_cmd.as_task())


@dataclass
class AardvarkCompareResult:
    summary: File  # <prefix>.summary.tsv
    labelled_truth: File  # <prefix>.truth.vcf.gz
    labelled_query: File  # <prefix>.query.vcf.gz
    runinfo: File  # <prefix>.cli_settings.json
    region_sequences: File  # <prefix>.region_sequences.tsv.gz
    region_summary: File  # <prefix>.region_summary.tsv.gz
    results: Dir  # everything aardvark wrote


async def aardvark_compare_result(results: Dir, prefix: str) -> AardvarkCompareResult:
    """The outputs ``aardvark_compare_cmd`` wrote to ``results``."""

    async def pick(name: str) -> File:
        found = await results.get_file(f"{prefix}.{name}")
        if found is None:
            raise FileNotFoundError(f"aardvark compare wrote no {prefix}.{name}")
        return found

    return AardvarkCompareResult(
        summary=await pick("summary.tsv"),
        labelled_truth=await pick("truth.vcf.gz"),
        labelled_query=await pick("query.vcf.gz"),
        runinfo=await pick("cli_settings.json"),
        region_sequences=await pick("region_sequences.tsv.gz"),
        region_summary=await pick("region_summary.tsv.gz"),
        results=results,
    )


async def aardvark_compare(
    query_vcf: File,
    query_tbi: File,
    truth_vcf: File,
    truth_tbi: File,
    regions_bed: File,
    fasta: File,
    prefix: str,
    stratification_tsv: File | None = None,
    stratification_beds: list[File] | None = None,
    args: str = "",
) -> AardvarkCompareResult:
    """Compare ``query_vcf``'s variants with ``truth_vcf``'s against ``fasta`` within ``regions_bed``.

    ``stratification_tsv`` lists ``<name>\\t<bed>`` strata, naming the
    ``stratification_beds`` by file name; both are needed to stratify.
    """
    results = await aardvark_compare_cmd(
        query_vcf=query_vcf,
        query_tbi=query_tbi,
        truth_vcf=truth_vcf,
        truth_tbi=truth_tbi,
        regions_bed=regions_bed,
        fasta=fasta,
        stratification_tsv=stratification_tsv,
        stratification_beds=stratification_beds or [],
        prefix=prefix,
        args=args,
    )
    return await aardvark_compare_result(results, prefix)
