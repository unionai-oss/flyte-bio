"""svanalyzer — SVanalyzer, structural variant comparison.

Exposes :func:`svanalyzer_svbenchmark`, which compares a VCF's SVs with a
truth VCF's against a reference (``svanalyzer benchmark``), optionally within
a BED's regions (``-includebed``), writing:

- ``<prefix>.falsenegatives.vcf.gz`` / ``<prefix>.falsepositives.vcf.gz``:
  the truth's SVs the calls miss, and the calls not in the truth (bgzipped,
  as upstream, with ``bgzip_args`` and the task's CPUs as threads);
- ``<prefix>.distances``: the distances between matched SVs;
- ``<prefix>.report``: recall, precision and F1;
- ``<prefix>.log``.

The VCFs, their indexes and the reference with its ``.fai`` are linked side by
side, and named relative to the working directory, as upstream.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

# Upstream's main.nf names the bare biocontainer (Docker Hub); this is its quay.io mirror.
SVANALYZER_IMAGE = "quay.io/biocontainers/svanalyzer:0.36--pl526_0"
# Upstream's process_single.
DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory="6Gi")

svanalyzer_svbenchmark_cmd = shell.create(
    name="svanalyzer_svbenchmark",
    image=SVANALYZER_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "vcf": File,
        "tbi": File,
        "truth_vcf": File,
        "truth_tbi": File,
        "bed": File | None,
        "fasta": File,
        "fai": File,
        "prefix": str,
        "args": str,
        "bgzip_args": str,
    },
    defaults={"args": "", "bgzip_args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        VCF=({{inputs.vcf}}); TBI=({{inputs.tbi}}); TRUTH=({{inputs.truth_vcf}}); TRUTH_TBI=({{inputs.truth_tbi}})
        FA=({{inputs.fasta}}); FAI=({{inputs.fai}})
        shopt -s nullglob; BED=({{inputs.bed}}); shopt -u nullglob
        P={{inputs.prefix}}
        ARGS={{inputs.args}}; BGZIP_ARGS={{inputs.bgzip_args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        for f in "${{VCF[0]}}" "${{TBI[0]}}" "${{TRUTH[0]}}" "${{TRUTH_TBI[0]}}" "${{FA[0]}}" "${{FAI[0]}}" \
            "${{BED[@]}}"; do
            ln -s "$f" "$(basename "$f")"
        done
        REGIONS=()
        if [ ${{#BED[@]}} -gt 0 ]; then REGIONS=(-includebed "$(basename "${{BED[0]}}")"); fi
        # svanalyzer is Perl: a die exits 255, which would be retried forever.
        svanalyzer benchmark \
            $ARGS \
            --ref "$(basename "${{FA[0]}}")" \
            --test "$(basename "${{VCF[0]}}")" \
            --truth "$(basename "${{TRUTH[0]}}")" \
            --prefix "$P" \
            "${{REGIONS[@]}}" || exit 1
        for kind in falsenegatives falsepositives; do
            bgzip $BGZIP_ARGS --threads $TASK_CPUS -c "$P.$kind.vcf" > "$OUT/$P.$kind.vcf.gz" || exit 1
        done
        mv "$P.distances" "$P.log" "$P.report" "$OUT"/
    """,
)

env = flyte.TaskEnvironment.from_task("svanalyzer", svanalyzer_svbenchmark_cmd.as_task())


@dataclass
class SvbenchmarkResult:
    fns: File  # <prefix>.falsenegatives.vcf.gz
    fps: File  # <prefix>.falsepositives.vcf.gz
    distances: File  # <prefix>.distances
    log: File  # <prefix>.log
    report: File  # <prefix>.report


async def svanalyzer_svbenchmark(
    vcf: File,
    tbi: File,
    truth_vcf: File,
    truth_tbi: File,
    fasta: File,
    fai: File,
    prefix: str,
    bed: File | None = None,
    args: str = "",
    bgzip_args: str = "",
) -> SvbenchmarkResult:
    """Compare ``vcf``'s SVs with ``truth_vcf``'s against ``fasta``, within ``bed`` when given."""
    results = await svanalyzer_svbenchmark_cmd(
        vcf=vcf,
        tbi=tbi,
        truth_vcf=truth_vcf,
        truth_tbi=truth_tbi,
        bed=bed,
        fasta=fasta,
        fai=fai,
        prefix=prefix,
        args=args,
        bgzip_args=bgzip_args,
    )

    async def pick(name: str) -> File:
        found = await results.get_file(f"{prefix}.{name}")
        if found is None:
            raise FileNotFoundError(f"svanalyzer benchmark wrote no {prefix}.{name}")
        return found

    return SvbenchmarkResult(
        fns=await pick("falsenegatives.vcf.gz"),
        fps=await pick("falsepositives.vcf.gz"),
        distances=await pick("distances"),
        log=await pick("log"),
        report=await pick("report"),
    )
