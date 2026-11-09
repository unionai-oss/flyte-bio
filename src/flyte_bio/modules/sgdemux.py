"""sgdemux — Singular Genomics G4 FASTQ demultiplexing.

Exposes :func:`sgdemux`, which demultiplexes a folder of (bgzipped) FASTQs,
given as a ``.tar.gz`` archive (extracted into ``input/``, dropping its
top-level directory when it has just one, as upstream's untar step does),
with a sample sheet (``Sample_Barcode`` and ``Sample_ID`` columns, plus
optional ``Lane``, ``Project``, ...).

sgdemux writes one FASTQ per sample and read, which can be many, so they come
back as one archive. The output Dir holds:

- ``fastq.tar.gz``: ``<sample>_S<n>_L00<l>_R1_001.fastq.gz``,
  ``..._R2_001.fastq.gz`` and ``Undetermined_S<n>_...``; the order of the
  reads within them depends on the thread count;
- ``metrics.tsv``, ``per_sample_metrics.tsv``, ``per_project_metrics.tsv``,
  ``most_frequent_unmatched.tsv`` and ``sample_barcode_hop_metrics.tsv``.

As upstream, the demultiplexing, compressor and writer threads are each the
task's CPUs.

sgdemux is Singular Genomics' software; its licence allows use with the
Singular G4 sequencing platform only.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

# Upstream's main.nf names the quay.io biocontainer directly.
SGDEMUX_IMAGE = "quay.io/biocontainers/sgdemux:1.1.1--ha982bd6_0"
# Upstream's process_high asks for 12 CPUs and 72 GB.
DEFAULT_RESOURCES = flyte.Resources(cpu=8, memory="16Gi")

sgdemux_cmd = shell.create(
    name="sgdemux",
    image=SGDEMUX_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"sample_sheet": File, "fastqs": File, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        set -e
        SHEET=({{inputs.sample_sheet}})
        ARCHIVE=({{inputs.fastqs}})
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        mkdir input "$P"
        # As upstream: drop the top-level directory only when the archive has just one.
        if [[ $(tar -tzf "${{ARCHIVE[0]}}" | grep / | cut -d/ -f1 | sort -u | wc -l) -eq 1 ]]; then
            tar -C input --strip-components 1 -xzf "${{ARCHIVE[0]}}"
        else
            tar -C input -xzf "${{ARCHIVE[0]}}"
        fi
        # sgdemux is a Rust tool: a panic exits 255, which the platform would retry forever.
        sgdemux \
            --sample-metadata "${{SHEET[0]}}" \
            --fastqs input \
            --output-dir "$P" \
            --demux-threads $TASK_CPUS \
            --compressor-threads $TASK_CPUS \
            --writer-threads $TASK_CPUS \
            $ARGS || exit 1
        (cd "$P" && ls | grep '\.fastq\.gz$' | tar -czf "$OUT/fastq.tar.gz" -T -)
        cp "$P"/*.tsv "$OUT"/
    """,
)

env = flyte.TaskEnvironment.from_task("sgdemux", sgdemux_cmd.as_task())


@dataclass
class SgdemuxResult:
    fastq: File  # fastq.tar.gz: per-sample and undetermined FASTQs
    metrics: File  # metrics.tsv
    per_sample_metrics: File  # per_sample_metrics.tsv
    per_project_metrics: File  # per_project_metrics.tsv
    most_frequent_unmatched: File  # most_frequent_unmatched.tsv
    sample_barcode_hop_metrics: File  # sample_barcode_hop_metrics.tsv
    results: Dir  # all of the above


async def sgdemux(sample_sheet: File, fastqs: File, prefix: str = "sample", args: str = "") -> SgdemuxResult:
    """Demultiplex the FASTQs in archive ``fastqs`` by the samples in ``sample_sheet``."""
    results = await sgdemux_cmd(sample_sheet=sample_sheet, fastqs=fastqs, prefix=prefix, args=args)

    async def pick(name: str) -> File:
        found = await results.get_file(name)
        if found is None:
            raise FileNotFoundError(f"sgdemux wrote no {name}")
        return found

    return SgdemuxResult(
        fastq=await pick("fastq.tar.gz"),
        metrics=await pick("metrics.tsv"),
        per_sample_metrics=await pick("per_sample_metrics.tsv"),
        per_project_metrics=await pick("per_project_metrics.tsv"),
        most_frequent_unmatched=await pick("most_frequent_unmatched.tsv"),
        sample_barcode_hop_metrics=await pick("sample_barcode_hop_metrics.tsv"),
        results=results,
    )
