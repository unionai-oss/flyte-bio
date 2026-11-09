"""bases2fastq — Element Biosciences AVITI bases-to-FASTQ conversion and demultiplexing.

Exposes :func:`bases2fastq`, which demultiplexes an Element run folder, given
as a ``.tar.gz`` archive, with an optional run manifest (``RunManifest.csv``;
without one, bases2fastq uses the run folder's ``RunManifest.json``). As
upstream's test, the archive is extracted into ``input_dir/``, dropping its
top-level directory when it has just one.

bases2fastq writes its output folder ``<prefix>/``. The output Dir holds:

- its top-level files: ``Metrics.csv``, ``RunManifest.json`` (the manifest
  as applied), ``UnassignedSequences.csv``, ``RunStats.json``,
  ``<run name>_QC.html``, ``multiqc_report.html``, ``IndexAssignment.csv``,
  ...;
- ``Samples.tar.gz``: the ``Samples/`` folder, one archive since it holds a
  folder per project and sample: ``Samples/<project>/<sample>/<sample>_R1.fastq.gz``,
  ``..._R2.fastq.gz`` and ``<sample>_stats.json``, plus per-project QC
  reports.

``-p`` (threads) is the task's CPUs.

bases2fastq is Element Biosciences' proprietary software, distributed under
Element's end-user licence; the container is Element's own image.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

# Element's own image, as upstream.
BASES2FASTQ_IMAGE = "docker.io/elembio/bases2fastq:2.3.0"
# Upstream's process_high asks for 12 CPUs and 72 GB; raise for big runs with
# bases2fastq_cmd.as_task().override(resources=...).
DEFAULT_RESOURCES = flyte.Resources(cpu=8, memory="16Gi", disk="50Gi")

bases2fastq_cmd = shell.create(
    name="bases2fastq",
    image=BASES2FASTQ_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"run_manifest": File | None, "run_dir": File, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        set -e
        shopt -s nullglob; MANIFEST=({{inputs.run_manifest}}); shopt -u nullglob
        RUN=({{inputs.run_dir}})
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        mkdir input_dir
        # As upstream: drop the top-level directory only when the archive has just one.
        if [[ $(tar -tzf "${{RUN[0]}}" | grep / | cut -d/ -f1 | sort -u | wc -l) -eq 1 ]]; then
            tar -C input_dir --strip-components 1 -xzf "${{RUN[0]}}"
        else
            tar -C input_dir -xzf "${{RUN[0]}}"
        fi
        EXTRA=()
        if [ ${{#MANIFEST[@]}} -gt 0 ]; then EXTRA+=(-r "${{MANIFEST[0]}}"); fi
        bases2fastq -p $TASK_CPUS "${{EXTRA[@]}}" $ARGS input_dir "$P"
        find "$P" -maxdepth 1 -type f -exec cp {{}} "$OUT"/ \;
        tar -C "$P" -czf "$OUT/Samples.tar.gz" Samples
    """,
)

env = flyte.TaskEnvironment.from_task("bases2fastq", bases2fastq_cmd.as_task())


@dataclass
class Bases2fastqResult:
    samples: File  # Samples.tar.gz: per-sample FASTQs and stats
    metrics: File  # Metrics.csv
    run_stats: File  # RunStats.json
    run_manifest: File  # RunManifest.json, the manifest as applied
    unassigned: File  # UnassignedSequences.csv
    results: Dir  # the output folder's top-level files and Samples.tar.gz


async def bases2fastq(
    run_dir: File, run_manifest: File | None = None, prefix: str = "sample", args: str = ""
) -> Bases2fastqResult:
    """Demultiplex the run folder archive ``run_dir`` (``.tar.gz``), with ``run_manifest`` if given."""
    results = await bases2fastq_cmd(run_manifest=run_manifest, run_dir=run_dir, prefix=prefix, args=args)

    async def pick(name: str) -> File:
        found = await results.get_file(name)
        if found is None:
            raise FileNotFoundError(f"bases2fastq wrote no {name}")
        return found

    return Bases2fastqResult(
        samples=await pick("Samples.tar.gz"),
        metrics=await pick("Metrics.csv"),
        run_stats=await pick("RunStats.json"),
        run_manifest=await pick("RunManifest.json"),
        unassigned=await pick("UnassignedSequences.csv"),
        results=results,
    )
