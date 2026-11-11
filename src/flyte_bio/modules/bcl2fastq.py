"""bcl2fastq — Illumina BCL-to-FASTQ conversion and demultiplexing.

Exposes :func:`bcl2fastq`, which demultiplexes an Illumina run folder,
given as a ``.tar.gz`` archive, with its sample sheet. As upstream, the
archive is extracted into ``<archive name without .tar.gz>/``, dropping its
top-level directory when it has just one, and ``args2`` / ``args3`` are
passed to ``tar`` before and after the archive name.

bcl2fastq writes one FASTQ per sample, read and lane, which can be many, so
they come back as one archive. The output Dir holds:

- ``fastq.tar.gz``: every FASTQ bcl2fastq wrote, under the paths it gave them
  (``<sample>_S<n>_L00<l>_R1_001.fastq.gz``, ``Undetermined_S0_...``, index
  reads ``..._I1_001.fastq.gz`` with ``--create-fastq-for-index-reads``, in
  ``<project>/`` folders when the sample sheet names projects);
- ``Reports.tar.gz``: the ``Reports/`` folder (an HTML page per project,
  sample and barcode);
- ``Stats/``: ``Stats.json``, ``ConversionStats.xml``,
  ``DemultiplexingStats.xml``, ``AdapterTrimming.txt`` and per-lane
  ``DemuxSummary`` / ``FastqSummary`` tables;
- ``InterOp/``: the run folder's InterOp ``*.bin`` metrics, as upstream copies
  them;
- ``RunInfo.xml``: the run folder's.

``--processing-threads`` is the task's CPUs. bcl2fastq's statistics (its
unknown-barcode counts, and so ``Stats.json``, ``ConversionStats.xml``,
``DemuxSummary`` and the reports) change with 4 or more threads on small
runs; the FASTQs don't.

bcl2fastq is Illumina's proprietary software, distributed under Illumina's
end-user licence; the container is upstream's build of it.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

# Upstream's own build of Illumina's bcl2fastq.
BCL2FASTQ_IMAGE = "quay.io/nf-core/bcl2fastq:2.20.0.422"
# Upstream's process_high asks for 12 CPUs and 72 GB; raise for big runs with
# bcl2fastq_cmd.override(resources=...).
DEFAULT_RESOURCES = flyte.Resources(cpu=8, memory="16Gi", disk="50Gi")

bcl2fastq_cmd = shell.create(
    name="bcl2fastq",
    image=BCL2FASTQ_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"samplesheet": File, "run_dir": File, "args": str, "args2": str, "args3": str},
    defaults={"args": "", "args2": "", "args3": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        set -e
        SHEET=({{inputs.samplesheet}})
        RUN=({{inputs.run_dir}})
        ARGS={{inputs.args}}
        ARGS2={{inputs.args2}}
        ARGS3={{inputs.args3}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        NAME=$(basename "${{RUN[0]}}" .tar.gz)
        mkdir -p "$NAME"
        # As upstream: drop the top-level directory only when the archive has just one.
        if [[ $(tar -taf "${{RUN[0]}}" | grep -o -P "^.*?\/" | uniq | wc -l) -eq 1 ]]; then
            tar -C "$NAME" --strip-components 1 -xaf $ARGS2 "${{RUN[0]}}" $ARGS3
        else
            tar -C "$NAME" -xaf $ARGS2 "${{RUN[0]}}" $ARGS3
        fi
        bcl2fastq $ARGS \
            --output-dir output \
            --runfolder-dir "$NAME" \
            --sample-sheet "${{SHEET[0]}}" \
            --processing-threads $TASK_CPUS
        (cd output && find . -name '*.fastq.gz' | sort | tar -czf "$OUT/fastq.tar.gz" -T -)
        tar -C output -czf "$OUT/Reports.tar.gz" Reports
        cp -r output/Stats "$OUT"/
        mkdir -p "$OUT/InterOp"
        cp "$NAME"/InterOp/*.bin "$OUT/InterOp"/
        cp "$NAME/RunInfo.xml" "$OUT"/
    """,
)

env = flyte.TaskEnvironment.from_task("bcl2fastq", bcl2fastq_cmd.as_task())


@dataclass
class Bcl2fastqResult:
    fastq: File  # fastq.tar.gz: every FASTQ written
    reports: File  # Reports.tar.gz
    stats: File  # Stats/Stats.json
    results: Dir  # everything above, plus Stats/, InterOp/ and RunInfo.xml


async def bcl2fastq(
    samplesheet: File, run_dir: File, args: str = "", args2: str = "", args3: str = ""
) -> Bcl2fastqResult:
    """Demultiplex the run folder archive ``run_dir`` (``.tar.gz``) with ``samplesheet``."""
    results = await bcl2fastq_cmd(samplesheet=samplesheet, run_dir=run_dir, args=args, args2=args2, args3=args3)

    async def pick(name: str) -> File:
        found = await results.get_file(name)
        if found is None:
            raise FileNotFoundError(f"bcl2fastq wrote no {name}")
        return found

    return Bcl2fastqResult(
        fastq=await pick("fastq.tar.gz"),
        reports=await pick("Reports.tar.gz"),
        stats=await pick("Stats/Stats.json"),
        results=results,
    )
