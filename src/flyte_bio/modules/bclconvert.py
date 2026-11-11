"""bclconvert — Illumina BCL Convert: BCL-to-FASTQ conversion and demultiplexing.

Exposes :func:`bclconvert`, which demultiplexes an Illumina run folder,
given as a ``.tar.gz`` archive, with its (v1 or v2) sample sheet. As
upstream, the archive is extracted into ``<archive name without .tar.gz>/``,
dropping its top-level directory when it has just one, and ``args2`` /
``args3`` are passed to ``tar`` before and after the archive name.

BCL Convert writes one FASTQ per sample, read and lane, which can be many, so
they come back as one archive. The output Dir holds:

- ``fastq.tar.gz``: every FASTQ BCL Convert wrote, under the paths it gave
  them (``<sample>_S<n>_L00<l>_R1_001.fastq.gz``, ``Undetermined_S0_...``,
  index reads ``..._I1_001.fastq.gz``, in ``<project>/`` folders when the
  sample sheet names projects); samples with no reads get empty files;
- ``Reports/``: ``Demultiplex_Stats.csv``, ``Quality_Metrics.csv``,
  ``Top_Unknown_Barcodes.csv``, ``fastq_list.csv``, ``RunInfo.xml``,
  ``SampleSheet.csv``, ...;
- ``Logs/``: ``Info.log``, ``Warnings.log``, ``Errors.log``,
  ``FastqComplete.txt``;
- ``InterOp/``: the run folder's InterOp ``*.bin`` metrics, as upstream
  copies them.

BCL Convert's conversion and compression threads are the task's CPUs.

BCL Convert is Illumina's proprietary software, distributed under Illumina's
end-user licence; the container is upstream's build of it.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

# Upstream's own build of Illumina's BCL Convert.
BCLCONVERT_IMAGE = "quay.io/nf-core/bclconvert:4.6.2"
# Upstream's process_high asks for 12 CPUs and 72 GB; raise for big runs with
# bclconvert_cmd.override(resources=...).
DEFAULT_RESOURCES = flyte.Resources(cpu=8, memory="16Gi", disk="50Gi")

bclconvert_cmd = shell.create(
    name="bclconvert",
    image=BCLCONVERT_IMAGE,
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
        bcl-convert $ARGS \
            --output-directory output \
            --bcl-input-directory "$NAME" \
            --sample-sheet "${{SHEET[0]}}" \
            --bcl-num-conversion-threads $TASK_CPUS \
            --bcl-num-compression-threads $TASK_CPUS
        mkdir -p output/InterOp
        cp -n "$NAME"/InterOp/*.bin output/InterOp/
        (cd output && find . -name '*.fastq.gz' | sort | tar -czf "$OUT/fastq.tar.gz" -T -)
        cp -r output/Reports output/Logs output/InterOp "$OUT"/
    """,
)

env = flyte.TaskEnvironment.from_task("bclconvert", bclconvert_cmd.as_task())


@dataclass
class BclconvertResult:
    fastq: File  # fastq.tar.gz: every FASTQ written
    demultiplex_stats: File  # Reports/Demultiplex_Stats.csv
    results: Dir  # fastq.tar.gz, Reports/, Logs/ and InterOp/


async def bclconvert(
    samplesheet: File, run_dir: File, args: str = "", args2: str = "", args3: str = ""
) -> BclconvertResult:
    """Demultiplex the run folder archive ``run_dir`` (``.tar.gz``) with ``samplesheet``."""
    results = await bclconvert_cmd(samplesheet=samplesheet, run_dir=run_dir, args=args, args2=args2, args3=args3)

    async def pick(name: str) -> File:
        found = await results.get_file(name)
        if found is None:
            raise FileNotFoundError(f"bcl-convert wrote no {name}")
        return found

    return BclconvertResult(
        fastq=await pick("fastq.tar.gz"),
        demultiplex_stats=await pick("Reports/Demultiplex_Stats.csv"),
        results=results,
    )
