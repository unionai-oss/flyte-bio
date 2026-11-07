"""fastqdl — download FASTQ files from ENA or SRA with fastq-dl.

Exposes :func:`fastqdl`, which downloads the reads of one accession: a run
(``ERR...``, ``SRR...``), an experiment (``ERX...``), a sample, a study or a
BioProject. Outputs keep fastq-dl's names:

- ``<run>.fastq.gz`` (single-end or unpaired reads) and/or
  ``<run>_1.fastq.gz`` / ``<run>_2.fastq.gz``, one set per run (fastq-dl
  names reads by run; with ``--group-by-experiment`` or
  ``--group-by-sample`` in ``args`` it merges them by experiment or sample);
- ``<prefix>-run-info.tsv``: the runs' metadata;
- ``<prefix>-run-mergers.tsv``: which runs were merged, when grouping.

fastq-dl downloads from ENA first and falls back to SRA (``--provider`` in
``args`` picks one), retrying failed downloads itself. It uses the task's CPUs.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

FASTQDL_IMAGE = "community.wave.seqera.io/library/fastq-dl:4.0.1--0758300a069ba945"
# Upstream labels fastq-dl process_medium. Raise disk for large studies.
DEFAULT_RESOURCES = flyte.Resources(cpu=4, memory="8Gi", disk="50Gi")

fastqdl_cmd = shell.create(
    name="fastqdl",
    image=FASTQDL_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"accession": str, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        ACCESSION={{inputs.accession}}
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"

        fastq-dl \
            $ARGS \
            --prefix "$P" \
            --accession "$ACCESSION" \
            --cpus $TASK_CPUS \
            --outdir . \
            || exit 1
        shopt -s nullglob
        mv ./*.fastq.gz ./*-run-info.tsv ./*-run-mergers.tsv "$OUT"/
    """,
)

env = flyte.TaskEnvironment.from_task("fastqdl", fastqdl_cmd.as_task())


@dataclass
class FastqdlResult:
    fastq: list[File]  # *.fastq.gz, sorted by name
    runinfo: File  # <prefix>-run-info.tsv
    runmergers: File | None  # <prefix>-run-mergers.tsv, when runs were merged
    results: Dir  # everything fastq-dl wrote


async def fastqdl(accession: str, prefix: str = "sample", args: str = "") -> FastqdlResult:
    """Download the reads of ``accession`` (a run, experiment, sample, study or BioProject)."""
    results = await fastqdl_cmd(accession=accession, prefix=prefix, args=args)
    files = {f.name: f async for f in results.walk()}
    runinfo = files.get(f"{prefix}-run-info.tsv")
    if runinfo is None:
        raise FileNotFoundError(f"fastq-dl wrote no {prefix}-run-info.tsv")
    fastq = [files[n] for n in sorted(files) if n.endswith(".fastq.gz")]
    if not fastq:
        raise FileNotFoundError("fastq-dl wrote no .fastq.gz")
    return FastqdlResult(
        fastq=fastq, runinfo=runinfo, runmergers=files.get(f"{prefix}-run-mergers.tsv"), results=results
    )
