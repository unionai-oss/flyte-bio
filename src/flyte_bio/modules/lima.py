"""lima — PacBio's barcode and primer demultiplexer: primer removal from CCS reads.

Exposes :func:`lima`, which runs ``lima`` on one sample's CCS reads (BAM,
FASTA or FASTQ, optionally gzipped) against a primer FASTA. lima writes, as
upstream, with ``<prefix>`` as the output stem:

- reports: ``<prefix>.lima.counts``, ``<prefix>.lima.report``,
  ``<prefix>.lima.summary`` and, depending on ``args``, ``<prefix>.lima.clips``
  and ``<prefix>.lima.guess`` (``--peek-guess``);
- reads, one set per primer pair: ``<prefix>.<5p>--<3p>.bam`` with its
  ``.bam.pbi`` index, and for FASTA/FASTQ input
  ``<prefix>.<5p>--<3p>.{fasta,fastq}[.gz]`` in the input's format;
- dataset files: ``*.xml`` and ``*.json``.

As upstream, the output's format follows the input's extension
(``.bam``, ``.fasta[.gz]``, ``.fastq[.gz]``), lima gets the task's CPUs as
``-j``, and the inputs are linked into a working directory under their own
names, so the BAM header's ``@PG`` command line names them as upstream's
does. ``args`` go to lima (the upstream pipeline uses ``--isoseq
--peek-guess``).
"""

import os
from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

LIMA_IMAGE = "quay.io/biocontainers/lima:2.12.0--h9ee0642_1"
# Upstream labels it process_low (2 CPUs). The output BAM's @PG command line records -j.
DEFAULT_RESOURCES = flyte.Resources(cpu=2, memory="8Gi")

READ_EXTENSIONS = ("bam", "fasta", "fasta.gz", "fastq", "fastq.gz")

lima_cmd = shell.create(
    name="lima",
    image=LIMA_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"ccs": File, "primers": File, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        CCS=({{inputs.ccs}}); PRIMERS=({{inputs.primers}})
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        NAME=$(basename "${{CCS[0]}}"); ln -s "${{CCS[0]}}" "$NAME"
        PRIMERS_NAME=$(basename "${{PRIMERS[0]}}"); ln -s "${{PRIMERS[0]}}" "$PRIMERS_NAME"
        OUT_EXT=""
        if [[ $NAME =~ bam$ ]]; then
            OUT_EXT="bam"
        elif [[ $NAME =~ fasta$ ]]; then
            OUT_EXT="fasta"
        elif [[ $NAME =~ fasta.gz$ ]]; then
            OUT_EXT="fasta.gz"
        elif [[ $NAME =~ fastq$ ]]; then
            OUT_EXT="fastq"
        elif [[ $NAME =~ fastq.gz$ ]]; then
            OUT_EXT="fastq.gz"
        fi
        lima "$NAME" "$PRIMERS_NAME" "$P.$OUT_EXT" -j $TASK_CPUS $ARGS || exit 1
        # Everything lima wrote; the inputs are symlinks.
        for f in *; do
            if [ -f "$f" ] && [ ! -L "$f" ]; then mv "$f" "$OUT"/; fi
        done
    """,
)

env = flyte.TaskEnvironment.from_task("lima", lima_cmd.as_task())


@dataclass
class LimaResult:
    counts: File  # <prefix>.lima.counts
    report: File  # <prefix>.lima.report
    summary: File  # <prefix>.lima.summary
    clips: File | None  # <prefix>.lima.clips
    guess: File | None  # <prefix>.lima.guess, with --peek-guess
    bams: list[File]  # <prefix>.<pair>.bam, one per primer pair
    pbis: list[File]  # their .bam.pbi indexes
    reads: list[File]  # <prefix>.<pair>.{fasta,fastq}[.gz], for FASTA/FASTQ input
    results: Dir  # everything lima wrote


async def lima(ccs: File, primers: File, prefix: str = "sample", args: str = "") -> LimaResult:
    """Remove ``primers`` from CCS reads ``ccs`` (BAM, FASTA or FASTQ, optionally gzipped)."""
    name = ccs.name or os.path.basename(ccs.path)
    if not name.endswith(READ_EXTENSIONS):
        raise ValueError(f"{name}: lima reads .bam, .fasta[.gz] or .fastq[.gz]")
    if name in {f"{prefix}.{ext}" for ext in READ_EXTENSIONS}:
        raise ValueError("Input and output names are the same, set another prefix")
    results = await lima_cmd(ccs=ccs, primers=primers, prefix=prefix, args=args)
    reports = {}
    for kind in ("counts", "report", "summary", "clips", "guess"):
        reports[kind] = await results.get_file(f"{prefix}.lima.{kind}")
    for kind in ("counts", "report", "summary"):
        if reports[kind] is None:
            raise FileNotFoundError(f"lima wrote no {prefix}.lima.{kind}")
    files = sorted([f async for f in results.walk()], key=lambda f: os.path.basename(f.path))
    bams = [f for f in files if f.path.endswith(".bam")]
    pbis = [f for f in files if f.path.endswith(".bam.pbi")]
    reads = [f for f in files if f.path.endswith(READ_EXTENSIONS[1:])]
    return LimaResult(bams=bams, pbis=pbis, reads=reads, results=results, **reports)
