"""fqtk — demultiplexing of FASTQs by inline or index-read sample barcodes.

Exposes :func:`fqtk`, which runs ``fqtk demux`` on a folder of FASTQs, given
as a ``.tar.gz`` archive (extracted into ``input/``, dropping its top-level
directory when it has just one, as upstream's untar step does), with a
sample metadata TSV (``sample_id`` and ``barcode`` columns).
``fastq_read_structures`` pairs each FASTQ's name in the archive with its
read structure (e.g. ``("x_I1.fastq.gz", "8B")``, ``("x_R1.fastq.gz",
"150T")``), in the order fqtk should read them.

fqtk writes one FASTQ per sample and read, which can be many, so they come
back as one archive. The output Dir holds:

- ``fastq.tar.gz``: ``<sample_id>.R1.fq.gz``, ``<sample_id>.R2.fq.gz``, ...
  (names follow ``--output-types``, template reads by default), with the
  unmatched reads as ``unmatched.R1.fq.gz``, ...;
- ``demux-metrics.txt``: per-sample read counts.

fqtk's ``--threads`` is the task's CPUs, at least 3 (fqtk's minimum).
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

FQTK_IMAGE = "community.wave.seqera.io/library/fqtk_gzip:8f4e4f93941dcf35"
# Upstream's process_high; fqtk's default is 8 threads.
DEFAULT_RESOURCES = flyte.Resources(cpu=8, memory="16Gi")

fqtk_cmd = shell.create(
    name="fqtk",
    image=FQTK_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"sample_sheet": File, "fastqs": File, "fastq_names": str, "read_structures": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        set -e
        SHEET=({{inputs.sample_sheet}})
        ARCHIVE=({{inputs.fastqs}})
        NAMES={{inputs.fastq_names}}
        STRUCTURES={{inputs.read_structures}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        THREADS=$TASK_CPUS
        if [ $THREADS -lt 3 ]; then THREADS=3; fi
        W=$(mktemp -d); cd "$W"
        mkdir input output
        # As upstream: drop the top-level directory only when the archive has just one.
        if [[ $(tar -tzf "${{ARCHIVE[0]}}" | grep / | cut -d/ -f1 | sort -u | wc -l) -eq 1 ]]; then
            tar -C input --strip-components 1 -xzf "${{ARCHIVE[0]}}"
        else
            tar -C input -xzf "${{ARCHIVE[0]}}"
        fi
        FASTQS=()
        for name in $NAMES; do FASTQS+=("input/$name"); done
        # fqtk is a Rust tool: a panic exits 255, which the platform would retry forever.
        fqtk demux \
            --inputs "${{FASTQS[@]}}" \
            --read-structures $STRUCTURES \
            --output output/ \
            --sample-metadata "${{SHEET[0]}}" \
            --threads $THREADS \
            $ARGS || exit 1
        (cd output && ls | grep '\.fq\.gz$' | tar -czf "$OUT/fastq.tar.gz" -T -)
        cp output/demux-metrics.txt "$OUT"/
    """,
)

env = flyte.TaskEnvironment.from_task("fqtk", fqtk_cmd.as_task())


@dataclass
class FqtkResult:
    fastq: File  # fastq.tar.gz: per-sample FASTQs, and the unmatched reads
    metrics: File  # demux-metrics.txt
    results: Dir  # both


async def fqtk(
    sample_sheet: File, fastqs: File, fastq_read_structures: list[tuple[str, str]], args: str = ""
) -> FqtkResult:
    """Demultiplex the FASTQs in archive ``fastqs`` by the samples' barcodes in ``sample_sheet``."""
    results = await fqtk_cmd(
        sample_sheet=sample_sheet,
        fastqs=fastqs,
        fastq_names=" ".join(name for name, _ in fastq_read_structures),
        read_structures=" ".join(structure for _, structure in fastq_read_structures),
        args=args,
    )

    async def pick(name: str) -> File:
        found = await results.get_file(name)
        if found is None:
            raise FileNotFoundError(f"fqtk wrote no {name}")
        return found

    return FqtkResult(fastq=await pick("fastq.tar.gz"), metrics=await pick("demux-metrics.txt"), results=results)
