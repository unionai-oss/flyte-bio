"""mirtrace — miRTrace quality control of small-RNA sequencing reads.

Exposes :func:`mirtrace_qc`, ``mirtrace qc`` over one or more FASTQ files for
one species (miRBase code, e.g. ``hsa``), writing, as upstream:

- ``mirtrace-report.html`` and ``mirtrace-results.json``;
- ``mirtrace-stats-*.tsv`` (contamination, length, complexity, phred,
  QC status, RNA type);
- the QC-passed reads, collapsed, as ``qc_passed_reads.all.collapsed/*.fasta``
  and ``qc_passed_reads.rnatype_unknown.collapsed/*.fasta``.

With a ``config`` (miRTrace's CSV of ``path,name,adapter,phred``), the reads
are listed by the config instead of the command line. The reads are linked
into the working directory under their names first, so the config's relative
paths (``./reads.fastq.gz``) resolve, whether miRTrace reads them against the
working directory or the config's.

The JVM is started directly rather than through the bundled launcher, which
sizes the heap from the host's memory and ignores the tool's exit status. As
upstream's wrapper does, it is told its wrapper name (``mirtrace``) for the
report; the heap is 80% of the task's memory.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import JAVA_HEAP_MB, TASK_LIMITS

MIRTRACE_IMAGE = "quay.io/biocontainers/mirtrace:1.0.1--0"
DEFAULT_RESOURCES = flyte.Resources(cpu=4, memory="8Gi")

mirtrace_qc_cmd = shell.create(
    name="mirtrace_qc",
    image=MIRTRACE_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"reads": list[File], "names": str, "config": File | None, "species": str},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        READS=({{inputs.reads}})
        shopt -s nullglob; CONFIG=({{inputs.config}}); shopt -u nullglob
        NAMES={{inputs.names}}
        SPECIES={{inputs.species}}
        OUT={{outputs.results}}
        D=$(dirname "${{READS[0]}}")
        W=$(mktemp -d); cd "$W"
        FILES=()
        while IFS= read -r NAME; do ln -s "$D/$NAME" "$NAME"; FILES+=("$NAME"); done <<< "$NAMES"
        if [ ${{#CONFIG[@]}} -gt 0 ]; then
            cp "${{CONFIG[0]}}" mirtrace_config.txt
            MODE=(--config mirtrace_config.txt)
        else
            MODE=("${{FILES[@]}}")
        fi

        JAR=$(dirname "$(readlink -f "$(command -v mirtrace)")")/mirtrace.jar
        mkdir results
        java -Xmx{JAVA_HEAP_MB}m -jar "$JAR" --mirtrace-wrapper-name mirtrace qc \
            --species "$SPECIES" --write-fasta --output-dir results --force "${{MODE[@]}}"
        mv results/* "$OUT"/
    """,
)

env = flyte.TaskEnvironment.from_task("mirtrace", mirtrace_qc_cmd.as_task())


def file_name(file: File) -> str:
    return file.path.rstrip("/").rsplit("/", 1)[-1]


@dataclass
class MirtraceQcResult:
    html: File  # mirtrace-report.html
    json: File  # mirtrace-results.json
    tsv: list[File]  # mirtrace-stats-*.tsv, by name
    all_fasta: list[File]  # qc_passed_reads.all.collapsed/*.fasta
    rnatype_unknown_fasta: list[File]  # qc_passed_reads.rnatype_unknown.collapsed/*.fasta
    results: Dir


async def mirtrace_qc(reads: list[File], species: str, config: File | None = None) -> MirtraceQcResult:
    """miRTrace QC of ``reads`` for ``species``; with ``config``, the config lists the samples."""
    names = [file_name(r) for r in reads]
    if not names:
        raise ValueError("mirtrace qc needs at least one reads file")
    if len(set(names)) != len(names):
        raise ValueError(f"mirtrace qc inputs need distinct names, got {names}")
    results = await mirtrace_qc_cmd(reads=reads, names="\n".join(names), config=config, species=species)
    files = sorted([f async for f in results.walk()], key=lambda f: f.path)

    def under(folder: str) -> list[File]:
        return [f for f in files if f"/{folder}/" in f.path and f.path.endswith((".fa", ".fasta"))]

    def one(name: str) -> File:
        found = [f for f in files if file_name(f) == name]
        if len(found) != 1:
            raise FileNotFoundError(f"mirtrace qc wrote no {name}")
        return found[0]

    return MirtraceQcResult(
        html=one("mirtrace-report.html"),
        json=one("mirtrace-results.json"),
        tsv=[f for f in files if file_name(f).endswith(".tsv")],
        all_fasta=under("qc_passed_reads.all.collapsed"),
        rnatype_unknown_fasta=under("qc_passed_reads.rnatype_unknown.collapsed"),
        results=results,
    )
