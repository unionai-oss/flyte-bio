"""kraken2 — taxonomic classification of reads against a Kraken2 database.

Exposes :func:`kraken2`, which classifies one sample's reads (single-end, or
paired-end with ``reads_2``) against ``db``, a Dir holding the database's
``hash.k2d``, ``opts.k2d`` and ``taxo.k2d``. Outputs keep upstream's names:

- the report: ``<prefix>.kraken2.report.txt`` (always);
- with ``save_output_fastqs``: the classified and unclassified reads,
  ``<prefix>.classified.fastq.gz`` / ``<prefix>.unclassified.fastq.gz``
  (single-end) or ``<prefix>.classified_{1,2}.fastq.gz`` /
  ``<prefix>.unclassified_{1,2}.fastq.gz`` (paired-end), compressed with pigz;
- with ``save_reads_assignment``: the per-read classifications,
  ``<prefix>.kraken2.classifiedreads.txt``.

kraken2 loads the whole database into memory, so a real database (the
standard one is tens of GiB) needs ``memory`` raised with
``kraken2_cmd.as_task().override(resources=...)``; the default fits small
databases like the test one.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

KRAKEN2_IMAGE = "community.wave.seqera.io/library/kraken2_coreutils_pigz:920ecc6b96e2ba71"
DEFAULT_RESOURCES = flyte.Resources(cpu=4, memory="8Gi")

kraken2_cmd = shell.create(
    name="kraken2",
    image=KRAKEN2_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "reads_1": File,
        "reads_2": File | None,
        "db": Dir,
        "prefix": str,
        "save_output_fastqs": bool,
        "save_reads_assignment": bool,
        "args": str,
    },
    defaults={"save_output_fastqs": False, "save_reads_assignment": False, "args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        R1=({{inputs.reads_1}})
        shopt -s nullglob; R2=({{inputs.reads_2}}); shopt -u nullglob
        DB={{inputs.db}}
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"

        READS=("${{R1[0]}}")
        EXTRA=()
        if [ ${{#R2[@]}} -gt 0 ]; then
            READS+=("${{R2[0]}}")
            EXTRA+=(--paired)
            CLASSIFIED="$P.classified#.fastq"
            UNCLASSIFIED="$P.unclassified#.fastq"
        else
            CLASSIFIED="$P.classified.fastq"
            UNCLASSIFIED="$P.unclassified.fastq"
        fi
        SAVE_FASTQS={{inputs.save_output_fastqs}}
        if [ $SAVE_FASTQS = true ]; then
            EXTRA+=(--classified-out "$CLASSIFIED" --unclassified-out "$UNCLASSIFIED")
        fi
        if [ {{inputs.save_reads_assignment}} = true ]; then
            EXTRA+=(--output "$P.kraken2.classifiedreads.txt")
        else
            EXTRA+=(--output /dev/null)
        fi

        kraken2 \
            --db "$DB" \
            --threads $TASK_CPUS \
            --report "$P.kraken2.report.txt" \
            --gzip-compressed \
            "${{EXTRA[@]}}" \
            $ARGS \
            "${{READS[@]}}"

        if [ $SAVE_FASTQS = true ]; then pigz -p $TASK_CPUS ./*.fastq; fi
        shopt -s nullglob
        mv ./*.kraken2.report.txt ./*.classifiedreads.txt ./*.fastq.gz "$OUT"/
    """,
)

env = flyte.TaskEnvironment.from_task("kraken2", kraken2_cmd.as_task())


@dataclass
class Kraken2Result:
    report: File
    classified_reads_fastq: list[File]  # with save_output_fastqs: one file, or mates 1 and 2
    unclassified_reads_fastq: list[File]  # with save_output_fastqs: one file, or mates 1 and 2
    classified_reads_assignment: File | None  # with save_reads_assignment
    results: Dir  # everything kraken2 wrote


async def kraken2(
    reads_1: File,
    db: Dir,
    reads_2: File | None = None,
    prefix: str = "sample",
    save_output_fastqs: bool = False,
    save_reads_assignment: bool = False,
    args: str = "",
) -> Kraken2Result:
    """Classify one sample's reads against ``db``; ``reads_2`` is None for single-end reads."""
    results = await kraken2_cmd(
        reads_1=reads_1,
        reads_2=reads_2,
        db=db,
        prefix=prefix,
        save_output_fastqs=save_output_fastqs,
        save_reads_assignment=save_reads_assignment,
        args=args,
    )

    async def pick(name: str) -> File:
        found = await results.get_file(name)
        if found is None:
            raise FileNotFoundError(f"kraken2 wrote no {name}")
        return found

    mates = ["_1", "_2"] if reads_2 is not None else [""]
    classified = []
    unclassified = []
    if save_output_fastqs:
        classified = [await pick(f"{prefix}.classified{mate}.fastq.gz") for mate in mates]
        unclassified = [await pick(f"{prefix}.unclassified{mate}.fastq.gz") for mate in mates]
    assignment = None
    if save_reads_assignment:
        assignment = await pick(f"{prefix}.kraken2.classifiedreads.txt")
    return Kraken2Result(
        report=await pick(f"{prefix}.kraken2.report.txt"),
        classified_reads_fastq=classified,
        unclassified_reads_fastq=unclassified,
        classified_reads_assignment=assignment,
        results=results,
    )
