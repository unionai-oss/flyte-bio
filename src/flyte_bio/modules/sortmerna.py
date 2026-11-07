"""sortmerna — filter ribosomal RNA reads out of FASTQ reads.

Exposes :func:`sortmerna`, which aligns one sample's reads (single- or
paired-end) against rRNA reference FASTAs and keeps the reads that don't
match. As upstream, ``args`` selects the indexing mode:

- ``--index 1``: only build the index (no reads);
- ``--index 0``: use the premade ``index`` (a Dir from an earlier run);
- neither: build the index, then filter.

Outputs keep upstream's names: ``<prefix>.non_rRNA.fastq.gz`` (single-end) or
``<prefix>_1.non_rRNA.fastq.gz`` / ``<prefix>_2.non_rRNA.fastq.gz``, the log
``<prefix>.sortmerna.log``, and the index directory ``idx`` (a few files per
reference, named by a hash of the reference).
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

SORTMERNA_IMAGE = "community.wave.seqera.io/library/sortmerna:4.3.7--b730cad73fc42b8e"
# Upstream labels it process_high; the full rRNA databases index in about 16 GiB.
DEFAULT_RESOURCES = flyte.Resources(cpu=4, memory="16Gi")

# The references are linked into the working directory and named relatively,
# as upstream stages them. A premade index is copied there as `idx`, where
# upstream stages it, and returned again as the index.
sortmerna_cmd = shell.create(
    name="sortmerna",
    image=SORTMERNA_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "reads_1": File | None,
        "reads_2": File | None,
        "fastas": list[File],
        "index": Dir | None,
        "prefix": str,
        "args": str,
    },
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        shopt -s nullglob; R1=({{inputs.reads_1}}); R2=({{inputs.reads_2}}); shopt -u nullglob
        FASTAS=({{inputs.fastas}})
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"

        REFS=()
        for F in "${{FASTAS[@]}}"; do
            ln -s "$F" "$(basename "$F")"
            REFS+=(--ref "$(basename "$F")")
        done
        INDEX_ONLY=false; SKIP_INDEX=false
        if [[ "$ARGS" == *"--index 1"* ]]; then INDEX_ONLY=true; fi
        if [[ "$ARGS" == *"--index 0"* ]]; then SKIP_INDEX=true; fi
        IDX=()
        INDEX_IN={{inputs.index}}
        if [ $SKIP_INDEX = true ] && [ -d "$INDEX_IN" ] && [ -n "$(ls -A "$INDEX_IN")" ]; then
            cp -rL "$INDEX_IN" idx
            chmod -R u+w idx
            IDX=(--idx-dir idx)
        fi
        READS=()
        if [ $INDEX_ONLY = false ]; then
            READS=(--reads "${{R1[0]}}")
            if [ ${{#R2[@]}} -gt 0 ]; then READS+=(--reads "${{R2[0]}}"); fi
            READS+=(--aligned rRNA_reads --fastx --other non_rRNA_reads)
            if [ ${{#R2[@]}} -gt 0 ]; then READS+=(--paired_in --out2); fi
        fi

        sortmerna "${{REFS[@]}}" "${{IDX[@]}}" "${{READS[@]}}" --threads $TASK_CPUS --workdir . $ARGS

        if [ $INDEX_ONLY = false ]; then
            if [ ${{#R2[@]}} -gt 0 ]; then
                mv non_rRNA_reads_fwd.f*q.gz "$OUT/${{P}}_1.non_rRNA.fastq.gz"
                mv non_rRNA_reads_rev.f*q.gz "$OUT/${{P}}_2.non_rRNA.fastq.gz"
            else
                mv non_rRNA_reads.f*q.gz "$OUT/$P.non_rRNA.fastq.gz"
            fi
            mv rRNA_reads.log "$OUT/$P.sortmerna.log"
        fi
        mv idx "$OUT"/idx
    """,
)

env = flyte.TaskEnvironment.from_task("sortmerna", sortmerna_cmd.as_task())


@dataclass
class SortmernaResult:
    reads: list[File]  # the non-rRNA reads: one file, or mates 1 and 2; empty with --index 1
    log: File | None  # <prefix>.sortmerna.log; None with --index 1
    index: Dir  # idx: the index sortmerna built or was given


async def sortmerna(
    fastas: list[File],
    reads_1: File | None = None,
    reads_2: File | None = None,
    index: Dir | None = None,
    prefix: str = "sample",
    args: str = "",
) -> SortmernaResult:
    """Remove rRNA reads from one sample; ``reads_2`` is None for single-end."""
    results = await sortmerna_cmd(
        reads_1=reads_1, reads_2=reads_2, fastas=fastas, index=index, prefix=prefix, args=args
    )

    async def pick(name: str) -> File:
        found = await results.get_file(name)
        if found is None:
            raise FileNotFoundError(f"sortmerna wrote no {name}")
        return found

    reads = []
    log = None
    if "--index 1" not in args:
        names = [f"{prefix}_1.non_rRNA.fastq.gz", f"{prefix}_2.non_rRNA.fastq.gz"] if reads_2 else []
        names = names or [f"{prefix}.non_rRNA.fastq.gz"]
        reads = [await pick(name) for name in names]
        log = await pick(f"{prefix}.sortmerna.log")
    return SortmernaResult(reads=reads, log=log, index=Dir.from_existing_remote(results.path.rstrip("/") + "/idx"))
