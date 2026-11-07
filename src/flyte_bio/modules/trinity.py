"""trinity — Trinity de novo transcriptome assembly from RNA-seq reads.

Exposes :func:`trinity`, which assembles one sample from single-end reads
(``--single``) or read pairs (``--left`` / ``--right``). As upstream, the
reads are one or more files: single-end files are all passed to
``--single``; paired files alternate R1, R2, R1, R2, …, so the even ones go
to ``--left`` and the odd ones to ``--right``, in the order given.
``--seqType`` is ``fa`` when the first file is named ``*.fa``/``*.fasta``
(optionally ``.gz``), else ``fq``. Pass further Trinity options in ``args``.

Outputs keep upstream's names:

- ``<prefix>.fa.gz``: the assembled transcripts (``Trinity.fasta``);
- ``<prefix>.log``: Trinity's standard output.

Trinity's working directory (``<prefix>_trinity``) is not returned: it is
large and only holds intermediates.

``--CPU`` follows the task's CPUs and ``--max_memory`` is 80% of its memory in
whole GiB, as upstream.
"""

import re
from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

TRINITY_IMAGE = "quay.io/biocontainers/trinity:2.15.2--pl5321hdcf5f25_1"

# Upstream labels Trinity process_high and process_high_memory (12 CPUs,
# 200 GB). These defaults fit small assemblies such as the tests'; Trinity
# wants roughly 1 GB of memory per million read pairs, so raise them with
# trinity_cmd.as_task().override(resources=...) for real data.
DEFAULT_RESOURCES = flyte.Resources(cpu=4, memory="16Gi")

FASTA_NAME = re.compile(r"(.*fasta(.gz)?$)|(.*fa(.gz)?$)")

# The reads are staged together (one directory), so the wrapper passes their
# order: comma-separated basenames for --left/--single and for --right.
trinity_cmd = shell.create(
    name="trinity",
    image=TRINITY_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "reads": list[File],
        "left": str,
        "right": str,
        "seq_type": str,
        "prefix": str,
        "args": str,
    },
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        READS=({{inputs.reads}})
        LEFT={{inputs.left}}
        RIGHT={{inputs.right}}
        SEQ_TYPE={{inputs.seq_type}}
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        for F in "${{READS[@]}}"; do ln -s "$F" "$(basename "$F")"; done

        if [ -n "$RIGHT" ]; then
            READS_ARGS=(--left "$LEFT" --right "$RIGHT")
        else
            READS_ARGS=(--single "$LEFT")
        fi

        # Note that Trinity needs the word 'trinity' in the outdir
        Trinity \
            --seqType "$SEQ_TYPE" \
            --max_memory $(( TASK_MEMORY_MB * 8 / 10 / 1024 ))G \
            "${{READS_ARGS[@]}}" \
            --output "${{P}}_trinity" \
            --CPU $TASK_CPUS \
            $ARGS \
            | tee "$P.log" || exit 1

        gzip -cf "${{P}}_trinity.Trinity.fasta" > "$OUT/$P.fa.gz" || exit 1
        mv "$P.log" "$OUT"/
    """,
)

env = flyte.TaskEnvironment.from_task("trinity", trinity_cmd.as_task())


@dataclass
class TrinityResult:
    transcript_fasta: File  # <prefix>.fa.gz
    log: File  # <prefix>.log
    results: Dir  # both of the above


async def trinity(
    reads: list[File],
    single_end: bool = False,
    prefix: str = "sample",
    args: str = "",
) -> TrinityResult:
    """Assemble one sample's ``reads``: single-end files, or R1, R2[, R1, R2, …] pairs.

    The reads are staged side by side, so their names must differ.
    """
    if not reads:
        raise ValueError("trinity needs at least one read file")
    names = [r.name for r in reads]
    if len(set(names)) != len(names):
        raise ValueError(f"trinity's read files must have distinct names, got {names}")
    if any("," in n for n in names):
        raise ValueError(f"Trinity takes comma-separated read lists, so names can't contain commas: {names}")
    if single_end:
        left, right = names, []
    else:
        if len(names) % 2:
            raise ValueError(f"paired-end trinity takes R1/R2 pairs, got {len(names)} files")
        left, right = names[0::2], names[1::2]
    results = await trinity_cmd(
        reads=reads,
        left=",".join(left),
        right=",".join(right),
        seq_type="fa" if FASTA_NAME.fullmatch(left[0]) else "fq",
        prefix=prefix,
        args=args,
    )
    return await read_trinity_results(results, prefix)


async def read_trinity_results(results: Dir, prefix: str) -> TrinityResult:
    """Pick the outputs out of ``trinity_cmd``'s results Dir (e.g. from an overridden task)."""
    fasta = await results.get_file(f"{prefix}.fa.gz")
    log = await results.get_file(f"{prefix}.log")
    if fasta is None or log is None:
        raise FileNotFoundError(f"trinity wrote no {prefix}.fa.gz or {prefix}.log")
    return TrinityResult(transcript_fasta=fasta, log=log, results=results)
