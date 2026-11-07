"""seqkit — FASTA/FASTQ statistics and sequence editing.

Exposes:

- :func:`seqkit_stats` — ``seqkit stats --tabular`` over one or more
  FASTA/FASTQ files (``--all`` unless ``args`` says otherwise), written to
  ``<prefix>.tsv``. As upstream, the files are passed under their own names,
  which seqkit records in the table's ``file`` column, in the caller's order.
- :func:`seqkit_replace` — ``seqkit replace`` over one FASTA/FASTQ file, with
  the pattern and replacement in ``args`` (e.g. ``-s -p A -r N``). The output
  is ``<prefix>.<out_ext>``; ``out_ext`` defaults, as upstream, to ``fasta``
  for FASTA-named inputs and ``fastq`` otherwise, plus ``.gz`` for a gzipped
  input. seqkit compresses by the output's extension.
- :func:`seqkit_seq` — ``seqkit seq`` over one FASTA/FASTQ file, with the
  transformation in ``args`` (e.g. ``--min-len 100``). As upstream, the output
  is ``<prefix>.fasta`` for a FASTA-named input and ``<prefix>.fastq``
  otherwise, gzipped by ``gzip -c <gzip_args>`` (``<prefix>.<ext>.gz``) when
  the input is; a prefix that would overwrite the input is an error.
- :func:`seqkit_rmdup` — ``seqkit rmdup`` over one FASTA/FASTQ file (by
  sequence ID unless ``args`` says ``-s`` or ``-n``). As upstream, the output is
  ``<prefix>.fasta`` or ``<prefix>.fastq`` by the input's name, plus ``.gz``
  (compressed by seqkit) for a gzipped input, with seqkit's messages (``[INFO]
  N duplicated records removed``) in ``<prefix>.log``; a prefix that would
  overwrite the input is an error.
- :func:`seqkit_split2` — ``seqkit split2`` of single-end reads, or of a pair
  (``--read1``/``--read2``), into parts by ``args`` (``--by-part 3``,
  ``--by-size 50``, ``--by-length 8K``). The parts keep seqkit's names
  (``<input name>.part_001.fastq.gz``, …) and are returned in name order. The
  parts come back as a Dir, which stages up to a few hundred files: split into
  fewer, bigger parts (``--by-part``) for larger inputs.

All pass the task's CPUs as ``--threads``. ``args`` is split on whitespace
without shell quote removal, so write ``-p A``, not ``-p 'A'``; patterns are
not globbed, so regex characters (``*``, ``?``, ``[``) are safe, but a pattern
can't contain spaces.
"""

import re
from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

SEQKIT_IMAGE = "community.wave.seqera.io/library/seqkit:2.13.0--3bacfb68370dd835"
DEFAULT_RESOURCES = flyte.Resources(cpu=2, memory="4Gi")

# Input names that upstream treats as FASTA: seqkit replace's and seqkit seq's lists differ.
FASTA_NAME = re.compile(r".+\.(fasta|fa|fas|fna|faa)(\.gz)?")
SEQ_FASTA_NAME = re.compile(r".+\.(fasta|fa|fas|fna|fsa)(\.gz)?")

# A list[File] is staged under the files' names and globbed alphabetically, so
# ``names`` carries the caller's order (one name per line) into the script.
seqkit_stats_cmd = shell.create(
    name="seqkit_stats",
    image=SEQKIT_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"reads": list[File], "names": str, "prefix": str, "args": str},
    defaults={"args": "--all"},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        IN=({{inputs.reads}})
        NAMES={{inputs.names}}
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        set -f
        W=$(mktemp -d); cd "$W"
        D=$(dirname "${{IN[0]}}")
        FILES=()
        while IFS= read -r NAME; do
            ln -s "$D/$NAME" "$NAME"
            FILES+=("$NAME")
        done <<< "$NAMES"
        seqkit stats --tabular --threads $TASK_CPUS $ARGS "${{FILES[@]}}" > "$OUT/$P.tsv"
    """,
)

seqkit_replace_cmd = shell.create(
    name="seqkit_replace",
    image=SEQKIT_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"fastx": File, "filename": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        IN=({{inputs.fastx}})
        NAME={{inputs.filename}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        set -f
        seqkit replace $ARGS --threads $TASK_CPUS -i "${{IN[0]}}" -o "$OUT/$NAME"
    """,
)

seqkit_seq_cmd = shell.create(
    name="seqkit_seq",
    image=SEQKIT_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"fastx": File, "filename": str, "args": str, "gzip_args": str},
    defaults={"args": "", "gzip_args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        set -o pipefail
        IN=({{inputs.fastx}})
        NAME={{inputs.filename}}
        ARGS={{inputs.args}}
        GZIP_ARGS={{inputs.gzip_args}}
        OUT={{outputs.results}}
        set -f
        if [[ $NAME == *.gz ]]; then
            seqkit seq --threads $TASK_CPUS $ARGS "${{IN[0]}}" | gzip -c $GZIP_ARGS > "$OUT/$NAME" || exit 1
        else
            seqkit seq --threads $TASK_CPUS $ARGS "${{IN[0]}}" > "$OUT/$NAME" || exit 1
        fi
    """,
)

seqkit_rmdup_cmd = shell.create(
    name="seqkit_rmdup",
    image=SEQKIT_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"fastx": File, "filename": str, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        IN=({{inputs.fastx}})
        NAME={{inputs.filename}}
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        set -f
        STATUS=0
        seqkit rmdup --threads $TASK_CPUS $ARGS "${{IN[0]}}" -o "$OUT/$NAME" 2> "$OUT/$P.log" || STATUS=1
        cat "$OUT/$P.log" >&2
        exit $STATUS
    """,
)

# The reads are separate inputs so a pair keeps its order (a list[File] is staged alphabetically).
seqkit_split2_cmd = shell.create(
    name="seqkit_split2",
    image=SEQKIT_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"read1": File, "read2": File | None, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        shopt -s nullglob; R1=({{inputs.read1}}); R2=({{inputs.read2}}); shopt -u nullglob
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        set -f
        W=$(mktemp -d)
        if [ ${{#R2[@]}} -gt 0 ]; then INPUT=(--read1 "${{R1[0]}}" --read2 "${{R2[0]}}"); else INPUT=("${{R1[0]}}"); fi
        seqkit split2 $ARGS --threads $TASK_CPUS "${{INPUT[@]}}" --out-dir "$W/parts" || exit 1
        set +f
        mv "$W/parts"/* "$OUT"/
    """,
)

env = flyte.TaskEnvironment.from_task(
    "seqkit",
    seqkit_stats_cmd.as_task(),
    seqkit_replace_cmd.as_task(),
    seqkit_seq_cmd.as_task(),
    seqkit_rmdup_cmd.as_task(),
    seqkit_split2_cmd.as_task(),
)


@dataclass
class SeqkitRmdupResult:
    fastx: File  # <prefix>.fasta / <prefix>.fastq, plus .gz for a gzipped input
    log: File  # <prefix>.log: seqkit's messages, e.g. "N duplicated records removed"


def file_name(file: File) -> str:
    return file.path.rstrip("/").rsplit("/", 1)[-1]


async def pick(results: Dir, name: str, tool: str) -> File:
    found = await results.get_file(name)
    if found is None:
        raise FileNotFoundError(f"{tool} wrote no {name}")
    return found


async def seqkit_stats(reads: list[File], prefix: str = "sample", args: str = "--all") -> File:
    """``seqkit stats`` of ``reads`` (one row per file, in list order): ``<prefix>.tsv``.

    An empty ``args`` means ``--all``, as upstream.
    """
    if not reads:
        raise ValueError("seqkit_stats needs at least one file")
    names = [file_name(f) for f in reads]
    if len(set(names)) != len(names):
        raise ValueError(f"seqkit_stats inputs need distinct names, got {names}")
    results = await seqkit_stats_cmd(reads=reads, names="\n".join(names), prefix=prefix, args=args or "--all")
    return await pick(results, f"{prefix}.tsv", "seqkit stats")


def replace_extension(name: str) -> str:
    """Upstream's default output extension for an input named ``name``."""
    extension = "fasta" if FASTA_NAME.fullmatch(name) else "fastq"
    return extension + (".gz" if name.endswith(".gz") else "")


async def seqkit_replace(fastx: File, prefix: str = "sample", out_ext: str = "", args: str = "") -> File:
    """``seqkit replace`` over ``fastx``: ``<prefix>.<out_ext>``, gzipped when ``out_ext`` ends in ``.gz``."""
    filename = f"{prefix}.{out_ext or replace_extension(file_name(fastx))}"
    results = await seqkit_replace_cmd(fastx=fastx, filename=filename, args=args)
    return await pick(results, filename, "seqkit replace")


async def seqkit_seq(fastx: File, prefix: str = "sample", args: str = "", gzip_args: str = "") -> File:
    """``seqkit seq`` over ``fastx``: ``<prefix>.fasta`` or ``<prefix>.fastq``, plus ``.gz`` for a gzipped input."""
    name = file_name(fastx)
    extension = "fasta" if SEQ_FASTA_NAME.fullmatch(name) else "fastq"
    filename = f"{prefix}.{extension}" + (".gz" if name.endswith(".gz") else "")
    if filename == name:
        raise ValueError("Input and output names are the same, set prefix to disambiguate!")
    results = await seqkit_seq_cmd(fastx=fastx, filename=filename, args=args, gzip_args=gzip_args)
    return await pick(results, filename, "seqkit seq")


async def seqkit_rmdup(fastx: File, prefix: str = "sample", args: str = "") -> SeqkitRmdupResult:
    """``seqkit rmdup`` over ``fastx``: the deduplicated reads and seqkit's log."""
    name = file_name(fastx)
    extension = "fasta" if SEQ_FASTA_NAME.fullmatch(name) else "fastq"
    filename = f"{prefix}.{extension}" + (".gz" if name.endswith(".gz") else "")
    if filename == name:
        raise ValueError('Input and output names are the same, use "prefix" to disambiguate!')
    results = await seqkit_rmdup_cmd(fastx=fastx, filename=filename, prefix=prefix, args=args)
    return SeqkitRmdupResult(
        fastx=await pick(results, filename, "seqkit rmdup"),
        log=await pick(results, f"{prefix}.log", "seqkit rmdup"),
    )


async def seqkit_split2(reads: list[File], args: str = "") -> list[File]:
    """``seqkit split2`` of one file (single-end) or a pair: the parts, in name order."""
    if len(reads) not in (1, 2):
        raise ValueError(f"seqkit_split2 takes one file or a pair, got {len(reads)}")
    read2 = reads[1] if len(reads) == 2 else None
    results = await seqkit_split2_cmd(read1=reads[0], read2=read2, args=args)
    return sorted(await results.list_files(), key=file_name)
