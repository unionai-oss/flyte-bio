"""seqtk — subsample and quality-trim FASTQ reads.

Exposes :func:`seqtk_sample`, which samples ``sample_size`` reads (a count,
or a fraction below 1) from each FASTQ file with ``seqtk sample``. As
upstream:

- each input ``<name>`` gives ``<prefix>_<name>``, with a ``.fq.gz`` suffix
  spelled ``.fastq.gz``, always gzipped (``gzip --no-name``);
- a seed is added (``-s100``) unless ``args`` sets one (`` -s<n>``), so
  paired files sampled with one seed stay paired.

and :func:`seqtk_trim`, which quality-trims each FASTQ file with ``seqtk
trimfq`` (its default: Phred error-probability trimming at 0.05; ``args``
passes ``-q``, ``-b``, ``-e`` and the rest). As upstream, each input
``<name>`` gives ``<prefix>_<name>``, gzipped (``gzip --no-name``).

:func:`seqtk_seq` converts or transforms one FASTA/FASTQ file with ``seqtk
seq`` (``args``: ``-A`` for FASTA output, ``-l``, ``-q`` and the rest). As
upstream, the output is ``<prefix>.seqtk-seq.<ext>.gz`` (``gzip -c``), where
``ext`` is ``fasta`` for a FASTA-named input or when ``args`` is exactly
``-a`` or ``-A``, and ``fastq`` otherwise.

:func:`seqtk_cutn` writes the runs of N (gaps) in a FASTA as a BED,
``seqtk cutN -g`` (``-n`` sets the shortest gap, in ``args``), to
``<prefix>.bed``.

:func:`seqtk_subseq` extracts the sequences (or, with a BED, the regions)
named in ``filter_list`` from a FASTA/FASTQ file, ``seqtk subseq``. As
upstream, the output is ``<prefix>.<ext>.gz`` (``gzip --no-name``), where
``ext`` is ``fq`` for a FASTQ-named input and ``fa`` otherwise.

:func:`seqtk_mergepe` interleaves paired-end FASTQ files into one,
``seqtk mergepe``: ``<prefix>.fastq.gz`` (``gzip -n``). As upstream, a single
(single-end) file is passed through unchanged under that name.
"""

import re

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

SEQTK_IMAGE = "community.wave.seqera.io/library/seqtk:1.4--3fa479d3961a5b8d"
# Upstream's process_single.
DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory="2Gi")

# Upstream's test for a seed already in ``args``.
HAS_SEED = re.compile(r".* -s ?[0-9]+.*")

# A list[File] is staged under the files' names and globbed alphabetically, so
# ``names`` carries the caller's order (one name per line) into the script.
seqtk_sample_cmd = shell.create(
    name="seqtk_sample",
    image=SEQTK_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"reads": list[File], "names": str, "sample_size": str, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        set -o pipefail
        IN=({inputs.reads})
        NAMES={inputs.names}
        SIZE={inputs.sample_size}
        P={inputs.prefix}
        ARGS={inputs.args}
        OUT={outputs.results}
        set -f
        D=$(dirname "${IN[0]}")
        while IFS= read -r NAME; do
            NEW=$(printf '%s\n' "$NAME" | sed 's/\.f\(ast\)\?q\.gz$/.fastq.gz/')
            seqtk sample $ARGS "$D/$NAME" "$SIZE" | gzip --no-name > "$OUT/${P}_$NEW" || exit 1
        done <<< "$NAMES"
    """,
)

# Upstream's process_low asks for 2 CPUs and 12 GB; seqtk trimfq is
# single-threaded and streams, so 2Gi is plenty.
seqtk_trim_cmd = shell.create(
    name="seqtk_trim",
    image=SEQTK_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"reads": list[File], "names": str, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        set -o pipefail
        IN=({inputs.reads})
        NAMES={inputs.names}
        P={inputs.prefix}
        ARGS={inputs.args}
        OUT={outputs.results}
        D=$(dirname "${IN[0]}")
        while IFS= read -r NAME; do
            seqtk trimfq $ARGS "$D/$NAME" | gzip --no-name > "$OUT/${P}_$NAME" || exit 1
        done <<< "$NAMES"
    """,
)

# Upstream's test for a FASTA input (or FASTA output) to name seqtk seq's output.
SEQ_FASTA_NAME = re.compile(r".+\.fasta|.+\.fasta.gz|.+\.fa|.+\.fa.gz|.+\.fas|.+\.fas.gz|.+\.fna|.+\.fna.gz")
SEQ_FASTA_ARGS = re.compile(r"\-[aA]")

seqtk_seq_cmd = shell.create(
    name="seqtk_seq",
    image=SEQTK_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"fastx": File, "filename": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        set -o pipefail
        IN=({inputs.fastx})
        NAME={inputs.filename}
        ARGS={inputs.args}
        OUT={outputs.results}
        set -f
        seqtk seq $ARGS "${IN[0]}" | gzip -c > "$OUT/$NAME" || exit 1
    """,
)

# Upstream's test for a FASTQ input, to name seqtk subseq's output.
SUBSEQ_FASTQ_NAME = re.compile(r".+\.fq|.+\.fq.gz|.+\.fastq|.+\.fastq.gz")

# Upstream labels cutN process_low; it is single-threaded and streams.
seqtk_cutn_cmd = shell.create(
    name="seqtk_cutn",
    image=SEQTK_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"fasta": File, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        IN=({inputs.fasta})
        P={inputs.prefix}
        ARGS={inputs.args}
        OUT={outputs.results}
        set -f
        seqtk cutN $ARGS -g "${IN[0]}" > "$OUT/$P.bed" || exit 1
    """,
)

seqtk_subseq_cmd = shell.create(
    name="seqtk_subseq",
    image=SEQTK_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"sequences": File, "filter_list": File, "filename": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        set -o pipefail
        IN=({inputs.sequences})
        LIST=({inputs.filter_list})
        NAME={inputs.filename}
        ARGS={inputs.args}
        OUT={outputs.results}
        set -f
        seqtk subseq $ARGS "${IN[0]}" "${LIST[0]}" | gzip --no-name > "$OUT/$NAME" || exit 1
    """,
)

seqtk_mergepe_cmd = shell.create(
    name="seqtk_mergepe",
    image=SEQTK_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"reads_1": File, "reads_2": File | None, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        set -o pipefail
        R1=({inputs.reads_1})
        shopt -s nullglob; R2=({inputs.reads_2}); shopt -u nullglob
        P={inputs.prefix}
        ARGS={inputs.args}
        OUT={outputs.results}
        if [ ${#R2[@]} -eq 0 ]; then
            cp "${R1[0]}" "$OUT/$P.fastq.gz"
        else
            seqtk mergepe $ARGS "${R1[0]}" "${R2[0]}" | gzip -n > "$OUT/$P.fastq.gz" || exit 1
        fi
    """,
)

env = flyte.TaskEnvironment.from_task(
    "seqtk",
    seqtk_sample_cmd.as_task(),
    seqtk_trim_cmd.as_task(),
    seqtk_seq_cmd.as_task(),
    seqtk_cutn_cmd.as_task(),
    seqtk_subseq_cmd.as_task(),
    seqtk_mergepe_cmd.as_task(),
)


def file_name(file: File) -> str:
    return file.path.rstrip("/").rsplit("/", 1)[-1]


def sampled_name(name: str) -> str:
    return re.sub(r"\.f(ast)?q\.gz$", ".fastq.gz", name)


async def seqtk_sample(
    reads: list[File], sample_size: int | float, prefix: str = "sample", args: str = ""
) -> list[File]:
    """``sample_size`` reads sampled from each of ``reads``: ``<prefix>_<name>``, in list order."""
    if not sample_size:
        raise ValueError("SEQTK/SAMPLE must have a sample_size value included")
    if not reads:
        raise ValueError("seqtk_sample needs at least one file")
    names = [file_name(f) for f in reads]
    if len(set(names)) != len(names):
        raise ValueError(f"seqtk_sample inputs need distinct names, got {names}")
    if not HAS_SEED.fullmatch(args):
        args += " -s100"
    results = await seqtk_sample_cmd(
        reads=reads, names="\n".join(names), sample_size=str(sample_size), prefix=prefix, args=args.strip()
    )
    outputs = []
    for name in names:
        out = f"{prefix}_{sampled_name(name)}"
        found = await results.get_file(out)
        if found is None:
            raise FileNotFoundError(f"seqtk sample wrote no {out}")
        outputs.append(found)
    return outputs


async def seqtk_trim(reads: list[File], prefix: str = "sample", args: str = "") -> list[File]:
    """Each of ``reads`` quality-trimmed with ``seqtk trimfq``: ``<prefix>_<name>``, gzipped, in list order."""
    if not reads:
        raise ValueError("seqtk_trim needs at least one file")
    names = [file_name(f) for f in reads]
    if len(set(names)) != len(names):
        raise ValueError(f"seqtk_trim inputs need distinct names, got {names}")
    results = await seqtk_trim_cmd(reads=reads, names="\n".join(names), prefix=prefix, args=args)
    outputs = []
    for name in names:
        found = await results.get_file(f"{prefix}_{name}")
        if found is None:
            raise FileNotFoundError(f"seqtk trimfq wrote no {prefix}_{name}")
        outputs.append(found)
    return outputs


async def pick(results: Dir, name: str, tool: str) -> File:
    found = await results.get_file(name)
    if found is None:
        raise FileNotFoundError(f"{tool} wrote no {name}")
    return found


async def seqtk_seq(fastx: File, prefix: str = "sample", args: str = "") -> File:
    """``seqtk seq`` over ``fastx``: ``<prefix>.seqtk-seq.fasta.gz`` or ``<prefix>.seqtk-seq.fastq.gz``."""
    fasta = SEQ_FASTA_NAME.fullmatch(file_name(fastx)) or SEQ_FASTA_ARGS.fullmatch(args)
    filename = f"{prefix}.seqtk-seq.{'fasta' if fasta else 'fastq'}.gz"
    results = await seqtk_seq_cmd(fastx=fastx, filename=filename, args=args)
    return await pick(results, filename, "seqtk seq")


async def seqtk_cutn(fasta: File, prefix: str = "sample", args: str = "") -> File:
    """The gaps (runs of N) in ``fasta`` as a BED, ``seqtk cutN -g``: ``<prefix>.bed``."""
    results = await seqtk_cutn_cmd(fasta=fasta, prefix=prefix, args=args)
    return await pick(results, f"{prefix}.bed", "seqtk cutN")


async def seqtk_subseq(sequences: File, filter_list: File, prefix: str = "sample", args: str = "") -> File:
    """The sequences or regions in ``filter_list`` from ``sequences``: ``<prefix>.fa.gz`` or ``<prefix>.fq.gz``."""
    ext = "fq" if SUBSEQ_FASTQ_NAME.fullmatch(file_name(sequences)) else "fa"
    filename = f"{prefix}.{ext}.gz"
    results = await seqtk_subseq_cmd(sequences=sequences, filter_list=filter_list, filename=filename, args=args)
    return await pick(results, filename, "seqtk subseq")


async def seqtk_mergepe(reads: list[File], prefix: str = "sample", args: str = "") -> File:
    """Paired ``reads`` interleaved (one file passed through): ``<prefix>.fastq.gz``."""
    if len(reads) not in (1, 2):
        raise ValueError(f"seqtk_mergepe needs one (single-end) or two (paired-end) files, got {len(reads)}")
    results = await seqtk_mergepe_cmd(
        reads_1=reads[0], reads_2=reads[1] if len(reads) == 2 else None, prefix=prefix, args=args
    )
    return await pick(results, f"{prefix}.fastq.gz", "seqtk mergepe")
