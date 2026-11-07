"""chromap — fast alignment and preprocessing of chromatin profiles (ATAC-seq, ChIP-seq, Hi-C, scATAC-seq).

Exposes:

- :func:`chromap_index` — Chromap's index of a FASTA, ``<prefix>.index``
  (the prefix defaults to the FASTA's name without its last extension, as
  upstream).
- :func:`chromap_chromap` — map one sample's reads (single-end lanes, or one
  pair of mates) and write the deduplicated fragments.

Outputs of :func:`chromap_chromap` keep upstream's names, by the format
``args`` selects:

- the default BED of fragments, ``<prefix>.bed.gz``;
- with ``--TagAlign``: ``<prefix>.tagAlign.gz``;
- with ``--pairs`` (Hi-C): ``<prefix>.pairs.gz``;
- with ``--SAM``: ``<prefix>.bam``, converted by ``samtools view -bh`` with
  ``args2``.

As upstream, text outputs are compressed with ``gzip -n``, and Chromap and
samtools both use the task's CPUs. Optional ``barcodes`` (single-cell reads,
one file per read file, joined with commas as upstream), ``whitelist`` (with
barcodes only), ``chr_order`` and ``pairs_chr_order`` add ``-b``,
``--barcode-whitelist``, ``--chr-order`` and ``--pairs-natural-chr-order``.
Read and barcode files are passed to Chromap in the caller's order.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

CHROMAP_IMAGE = "community.wave.seqera.io/library/chromap:0.3.2--4ec4bca51cd82195"
CHROMAP_SAMTOOLS_IMAGE = "community.wave.seqera.io/library/chromap_samtools:b975c17adf0096ba"
# Upstream's process_medium has more; a mammalian genome's index and mapping need about 20Gi.
DEFAULT_RESOURCES = flyte.Resources(cpu=6, memory="16Gi")

chromap_index_cmd = shell.create(
    name="chromap_index",
    image=CHROMAP_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"fasta": File, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        FA=({{inputs.fasta}})
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        chromap -i $ARGS -t $TASK_CPUS -r "${{FA[0]}}" -o {{outputs.results}}/"$P.index"
    """,
)

# A list[File] is staged under the files' names and globbed alphabetically, so
# the names ("/"-separated) carry the caller's order into the script.
chromap_chromap_cmd = shell.create(
    name="chromap_chromap",
    image=CHROMAP_SAMTOOLS_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "reads": list[File],
        "read_names": str,
        "single_end": bool,
        "fasta": File,
        "index": File,
        "barcodes": list[File],
        "barcode_names": str,
        "whitelist": File | None,
        "chr_order": File | None,
        "pairs_chr_order": File | None,
        "prefix": str,
        "args": str,
        "args2": str,
    },
    defaults={"args": "", "args2": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        READS=({{inputs.reads}})
        FA=({{inputs.fasta}})
        INDEX=({{inputs.index}})
        shopt -s nullglob
        BARCODES=({{inputs.barcodes}}); WHITELIST=({{inputs.whitelist}})
        CHR_ORDER=({{inputs.chr_order}}); PAIRS_CHR_ORDER=({{inputs.pairs_chr_order}})
        shopt -u nullglob
        READ_NAMES={{inputs.read_names}}
        BARCODE_NAMES={{inputs.barcode_names}}
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        ARGS2={{inputs.args2}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"

        in_order() {{
            local DIR=$1 NAMES=$2 NAME
            local -a LIST
            IFS=/ read -ra LIST <<< "$NAMES"
            ORDERED=()
            for NAME in "${{LIST[@]}}"; do ORDERED+=("$DIR/$NAME"); done
        }}
        in_order "$(dirname "${{READS[0]}}")" "$READ_NAMES"; READS=("${{ORDERED[@]}}")

        EXT=bed
        if [[ $ARGS == *--SAM* ]]; then EXT=sam
        elif [[ $ARGS == *--TagAlign* ]]; then EXT=tagAlign
        elif [[ $ARGS == *--pairs* ]]; then EXT=pairs
        fi
        EXTRA=()
        if [ ${{#BARCODES[@]}} -gt 0 ]; then
            in_order "$(dirname "${{BARCODES[0]}}")" "$BARCODE_NAMES"
            EXTRA+=(-b "$(IFS=,; echo "${{ORDERED[*]}}")")
            if [ ${{#WHITELIST[@]}} -gt 0 ]; then EXTRA+=(--barcode-whitelist "${{WHITELIST[0]}}"); fi
        fi
        if [ ${{#CHR_ORDER[@]}} -gt 0 ]; then EXTRA+=(--chr-order "${{CHR_ORDER[0]}}"); fi
        if [ ${{#PAIRS_CHR_ORDER[@]}} -gt 0 ]; then EXTRA+=(--pairs-natural-chr-order "${{PAIRS_CHR_ORDER[0]}}"); fi

        if [ {{inputs.single_end}} = true ]; then
            MATES=(-1 "$(IFS=,; echo "${{READS[*]}}")")
        else
            MATES=(-1 "${{READS[0]}}" -2 "${{READS[1]}}")
        fi
        chromap $ARGS "${{EXTRA[@]}}" -t $TASK_CPUS -x "${{INDEX[0]}}" -r "${{FA[0]}}" "${{MATES[@]}}" -o "$P.$EXT"

        if [ $EXT = sam ]; then
            samtools view $ARGS2 -@ $TASK_CPUS -bh -o "$OUT/$P.bam" "$P.sam"
        else
            gzip -n "$P.$EXT"
            mv "$P.$EXT.gz" "$OUT"/
        fi
    """,
)

chromap_index_env = flyte.TaskEnvironment.from_task("chromap_index", chromap_index_cmd.as_task())
chromap_samtools_env = flyte.TaskEnvironment.from_task("chromap_samtools", chromap_chromap_cmd.as_task())
env = flyte.TaskEnvironment(name="chromap", depends_on=[chromap_index_env, chromap_samtools_env])


def file_name(file: File) -> str:
    return file.path.rstrip("/").rsplit("/", 1)[-1]


def ordered_names(files: list[File], label: str) -> str:
    names = [file_name(f) for f in files]
    if len(set(names)) != len(names):
        raise ValueError(f"chromap {label} need distinct names, got {names}")
    return "/".join(names)


async def chromap_index(fasta: File, prefix: str = "", args: str = "") -> File:
    """Chromap index of ``fasta``: ``<prefix>.index``, by default named after the FASTA."""
    prefix = prefix or file_name(fasta).rsplit(".", 1)[0]
    results = await chromap_index_cmd(fasta=fasta, prefix=prefix, args=args)
    index = await results.get_file(f"{prefix}.index")
    if index is None:
        raise FileNotFoundError(f"chromap wrote no {prefix}.index")
    return index


def output_extension(args: str) -> str:
    """The output's extension, by the format ``args`` selects, as upstream."""
    if "--SAM" in args:
        return "bam"
    if "--TagAlign" in args:
        return "tagAlign.gz"
    if "--pairs" in args:
        return "pairs.gz"
    return "bed.gz"


@dataclass
class ChromapResult:
    output: File  # <prefix>.bed.gz, .tagAlign.gz, .pairs.gz or .bam, by args
    results: Dir  # everything the task wrote


async def chromap_chromap(
    reads: list[File],
    fasta: File,
    index: File,
    prefix: str = "sample",
    single_end: bool = False,
    barcodes: list[File] | None = None,
    whitelist: File | None = None,
    chr_order: File | None = None,
    pairs_chr_order: File | None = None,
    args: str = "",
    args2: str = "",
) -> ChromapResult:
    """Map one sample's ``reads`` against ``index``; ``args`` go to Chromap, ``args2`` to samtools.

    ``reads`` are single-end read files (several lanes are joined) or, with
    ``single_end=False``, mates 1 and 2.
    """
    if not reads:
        raise ValueError("chromap_chromap needs reads")
    if not single_end and len(reads) != 2:
        raise ValueError(f"chromap_chromap needs two mates for paired-end reads, got {len(reads)} files")
    barcodes = barcodes or []
    results = await chromap_chromap_cmd(
        reads=reads,
        read_names=ordered_names(reads, "read files"),
        single_end=single_end,
        fasta=fasta,
        index=index,
        barcodes=barcodes,
        barcode_names=ordered_names(barcodes, "barcode files"),
        whitelist=whitelist,
        chr_order=chr_order,
        pairs_chr_order=pairs_chr_order,
        prefix=prefix,
        args=args,
        args2=args2,
    )
    name = f"{prefix}.{output_extension(args)}"
    output = await results.get_file(name)
    if output is None:
        raise FileNotFoundError(f"chromap wrote no {name}")
    return ChromapResult(output=output, results=results)
