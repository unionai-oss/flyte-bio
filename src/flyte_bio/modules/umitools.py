"""umitools — UMI handling with UMI-tools.

Exposes:

- :func:`umitools_extract` — move the UMI from the read sequence into the
  read name (one or two mates).
- :func:`umitools_dedup` — deduplicate a coordinate-sorted, indexed BAM by
  UMI and mapping position, optionally with UMI statistics.
- :func:`umitools_prepareforrsem` — make a name-sorted, deduplicated
  paired-end BAM consistent again for RSEM/Salmon (dedup can leave mates
  without their partner and the flags out of step).

UMI patterns are regexes (``--bc-pattern='^(?P<umi_1>CGA.{8}).*'``) passed
in ``args``, which the scripts expand unquoted so it can carry several flags;
globbing is switched off before that, or bash would try to match the
pattern's ``?`` and ``*`` against file names.

Each tool writes into one output Dir under the names upstream uses
(``<prefix>.umi_extract[_1|_2].fastq.gz``, ``<prefix>.bam``,
``<prefix>.log`` …); the wrappers pick the files out by name.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

UMITOOLS_IMAGE = "community.wave.seqera.io/library/umi_tools_future_matplotlib_numpy_pruned:1ee668bafc8c9f81"

# UMI-tools is single-threaded; dedup and prepare-for-rsem hold the reads of a
# position (or a read name) in memory, so they get the memory.
EXTRACT_RESOURCES = flyte.Resources(cpu=1, memory="6Gi")
DEDUP_RESOURCES = flyte.Resources(cpu=1, memory="16Gi")

umitools_extract_cmd = shell.create(
    name="umitools_extract",
    image=UMITOOLS_IMAGE,
    resources=EXTRACT_RESOURCES,
    inputs={"reads_1": File, "reads_2": File | None, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        R1=({inputs.reads_1})
        shopt -s nullglob; R2=({inputs.reads_2}); shopt -u nullglob
        set -f
        P={inputs.prefix}
        ARGS={inputs.args}
        OUT={outputs.results}
        if [ ${#R2[@]} -eq 0 ]; then
            umi_tools extract -I "${R1[0]}" -S "$OUT/$P.umi_extract.fastq.gz" \
                $ARGS > "$OUT/$P.umi_extract.log"
        else
            umi_tools extract -I "${R1[0]}" --read2-in="${R2[0]}" \
                -S "$OUT/${P}.umi_extract_1.fastq.gz" --read2-out="$OUT/${P}.umi_extract_2.fastq.gz" \
                $ARGS > "$OUT/$P.umi_extract.log"
        fi
    """,
)

# As upstream: a fixed --random-seed unless `args` sets one, and matplotlib
# kept out of /tmp. umi_tools reads the index next to the BAM, so both are
# linked side by side first.
umitools_dedup_cmd = shell.create(
    name="umitools_dedup",
    image=UMITOOLS_IMAGE,
    resources=DEDUP_RESOURCES,
    inputs={"bam": File, "bai": File, "prefix": str, "paired": bool, "output_stats": bool, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        BAM=({inputs.bam}); BAI=({inputs.bai})
        set -f
        P={inputs.prefix}
        ARGS={inputs.args}
        OUT={outputs.results}
        case " $ARGS " in *--random-seed*) ;; *) ARGS="$ARGS --random-seed=100" ;; esac
        PAIRED=(); [ {inputs.paired} = true ] && PAIRED=(--paired)
        STATS=(); [ {inputs.output_stats} = true ] && STATS=(--output-stats "$OUT/$P")
        W=$(mktemp -d); cd "$W"
        ln -s "${BAM[0]}" in.bam
        ln -s "${BAI[0]}" in.bam.bai
        mkdir .tmp && chmod 777 .tmp
        MPLCONFIGDIR=.tmp TMPDIR=.tmp PYTHONHASHSEED=0 umi_tools dedup \
            -I in.bam -S "$OUT/$P.bam" -L "$OUT/$P.log" \
            "${STATS[@]}" "${PAIRED[@]}" $ARGS
    """,
)

umitools_prepareforrsem_cmd = shell.create(
    name="umitools_prepareforrsem",
    image=UMITOOLS_IMAGE,
    resources=DEDUP_RESOURCES,
    inputs={"bam": File, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        BAM=({inputs.bam})
        set -f
        P={inputs.prefix}
        ARGS={inputs.args}
        OUT={outputs.results}
        umi_tools prepare-for-rsem --stdin="${BAM[0]}" --stdout="$OUT/$P.bam" \
            --log="$OUT/$P.prepare_for_rsem.log" $ARGS
    """,
)


env = flyte.TaskEnvironment.from_task(
    "umitools",
    umitools_extract_cmd.as_task(),
    umitools_dedup_cmd.as_task(),
    umitools_prepareforrsem_cmd.as_task(),
)


async def pick(results: Dir, name: str, tool: str) -> File:
    found = await results.get_file(name)
    if found is None:
        raise FileNotFoundError(f"umi_tools {tool} wrote no {name}")
    return found


@dataclass
class UmiExtractResult:
    reads_1: File
    reads_2: File | None  # None for single-end
    log: File


async def umitools_extract(
    reads_1: File,
    reads_2: File | None = None,
    prefix: str = "umi",
    args: str = "",
) -> UmiExtractResult:
    """Extract UMIs into the read names; ``reads_2`` is None for single-end."""
    results = await umitools_extract_cmd(reads_1=reads_1, reads_2=reads_2, prefix=prefix, args=args)
    log = await pick(results, f"{prefix}.umi_extract.log", "extract")
    if reads_2 is None:
        return UmiExtractResult(
            reads_1=await pick(results, f"{prefix}.umi_extract.fastq.gz", "extract"), reads_2=None, log=log
        )
    return UmiExtractResult(
        reads_1=await pick(results, f"{prefix}.umi_extract_1.fastq.gz", "extract"),
        reads_2=await pick(results, f"{prefix}.umi_extract_2.fastq.gz", "extract"),
        log=log,
    )


@dataclass
class UmiDedupResult:
    bam: File
    log: File
    # With output_stats only:
    edit_distance: File | None = None
    per_umi: File | None = None
    per_umi_per_position: File | None = None


async def umitools_dedup(
    bam: File,
    bai: File,
    prefix: str,
    paired: bool,
    output_stats: bool = False,
    args: str = "",
) -> UmiDedupResult:
    """Deduplicate ``bam`` (coordinate-sorted, with index ``bai``) by UMI."""
    results = await umitools_dedup_cmd(
        bam=bam, bai=bai, prefix=prefix, paired=paired, output_stats=output_stats, args=args
    )
    result = UmiDedupResult(
        bam=await pick(results, f"{prefix}.bam", "dedup"),
        log=await pick(results, f"{prefix}.log", "dedup"),
    )
    if output_stats:
        result.edit_distance = await pick(results, f"{prefix}_edit_distance.tsv", "dedup")
        result.per_umi = await pick(results, f"{prefix}_per_umi.tsv", "dedup")
        result.per_umi_per_position = await pick(results, f"{prefix}_per_umi_per_position.tsv", "dedup")
    return result


@dataclass
class PrepareForRsemResult:
    bam: File
    log: File


async def umitools_prepareforrsem(bam: File, prefix: str, args: str = "") -> PrepareForRsemResult:
    """Fix up a name-sorted, UMI-deduplicated paired-end BAM for RSEM/Salmon."""
    results = await umitools_prepareforrsem_cmd(bam=bam, prefix=prefix, args=args)
    return PrepareForRsemResult(
        bam=await pick(results, f"{prefix}.bam", "prepare-for-rsem"),
        log=await pick(results, f"{prefix}.prepare_for_rsem.log", "prepare-for-rsem"),
    )
