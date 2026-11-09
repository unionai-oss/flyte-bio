"""pairtools — Hi-C contact pairs: parse, sort, deduplicate, filter and annotate.

Exposes:

- :func:`pairtools_parse` — parse a name-sorted SAM/BAM of Hi-C read pairs
  into ``<prefix>.pairsam.gz``, with its statistics in ``<prefix>.pairsam.stat``.
- :func:`pairtools_sort` — sort a pairs/pairsam file into ``<prefix>.pairs.gz``,
  on the task's CPUs (``--nproc``) with half its memory as the sort buffer, as
  upstream.
- :func:`pairtools_dedup` — remove duplicate pairs, writing
  ``<prefix>.pairs.gz`` and its statistics ``<prefix>.pairs.stat``.
- :func:`pairtools_flip` — flip pairs so the first side is upstream of the
  second (by the chrom sizes' order), as ``<prefix>.flip.gz``.
- :func:`pairtools_restrict` — annotate pairs with the restriction fragments
  (a BED) each side falls in, as ``<prefix>.pairs.gz``.
- :func:`pairtools_select` — split pairs by a Python ``condition`` on their
  fields into ``<prefix>.selected.pairs.gz`` and
  ``<prefix>.unselected.pairs.gz``.

``args`` is upstream's ``ext.args`` and is split on whitespace. Upstream's
select passes its ``ext.args`` as the (quoted) condition; here the condition is
its own input, so it may contain spaces, and ``args`` carries other options.

pairtools records its command line in the output headers
(``#samheader: @PG ... CL:...``), so outputs are compared by content rather
than md5.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

PAIRTOOLS_IMAGE = "community.wave.seqera.io/library/pairtools:1.1.3--f7f48429aa47aeaf"

# Upstream's process_low / process_medium / process_high, scaled down: enough
# for test data. Genome-wide sorts want more memory (the sort buffer is half).
LOW_RESOURCES = flyte.Resources(cpu=2, memory="4Gi")
MEDIUM_RESOURCES = flyte.Resources(cpu=2, memory="8Gi")
HIGH_RESOURCES = flyte.Resources(cpu=4, memory="16Gi")

pairtools_parse_cmd = shell.create(
    name="pairtools_parse",
    image=PAIRTOOLS_IMAGE,
    resources=LOW_RESOURCES,
    inputs={"alignments": File, "chromsizes": File, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        ALN=({inputs.alignments}); SIZES=({inputs.chromsizes})
        P={inputs.prefix}
        ARGS={inputs.args}
        OUT={outputs.results}
        pairtools parse -c "${SIZES[0]}" $ARGS \
            --output-stats "$OUT/$P.pairsam.stat" -o "$OUT/$P.pairsam.gz" "${ALN[0]}" || exit 1
    """,
)

pairtools_sort_cmd = shell.create(
    name="pairtools_sort",
    image=PAIRTOOLS_IMAGE,
    resources=HIGH_RESOURCES,
    inputs={"pairs": File, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        PAIRS=({{inputs.pairs}})
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        BUFFER=$(( TASK_MEMORY_MB / 1024 / 2 ))
        if [ $BUFFER -lt 1 ]; then BUFFER=1; fi
        pairtools sort $ARGS --nproc $TASK_CPUS --memory ${{BUFFER}}G \
            -o {{outputs.results}}/"$P.pairs.gz" "${{PAIRS[0]}}" || exit 1
    """,
)

pairtools_dedup_cmd = shell.create(
    name="pairtools_dedup",
    image=PAIRTOOLS_IMAGE,
    resources=HIGH_RESOURCES,
    inputs={"pairs": File, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        PAIRS=({inputs.pairs})
        P={inputs.prefix}
        ARGS={inputs.args}
        OUT={outputs.results}
        pairtools dedup $ARGS -o "$OUT/$P.pairs.gz" --output-stats "$OUT/$P.pairs.stat" "${PAIRS[0]}" || exit 1
    """,
)

pairtools_flip_cmd = shell.create(
    name="pairtools_flip",
    image=PAIRTOOLS_IMAGE,
    resources=LOW_RESOURCES,
    inputs={"pairs": File, "chromsizes": File, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        PAIRS=({inputs.pairs}); SIZES=({inputs.chromsizes})
        ARGS={inputs.args}
        pairtools flip -c "${SIZES[0]}" $ARGS -o {outputs.results}/{inputs.prefix}.flip.gz "${PAIRS[0]}" || exit 1
    """,
)

pairtools_restrict_cmd = shell.create(
    name="pairtools_restrict",
    image=PAIRTOOLS_IMAGE,
    resources=HIGH_RESOURCES,
    inputs={"pairs": File, "fragments": File, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        PAIRS=({inputs.pairs}); FRAG=({inputs.fragments})
        ARGS={inputs.args}
        pairtools restrict -f "${FRAG[0]}" $ARGS -o {outputs.results}/{inputs.prefix}.pairs.gz "${PAIRS[0]}" \
            || exit 1
    """,
)

pairtools_select_cmd = shell.create(
    name="pairtools_select",
    image=PAIRTOOLS_IMAGE,
    resources=MEDIUM_RESOURCES,
    inputs={"pairs": File, "condition": str, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        PAIRS=({inputs.pairs})
        P={inputs.prefix}
        ARGS={inputs.args}
        OUT={outputs.results}
        set -f
        pairtools select $ARGS {inputs.condition} -o "$OUT/$P.selected.pairs.gz" \
            --output-rest "$OUT/$P.unselected.pairs.gz" "${PAIRS[0]}" || exit 1
    """,
)

env = flyte.TaskEnvironment.from_task(
    "pairtools",
    pairtools_parse_cmd.as_task(),
    pairtools_sort_cmd.as_task(),
    pairtools_dedup_cmd.as_task(),
    pairtools_flip_cmd.as_task(),
    pairtools_restrict_cmd.as_task(),
    pairtools_select_cmd.as_task(),
)


async def pick(results: Dir, name: str) -> File:
    found = await results.get_file(name)
    if found is None:
        raise FileNotFoundError(f"pairtools wrote no {name}")
    return found


@dataclass
class PairsWithStats:
    pairs: File  # <prefix>.pairsam.gz (parse) or <prefix>.pairs.gz (dedup)
    stat: File  # <prefix>.pairsam.stat (parse) or <prefix>.pairs.stat (dedup)
    results: Dir


@dataclass
class SelectResult:
    selected: File  # <prefix>.selected.pairs.gz
    unselected: File  # <prefix>.unselected.pairs.gz
    results: Dir


async def pairtools_parse(alignments: File, chromsizes: File, prefix: str, args: str = "") -> PairsWithStats:
    """Parse ``alignments`` (SAM/BAM) into ``<prefix>.pairsam.gz`` and ``<prefix>.pairsam.stat``."""
    results = await pairtools_parse_cmd(alignments=alignments, chromsizes=chromsizes, prefix=prefix, args=args)
    return PairsWithStats(
        pairs=await pick(results, f"{prefix}.pairsam.gz"),
        stat=await pick(results, f"{prefix}.pairsam.stat"),
        results=results,
    )


async def pairtools_sort(pairs: File, prefix: str, args: str = "") -> File:
    """``<prefix>.pairs.gz``: ``pairs`` sorted."""
    return await pick(await pairtools_sort_cmd(pairs=pairs, prefix=prefix, args=args), f"{prefix}.pairs.gz")


async def pairtools_dedup(pairs: File, prefix: str, args: str = "") -> PairsWithStats:
    """Deduplicate sorted ``pairs`` into ``<prefix>.pairs.gz`` and ``<prefix>.pairs.stat``."""
    results = await pairtools_dedup_cmd(pairs=pairs, prefix=prefix, args=args)
    return PairsWithStats(
        pairs=await pick(results, f"{prefix}.pairs.gz"),
        stat=await pick(results, f"{prefix}.pairs.stat"),
        results=results,
    )


async def pairtools_flip(pairs: File, chromsizes: File, prefix: str, args: str = "") -> File:
    """``<prefix>.flip.gz``: ``pairs`` flipped to upper-triangular order of ``chromsizes``."""
    results = await pairtools_flip_cmd(pairs=pairs, chromsizes=chromsizes, prefix=prefix, args=args)
    return await pick(results, f"{prefix}.flip.gz")


async def pairtools_restrict(pairs: File, fragments: File, prefix: str, args: str = "") -> File:
    """``<prefix>.pairs.gz``: ``pairs`` annotated with their restriction ``fragments``."""
    results = await pairtools_restrict_cmd(pairs=pairs, fragments=fragments, prefix=prefix, args=args)
    return await pick(results, f"{prefix}.pairs.gz")


async def pairtools_select(pairs: File, condition: str, prefix: str, args: str = "") -> SelectResult:
    """Split ``pairs`` by ``condition`` into ``<prefix>.selected.pairs.gz`` and ``<prefix>.unselected.pairs.gz``."""
    results = await pairtools_select_cmd(pairs=pairs, condition=condition, prefix=prefix, args=args)
    return SelectResult(
        selected=await pick(results, f"{prefix}.selected.pairs.gz"),
        unselected=await pick(results, f"{prefix}.unselected.pairs.gz"),
        results=results,
    )
