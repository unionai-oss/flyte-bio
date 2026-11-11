"""tidk — telomere identification: find, count and plot telomeric repeats.

Exposes:

- :func:`tidk_explore` — the repeats of a FASTA that look telomeric
  (``tidk explore``; repeat lengths as ``args``, e.g. ``--minimum 5 --maximum
  30``): ``<prefix>.tidk.explore.tsv`` and, as upstream, when it found any,
  the most frequent repeat in ``<prefix>.top.sequence.txt``.
- :func:`tidk_search` — counts of a given telomeric repeat (``TTAGGG``) in
  windows along each sequence (``tidk search``): ``<prefix>.tsv``, or
  ``<prefix>.bedgraph`` with ``--extension bedgraph`` in ``args``.
- :func:`tidk_plot` — a plot of :func:`tidk_search`'s TSV: ``<prefix>.svg``.

tidk is a Rust tool: a panic exits with 101, and other failures may exit 255,
which the platform retries forever, so every failure exits 1.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

TIDK_IMAGE = "community.wave.seqera.io/library/tidk:0.2.7--e781b1e8244cc913"
# Upstream's process_single.
DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory="6Gi")

tidk_explore_cmd = shell.create(
    name="tidk_explore",
    image=TIDK_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"fasta": File, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        FA=({inputs.fasta})
        P={inputs.prefix}
        ARGS={inputs.args}
        OUT={outputs.results}
        W=$(mktemp -d); cd "$W"
        ln -s "${FA[0]}" "$(basename "${FA[0]}")"
        tidk explore $ARGS "$(basename "${FA[0]}")" > "$OUT/$P.tidk.explore.tsv" || exit 1
        # As upstream: the top repeat is the first column of the table's first row.
        if [ "$(wc -l < "$OUT/$P.tidk.explore.tsv")" -gt 1 ]; then
            sed -n 2p "$OUT/$P.tidk.explore.tsv" | awk '{print $1;}' > "$OUT/$P.top.sequence.txt"
        else
            echo "No sequence identified"
        fi
    """,
)

tidk_search_cmd = shell.create(
    name="tidk_search",
    image=TIDK_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"fasta": File, "string": str, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        FA=({inputs.fasta})
        STRING={inputs.string}
        P={inputs.prefix}
        ARGS={inputs.args}
        OUT={outputs.results}
        W=$(mktemp -d); cd "$W"
        ln -s "${FA[0]}" "$(basename "${FA[0]}")"
        tidk search --string "$STRING" --output "$P" --dir tidk $ARGS "$(basename "${FA[0]}")" || exit 1
        for EXT in tsv bedgraph; do
            if [ -f "tidk/${P}_telomeric_repeat_windows.$EXT" ]; then
                mv "tidk/${P}_telomeric_repeat_windows.$EXT" "$OUT/$P.$EXT"
            else
                echo "$EXT file was not produced"
            fi
        done
    """,
)

tidk_plot_cmd = shell.create(
    name="tidk_plot",
    image=TIDK_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"tsv": File, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        TSV=({inputs.tsv})
        P={inputs.prefix}
        ARGS={inputs.args}
        W=$(mktemp -d); cd "$W"
        ln -s "${TSV[0]}" "$(basename "${TSV[0]}")"
        tidk plot --output "$P" $ARGS --tsv "$(basename "${TSV[0]}")" || exit 1
        mv "$P.svg" {outputs.results}/
    """,
)

env = flyte.TaskEnvironment.from_task(
    "tidk", tidk_explore_cmd.as_task(), tidk_search_cmd.as_task(), tidk_plot_cmd.as_task()
)


@dataclass
class TidkExploreResult:
    explore_tsv: File  # <prefix>.tidk.explore.tsv
    top_sequence: File | None  # <prefix>.top.sequence.txt, when a repeat was found


@dataclass
class TidkSearchResult:
    tsv: File | None  # <prefix>.tsv (the default)
    bedgraph: File | None  # <prefix>.bedgraph (with --extension bedgraph)


async def tidk_explore(fasta: File, prefix: str = "sample", args: str = "") -> TidkExploreResult:
    """Candidate telomeric repeats of ``fasta``, and the most frequent one when there are any."""
    results = await tidk_explore_cmd(fasta=fasta, prefix=prefix, args=args)
    tsv = await results.get_file(f"{prefix}.tidk.explore.tsv")
    if tsv is None:
        raise FileNotFoundError(f"tidk explore wrote no {prefix}.tidk.explore.tsv")
    return TidkExploreResult(explore_tsv=tsv, top_sequence=await results.get_file(f"{prefix}.top.sequence.txt"))


async def tidk_search(fasta: File, string: str, prefix: str = "sample", args: str = "") -> TidkSearchResult:
    """Windowed counts of the telomeric repeat ``string`` along each sequence of ``fasta``."""
    results = await tidk_search_cmd(fasta=fasta, string=string, prefix=prefix, args=args)
    return TidkSearchResult(
        tsv=await results.get_file(f"{prefix}.tsv"), bedgraph=await results.get_file(f"{prefix}.bedgraph")
    )


async def tidk_plot(tsv: File, prefix: str = "sample", args: str = "") -> File:
    """A plot of :func:`tidk_search`'s windowed counts: ``<prefix>.svg``."""
    results = await tidk_plot_cmd(tsv=tsv, prefix=prefix, args=args)
    found = await results.get_file(f"{prefix}.svg")
    if found is None:
        raise FileNotFoundError(f"tidk plot wrote no {prefix}.svg")
    return found
