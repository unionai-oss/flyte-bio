"""gappa — analyse phylogenetic placements (jplace files) with gappa.

Exposes four of gappa's commands; each runs on the task's CPUs (``--threads``)
and names its outputs ``<prefix>.<gappa's name>``, as upstream (gappa's
``--file-prefix <prefix>.``). gappa reads gzipped jplace files itself.

- :func:`gappa_examine_assign` — taxonomic assignment of the placements in
  ``jplace`` from the reference sequences' taxonomy (``taxonomy``, a
  two-column TSV of name and ``;``-separated lineage; gzipped files are
  decompressed first, as gappa can't read them). Writes
  ``<prefix>.profile.tsv`` and ``<prefix>.labelled_tree.newick``, and with
  ``args`` ``--per-query-results``, ``--krona`` and ``--sativa`` also
  ``<prefix>.per_query.tsv``, ``<prefix>.krona.profile`` and
  ``<prefix>.sativa.tsv``.
- :func:`gappa_examine_graft` — the reference tree with each placement
  grafted on as a new leaf: ``<prefix>.<jplace name>.newick``.
- :func:`gappa_examine_heat_tree` — the reference tree coloured by placement
  mass. The tree formats come from ``args`` (``--write-newick-tree`` and so
  on: ``<prefix>.tree.{newick,nexus,phyloxml,svg}``). gappa's log goes to
  ``<prefix>.log``, and its colour scale (the log's ``At ...`` lines) to
  ``<prefix>.colours.txt``.
- :func:`gappa_edit_merge` — several jplace files on one reference tree
  merged into one: ``<prefix>.merge.jplace``.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

# Upstream's Docker image: gappa 0.9.0.
GAPPA_IMAGE = "quay.io/biocontainers/gappa:0.9.0--h077b44d_0"
# Upstream's process_low (graft, heat-tree, merge) and process_medium (assign),
# sized down: gappa handles the placements of a run in well under a gigabyte.
LOW_RESOURCES = flyte.Resources(cpu=2, memory="4Gi")
MEDIUM_RESOURCES = flyte.Resources(cpu=4, memory="8Gi")

gappa_examine_assign_cmd = shell.create(
    name="gappa_examine_assign",
    image=GAPPA_IMAGE,
    resources=MEDIUM_RESOURCES,
    inputs={"jplace": File, "taxonomy": File, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        JPLACE=({{inputs.jplace}})
        TAX=({{inputs.taxonomy}})
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        TAXFILE="${{TAX[0]}}"
        case "$TAXFILE" in
            *.gz)
                TAXFILE=$(mktemp -d)/$(basename "$TAXFILE" .gz)
                gzip -cd "${{TAX[0]}}" > "$TAXFILE" || exit 1
                ;;
        esac
        cd "$OUT"
        gappa \
            examine assign \
            $ARGS \
            --threads $TASK_CPUS \
            --jplace-path "${{JPLACE[0]}}" \
            --taxon-file "$TAXFILE" \
            --file-prefix "$P." || exit 1
    """,
)

gappa_examine_graft_cmd = shell.create(
    name="gappa_examine_graft",
    image=GAPPA_IMAGE,
    resources=LOW_RESOURCES,
    inputs={"jplace": File, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        JPLACE=({{inputs.jplace}})
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        cd "$OUT"
        gappa \
            examine \
            graft \
            $ARGS \
            --threads $TASK_CPUS \
            --file-prefix "$P." \
            --jplace-path "${{JPLACE[0]}}" || exit 1
    """,
)

gappa_examine_heat_tree_cmd = shell.create(
    name="gappa_examine_heat_tree",
    image=GAPPA_IMAGE,
    resources=LOW_RESOURCES,
    inputs={"jplace": File, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        JPLACE=({{inputs.jplace}})
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        cd "$OUT"
        gappa \
            examine \
            heat-tree \
            --threads $TASK_CPUS \
            --file-prefix "$P." \
            --jplace-path "${{JPLACE[0]}}" \
            $ARGS \
            --log-file "$P.log" || exit 1
        grep '^ *At' "$P.log" > "$P.colours.txt" || exit 1
    """,
)

# A list[File] is staged under the files' names and globbed alphabetically, so
# ``names`` carries the caller's order (one name per line) into the script.
gappa_edit_merge_cmd = shell.create(
    name="gappa_edit_merge",
    image=GAPPA_IMAGE,
    resources=LOW_RESOURCES,
    inputs={"jplaces": list[File], "names": str, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        IN=({{inputs.jplaces}})
        NAMES={{inputs.names}}
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        D=$(dirname "${{IN[0]}}")
        PATHS=()
        while IFS= read -r NAME; do PATHS+=("$D/$NAME"); done <<< "$NAMES"
        cd "$OUT"
        gappa \
            edit \
            merge \
            $ARGS \
            --threads $TASK_CPUS \
            --file-prefix "$P." \
            --jplace-path "${{PATHS[@]}}" || exit 1
    """,
)

env = flyte.TaskEnvironment.from_task(
    "gappa",
    gappa_examine_assign_cmd.as_task(),
    gappa_examine_graft_cmd.as_task(),
    gappa_examine_heat_tree_cmd.as_task(),
    gappa_edit_merge_cmd.as_task(),
)


def file_name(file: File) -> str:
    return file.path.rstrip("/").rsplit("/", 1)[-1]


async def pick(results: Dir, name: str, tool: str) -> File:
    found = await results.get_file(name)
    if found is None:
        raise FileNotFoundError(f"{tool} wrote no {name}")
    return found


@dataclass
class GappaAssignResult:
    profile: File  # <prefix>.profile.tsv
    labelled_tree: File  # <prefix>.labelled_tree.newick
    per_query: File | None  # <prefix>.per_query.tsv, with --per-query-results
    krona: File | None  # <prefix>.krona.profile, with --krona
    sativa: File | None  # <prefix>.sativa.tsv, with --sativa
    results: Dir  # everything gappa wrote


async def gappa_examine_assign(
    jplace: File, taxonomy: File, prefix: str = "sample", args: str = ""
) -> GappaAssignResult:
    """Taxonomic assignment of ``jplace``'s placements from the references' ``taxonomy``."""
    results = await gappa_examine_assign_cmd(jplace=jplace, taxonomy=taxonomy, prefix=prefix, args=args)
    tool = "gappa examine assign"
    return GappaAssignResult(
        profile=await pick(results, f"{prefix}.profile.tsv", tool),
        labelled_tree=await pick(results, f"{prefix}.labelled_tree.newick", tool),
        per_query=await results.get_file(f"{prefix}.per_query.tsv"),
        krona=await results.get_file(f"{prefix}.krona.profile"),
        sativa=await results.get_file(f"{prefix}.sativa.tsv"),
        results=results,
    )


async def gappa_examine_graft(jplace: File, prefix: str = "sample", args: str = "") -> File:
    """The reference tree with ``jplace``'s placements grafted on: ``<prefix>.<jplace name>.newick``."""
    results = await gappa_examine_graft_cmd(jplace=jplace, prefix=prefix, args=args)
    trees = [f for f in await results.list_files() if f.name.endswith(".newick")]
    if len(trees) != 1:
        raise FileNotFoundError(f"gappa examine graft wrote {len(trees)} newick trees, expected 1")
    return trees[0]


@dataclass
class GappaHeatTreeResult:
    newick: File | None  # <prefix>.tree.newick, with --write-newick-tree
    nexus: File | None  # <prefix>.tree.nexus, with --write-nexus-tree
    phyloxml: File | None  # <prefix>.tree.phyloxml, with --write-phyloxml-tree
    svg: File | None  # <prefix>.tree.svg, with --write-svg-tree
    colours: File  # <prefix>.colours.txt, the colour scale
    log: File  # <prefix>.log
    results: Dir  # everything gappa wrote


async def gappa_examine_heat_tree(jplace: File, prefix: str = "sample", args: str = "") -> GappaHeatTreeResult:
    """The reference tree coloured by ``jplace``'s placement mass, in the formats ``args`` asks for."""
    results = await gappa_examine_heat_tree_cmd(jplace=jplace, prefix=prefix, args=args)
    tool = "gappa examine heat-tree"
    return GappaHeatTreeResult(
        newick=await results.get_file(f"{prefix}.tree.newick"),
        nexus=await results.get_file(f"{prefix}.tree.nexus"),
        phyloxml=await results.get_file(f"{prefix}.tree.phyloxml"),
        svg=await results.get_file(f"{prefix}.tree.svg"),
        colours=await pick(results, f"{prefix}.colours.txt", tool),
        log=await pick(results, f"{prefix}.log", tool),
        results=results,
    )


async def gappa_edit_merge(jplaces: list[File], prefix: str = "sample", args: str = "") -> File:
    """``jplaces`` (on one reference tree) merged into one: ``<prefix>.merge.jplace``."""
    if not jplaces:
        raise ValueError("gappa_edit_merge needs at least one jplace file")
    names = [file_name(f) for f in jplaces]
    if len(set(names)) != len(names):
        raise ValueError(f"gappa_edit_merge inputs need distinct names, got {names}")
    results = await gappa_edit_merge_cmd(jplaces=jplaces, names="\n".join(names), prefix=prefix, args=args)
    merged = [f for f in await results.list_files() if ".merge.jplace" in f.name]
    if len(merged) != 1:
        raise FileNotFoundError(f"gappa edit merge wrote {len(merged)} merged jplace files, expected 1")
    return merged[0]
