"""epang — phylogenetic placement of query sequences on a reference tree with EPA-ng.

Exposes :func:`epang_place`, which places the aligned ``query`` sequences on
``reference_tree`` given its alignment, ``reference_alignment`` (``epa-ng
--query --ref-msa --tree``), or on a reference from a binary file written by
an earlier ``--dump-binary`` run (``binary``; it can't be combined with a
reference alignment or tree). ``bfast`` passes a bfast file (``--bfast``).
The model and other options (``--model LG``) go in ``args``. As upstream:

- epa-ng reads gzipped alignments itself, but not a gzipped tree, so a tree
  named ``*.gz`` is decompressed first;
- epa-ng's own outputs are kept under their names (``epa_result.jplace.gz``,
  gzipped; ``epa_info.log``; ``epa_binary_file`` with ``--dump-binary``),
  with the placements and log also copied to ``<prefix>.epa_result.jplace.gz``
  and ``<prefix>.epa_info.log``. A run that places nothing (such as
  ``--dump-binary``) writes no jplace.

Placement runs on the task's CPUs (``--threads``).

:func:`epang_split` splits an alignment of reference and query sequences,
``full_alignment``, into the two, given the reference alignment
(``epa-ng --split``): ``<prefix>.query.fasta.gz`` and
``<prefix>.reference.fasta.gz``.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

# Upstream's Docker image: EPA-ng 0.3.8.
EPANG_IMAGE = "quay.io/biocontainers/epa-ng:0.3.8--h9a82719_1"
# Upstream labels placement process_high; the test data places in seconds.
PLACE_RESOURCES = flyte.Resources(cpu=4, memory="8Gi")
# Upstream's process_single.
SPLIT_RESOURCES = flyte.Resources(cpu=1, memory="2Gi")

epang_place_cmd = shell.create(
    name="epang_place",
    image=EPANG_IMAGE,
    resources=PLACE_RESOURCES,
    inputs={
        "query": File | None,
        "reference_alignment": File | None,
        "reference_tree": File | None,
        "bfast": File | None,
        "binary": File | None,
        "prefix": str,
        "args": str,
    },
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        shopt -s nullglob
        QUERY=({{inputs.query}})
        REFALN=({{inputs.reference_alignment}})
        TREE=({{inputs.reference_tree}})
        BFAST=({{inputs.bfast}})
        BINARY=({{inputs.binary}})
        shopt -u nullglob
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"

        OPTS=()
        if [ ${{#QUERY[@]}} -gt 0 ]; then OPTS+=(--query "${{QUERY[0]}}"); fi
        if [ ${{#REFALN[@]}} -gt 0 ]; then OPTS+=(--ref-msa "${{REFALN[0]}}"); fi
        TREEFILE=""
        if [ ${{#TREE[@]}} -gt 0 ]; then
            # epa-ng can't parse a gzipped tree ("Treeparsing failed!").
            case "${{TREE[0]}}" in
                *.gz)
                    TREEFILE=$(mktemp -d)/$(basename "${{TREE[0]}}" .gz)
                    gzip -cd "${{TREE[0]}}" > "$TREEFILE" || exit 1
                    OPTS+=(--tree "$TREEFILE")
                    ;;
                *) OPTS+=(--tree "${{TREE[0]}}") ;;
            esac
        fi
        if [ ${{#BFAST[@]}} -gt 0 ]; then OPTS+=(--bfast "${{BFAST[0]}}"); fi
        if [ ${{#BINARY[@]}} -gt 0 ]; then OPTS+=(--binary "${{BINARY[0]}}"); fi

        epa-ng \
            $ARGS \
            --threads $TASK_CPUS \
            "${{OPTS[@]}}" || exit 1

        if [ -e epa_result.jplace ]; then
            gzip epa_result.jplace
            cp epa_result.jplace.gz "$P.epa_result.jplace.gz"
        fi
        if [ -e epa_info.log ]; then cp epa_info.log "$P.epa_info.log"; fi
        mv ./* "$OUT"/
    """,
)

epang_split_cmd = shell.create(
    name="epang_split",
    image=EPANG_IMAGE,
    resources=SPLIT_RESOURCES,
    inputs={"reference_alignment": File, "full_alignment": File, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        REFALN=({inputs.reference_alignment})
        FULL=({inputs.full_alignment})
        P={inputs.prefix}
        ARGS={inputs.args}
        OUT={outputs.results}
        W=$(mktemp -d); cd "$W"
        epa-ng $ARGS --split "${REFALN[0]}" "${FULL[0]}" || exit 1
        gzip -c query.fasta > "$OUT/$P.query.fasta.gz"
        gzip -c reference.fasta > "$OUT/$P.reference.fasta.gz"
    """,
)

env = flyte.TaskEnvironment.from_task("epang", epang_place_cmd.as_task(), epang_split_cmd.as_task())


async def pick(results: Dir, name: str, tool: str) -> File:
    found = await results.get_file(name)
    if found is None:
        raise FileNotFoundError(f"{tool} wrote no {name}")
    return found


@dataclass
class EpangPlaceResult:
    jplace: File | None  # <prefix>.epa_result.jplace.gz; None when the run places nothing
    log: File  # <prefix>.epa_info.log
    results: Dir  # everything epa-ng wrote, plus the prefixed copies


async def epang_place(
    query: File | None = None,
    reference_alignment: File | None = None,
    reference_tree: File | None = None,
    bfast: File | None = None,
    binary: File | None = None,
    prefix: str = "sample",
    args: str = "",
) -> EpangPlaceResult:
    """Place ``query`` on ``reference_tree`` (with ``reference_alignment``), or on the reference in ``binary``."""
    if binary is not None and (reference_alignment is not None or reference_tree is not None):
        raise ValueError("[EPANG] Cannot supply both binary and reference MSA or reference tree. Check input")
    results = await epang_place_cmd(
        query=query,
        reference_alignment=reference_alignment,
        reference_tree=reference_tree,
        bfast=bfast,
        binary=binary,
        prefix=prefix,
        args=args,
    )
    return EpangPlaceResult(
        jplace=await results.get_file(f"{prefix}.epa_result.jplace.gz"),
        log=await pick(results, f"{prefix}.epa_info.log", "epa-ng"),
        results=results,
    )


@dataclass
class EpangSplitResult:
    query: File  # <prefix>.query.fasta.gz
    reference: File  # <prefix>.reference.fasta.gz


async def epang_split(
    reference_alignment: File, full_alignment: File, prefix: str = "sample", args: str = ""
) -> EpangSplitResult:
    """``full_alignment`` split into its query and reference sequences, given ``reference_alignment``."""
    results = await epang_split_cmd(
        reference_alignment=reference_alignment, full_alignment=full_alignment, prefix=prefix, args=args
    )
    return EpangSplitResult(
        query=await pick(results, f"{prefix}.query.fasta.gz", "epa-ng --split"),
        reference=await pick(results, f"{prefix}.reference.fasta.gz", "epa-ng --split"),
    )
