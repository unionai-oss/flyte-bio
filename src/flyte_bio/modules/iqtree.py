"""iqtree — maximum-likelihood phylogenies with IQ-TREE.

Exposes :func:`iqtree`, which runs IQ-TREE on one or more alignments
(``-s``, comma-joined) and/or input trees, with each of upstream's optional
inputs mapped to its option:

- ``tree`` → ``-t`` (starting tree, or trees to evaluate),
- ``tree_te`` → ``-te`` (fixed tree topology),
- ``lmclust`` → ``-lmclust`` (likelihood-mapping clusters),
- ``mdef`` → ``-mdef`` (NEXUS model definitions),
- ``partitions_equal`` → ``-q``, ``partitions_proportional`` → ``-spp``,
  ``partitions_unlinked`` → ``-sp`` (partition files),
- ``guide_tree`` → ``-ft``, ``sitefreq_in`` → ``-fs``,
  ``constraint_tree`` → ``-g``, ``trees_z`` → ``-z``,
  ``suptree`` → ``-sup``, ``trees_rf`` → ``-rf``.

The model and everything else go in ``args`` (e.g. ``-m LG -bb 1000``). As
upstream, IQ-TREE picks its thread count (``-nt AUTO``) up to the task's
CPUs (``-ntmax``), and gets the task's memory (``-mem``) unless a partition
model is in use, which IQ-TREE doesn't allow with ``-mem``. Alignments are
passed in name order.

Outputs keep IQ-TREE's names, ``<prefix>.<suffix>`` (``treefile``,
``iqtree``, ``mldist``, ``ufboot``, ``contree``, ``splits.nex``, ``log``, …),
plus ``GTRPMIX.nex`` for a profile-mixture exchangeability matrix; which ones
exist depends on the options.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

IQTREE_IMAGE = "quay.io/biocontainers/iqtree:3.1.3--h8471819_0"
# Upstream labels IQ-TREE process_medium; its tests run on 4 CPUs.
DEFAULT_RESOURCES = flyte.Resources(cpu=4, memory="8Gi")

# Each optional input file, with the IQ-TREE option it's passed to.
OPTIONAL_INPUTS = {
    "tree": "-t",
    "tree_te": "-te",
    "lmclust": "-lmclust",
    "mdef": "-mdef",
    "partitions_equal": "-q",
    "partitions_proportional": "-spp",
    "partitions_unlinked": "-sp",
    "guide_tree": "-ft",
    "sitefreq_in": "-fs",
    "constraint_tree": "-g",
    "trees_z": "-z",
    "suptree": "-sup",
    "trees_rf": "-rf",
}

OPTIONAL_SCRIPT = "\n".join(
    f"""        shopt -s nullglob; F=({{inputs.{name}}}); shopt -u nullglob
        if [ ${{#F[@]}} -gt 0 ]; then OPTS+=({option} "${{F[0]}}"); fi"""
    for name, option in OPTIONAL_INPUTS.items()
)

iqtree_cmd = shell.create(
    name="iqtree",
    image=IQTREE_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "alignment": list[File] | None,
        **dict.fromkeys(OPTIONAL_INPUTS, File | None),
        "prefix": str,
        "args": str,
    },
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        shopt -s nullglob; ALIGNMENT=({{inputs.alignment}}); shopt -u nullglob
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        OPTS=()
        if [ ${{#ALIGNMENT[@]}} -gt 0 ]; then
            JOINED=$(IFS=,; echo "${{ALIGNMENT[*]}}")
            OPTS+=(-s "$JOINED")
        fi
{OPTIONAL_SCRIPT}
        PARTITIONED=false
        for O in "${{OPTS[@]}}"; do
            case "$O" in -q|-spp|-sp) PARTITIONED=true ;; esac
        done
        MEMORY=(-mem "${{TASK_MEMORY_MB}}MB")
        if [ $PARTITIONED = true ]; then MEMORY=(); fi
        W=$(mktemp -d); cd "$W"
        iqtree \
            $ARGS \
            "${{OPTS[@]}}" \
            -pre "$P" \
            -nt AUTO \
            -ntmax $TASK_CPUS \
            "${{MEMORY[@]}}" || exit 1
        mv ./* "$OUT"/
    """,
)

env = flyte.TaskEnvironment.from_task("iqtree", iqtree_cmd.as_task())


@dataclass
class IqtreeResult:
    phylogeny: File | None  # <prefix>.treefile
    report: File | None  # <prefix>.iqtree
    mldist: File | None  # <prefix>.mldist
    lmap_svg: File | None  # <prefix>.lmap.svg
    lmap_eps: File | None  # <prefix>.lmap.eps
    lmap_quartetlh: File | None  # <prefix>.lmap.quartetlh
    sitefreq_out: File | None  # <prefix>.sitefreq
    bootstrap: File | None  # <prefix>.ufboot
    state: File | None  # <prefix>.state
    contree: File | None  # <prefix>.contree
    nex: list[File]  # *.nex: <prefix>.splits.nex, GTRPMIX.nex, ...
    best_scheme: File | None  # <prefix>.best_scheme
    splits: File | None  # <prefix>.splits
    suptree: File | None  # <prefix>.suptree
    alninfo: File | None  # <prefix>.alninfo
    partlh: File | None  # <prefix>.partlh
    siteprob: File | None  # <prefix>.siteprob
    sitelh: File | None  # <prefix>.sitelh
    treels: File | None  # <prefix>.treels
    rate: File | None  # <prefix>.rate
    mlrate: File | None  # <prefix>.mlrate
    exch_matrix: File | None  # GTRPMIX.nex
    log: File  # <prefix>.log
    results: Dir  # everything IQ-TREE wrote


async def iqtree(
    alignment: File | list[File] | None = None,
    tree: File | None = None,
    tree_te: File | None = None,
    lmclust: File | None = None,
    mdef: File | None = None,
    partitions_equal: File | None = None,
    partitions_proportional: File | None = None,
    partitions_unlinked: File | None = None,
    guide_tree: File | None = None,
    sitefreq_in: File | None = None,
    constraint_tree: File | None = None,
    trees_z: File | None = None,
    suptree: File | None = None,
    trees_rf: File | None = None,
    prefix: str = "sample",
    args: str = "",
) -> IqtreeResult:
    """Run IQ-TREE on ``alignment`` (one or more files) and/or the given trees, with ``args``."""
    alignments = [alignment] if isinstance(alignment, File) else list(alignment or [])
    results = await iqtree_cmd(
        alignment=alignments or None,
        tree=tree,
        tree_te=tree_te,
        lmclust=lmclust,
        mdef=mdef,
        partitions_equal=partitions_equal,
        partitions_proportional=partitions_proportional,
        partitions_unlinked=partitions_unlinked,
        guide_tree=guide_tree,
        sitefreq_in=sitefreq_in,
        constraint_tree=constraint_tree,
        trees_z=trees_z,
        suptree=suptree,
        trees_rf=trees_rf,
        prefix=prefix,
        args=args,
    )
    return await read_iqtree_results(results, prefix)


async def read_iqtree_results(results: Dir, prefix: str) -> IqtreeResult:
    """The outputs in an :data:`iqtree_cmd` results Dir written with ``prefix``."""
    files = {f.name: f for f in await results.list_files()}

    def output(suffix: str) -> File | None:
        return files.get(f"{prefix}.{suffix}")

    log = output("log")
    if log is None:
        raise FileNotFoundError(f"iqtree wrote no {prefix}.log")
    return IqtreeResult(
        phylogeny=output("treefile"),
        report=output("iqtree"),
        mldist=output("mldist"),
        lmap_svg=output("lmap.svg"),
        lmap_eps=output("lmap.eps"),
        lmap_quartetlh=output("lmap.quartetlh"),
        sitefreq_out=output("sitefreq"),
        bootstrap=output("ufboot"),
        state=output("state"),
        contree=output("contree"),
        nex=[files[name] for name in sorted(files) if name.endswith(".nex")],
        best_scheme=output("best_scheme"),
        splits=output("splits"),
        suptree=output("suptree"),
        alninfo=output("alninfo"),
        partlh=output("partlh"),
        siteprob=output("siteprob"),
        sitelh=output("sitelh"),
        treels=output("treels"),
        rate=output("rate"),
        mlrate=output("mlrate"),
        exch_matrix=files.get("GTRPMIX.nex"),
        log=log,
        results=results,
    )
