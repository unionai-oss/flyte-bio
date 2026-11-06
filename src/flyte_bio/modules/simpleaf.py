"""simpleaf — single-cell quantification with alevin-fry.

Exposes:

- :func:`simpleaf_index` — a piscem index of the spliced and intronic
  transcriptome ("spliceu") built from a genome FASTA + GTF, or a plain
  index of a transcript FASTA. Written to ``<prefix>/index`` and, from a
  genome, the reference it was built from to ``<prefix>/ref`` (with
  ``t2g_3col.tsv``), as upstream. simpleaf's other files there
  (``index_info.json``, the log) are dropped, as upstream emits only these
  two: a Dir input is staged by key prefix, so ``index_info.json`` beside
  ``index`` would be pulled into it.
- :func:`simpleaf_quant` — map one sample's reads and quantify them per cell:
  ``<prefix>/af_map`` (the mapping) and ``<prefix>/af_quant`` (with the count
  matrix in ``alevin/``, also as ``quants.h5ad``).

simpleaf keeps its configuration in ``ALEVIN_FRY_HOME``, so both tasks run
in a scratch working directory.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

SIMPLEAF_IMAGE = "quay.io/biocontainers/simpleaf:0.19.5--ha6fb395_0"

SIMPLEAF_CPUS = 4
INDEX_RESOURCES = flyte.Resources(cpu=SIMPLEAF_CPUS, memory="24Gi")
QUANT_RESOURCES = flyte.Resources(cpu=SIMPLEAF_CPUS, memory="16Gi")

CELL_FILTERS = ("unfiltered-pl", "explicit-pl", "knee", "expect-cells", "forced-cells")

# As upstream: a transcript FASTA takes precedence over genome FASTA + GTF.
simpleaf_index_cmd = shell.create(
    name="simpleaf_index",
    image=SIMPLEAF_IMAGE,
    resources=INDEX_RESOURCES,
    inputs={
        "genome_fasta": File | None,
        "gtf": File | None,
        "transcript_fasta": File | None,
        "prefix": str,
        "args": str,
    },
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        shopt -s nullglob
        FA=({{inputs.genome_fasta}}); GTF=({{inputs.gtf}}); TX=({{inputs.transcript_fasta}})
        shopt -u nullglob
        ARGS={{inputs.args}}
        if [ ${{#TX[@]}} -gt 0 ]; then SEQ=(--ref-seq "${{TX[0]}}")
        else SEQ=(--fasta "${{FA[0]}}" --gtf "${{GTF[0]}}"); fi
        W=$(mktemp -d); cd "$W"
        export ALEVIN_FRY_HOME=.
        ulimit -n 2048
        simpleaf set-paths
        simpleaf index --threads $TASK_CPUS "${{SEQ[@]}}" $ARGS -o out
        mkdir -p {{outputs.results}}/{{inputs.prefix}}
        for d in index ref; do if [ -d out/$d ]; then mv out/$d {{outputs.results}}/{{inputs.prefix}}/; fi; done
    """,
)

# Reads come as two comma-separated lists, lane by lane. The order is passed
# as file names (list inputs are staged by name, not in list order).
simpleaf_quant_cmd = shell.create(
    name="simpleaf_quant",
    image=SIMPLEAF_IMAGE,
    resources=QUANT_RESOURCES,
    inputs={
        "reads_1": list[File],
        "reads_2": list[File],
        "order_1": str,
        "order_2": str,
        "index": Dir,
        "t2g": File | None,
        "chemistry": str,
        "resolution": str,
        "cell_filter": str,
        "permit_list": File | None,
        "prefix": str,
        "args": str,
    },
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        R1=({{inputs.reads_1}}); R2=({{inputs.reads_2}})
        shopt -s nullglob; T2G=({{inputs.t2g}}); PL=({{inputs.permit_list}}); shopt -u nullglob
        DIR1=$(dirname "${{R1[0]}}"); DIR2=$(dirname "${{R2[0]}}")
        READS1=$(echo {{inputs.order_1}} | tr ',' '\n' | sed "s#^#$DIR1/#" | paste -sd, -)
        READS2=$(echo {{inputs.order_2}} | tr ',' '\n' | sed "s#^#$DIR2/#" | paste -sd, -)
        T2G_ARG=(); [ ${{#T2G[@]}} -gt 0 ] && T2G_ARG=(--t2g-map "${{T2G[0]}}")
        FILTER=({{inputs.cell_filter}})
        [ ${{#PL[@]}} -gt 0 ] && FILTER+=("${{PL[0]}}")
        ARGS={{inputs.args}}
        INDEX={{inputs.index}}
        W=$(mktemp -d); cd "$W"
        export ALEVIN_FRY_HOME=.
        simpleaf set-paths
        simpleaf quant "${{T2G_ARG[@]}}" --chemistry {{inputs.chemistry}} --index "$INDEX" \
            --reads1 "$READS1" --reads2 "$READS2" --resolution {{inputs.resolution}} \
            --output {{outputs.results}}/{{inputs.prefix}} --threads $TASK_CPUS --anndata-out \
            ${{FILTER[@]}} $ARGS
    """,
)


env = flyte.TaskEnvironment.from_task(
    "simpleaf",
    simpleaf_index_cmd.as_task(),
    simpleaf_quant_cmd.as_task(),
)


def subdir(parent: Dir, rel: str) -> Dir:
    return Dir.from_existing_remote(parent.path.rstrip("/") + "/" + rel)


@dataclass
class SimpleafIndex:
    index: Dir  # <prefix>/index
    ref: Dir | None  # <prefix>/ref, when built from a genome
    t2g: File | None  # <prefix>/ref/t2g_3col.tsv (or t2g.tsv), when built from a genome


async def simpleaf_index(
    prefix: str,
    genome_fasta: File | None = None,
    gtf: File | None = None,
    transcript_fasta: File | None = None,
    args: str = "",
) -> SimpleafIndex:
    """A simpleaf index of ``transcript_fasta``, or else of the spliced + intronic transcriptome of FASTA + GTF."""
    if transcript_fasta is None and (genome_fasta is None or gtf is None):
        raise ValueError("simpleaf_index needs a transcript_fasta, or a genome_fasta and gtf")
    results = await simpleaf_index_cmd(
        genome_fasta=genome_fasta, gtf=gtf, transcript_fasta=transcript_fasta, prefix=prefix, args=args
    )
    index = subdir(results, f"{prefix}/index")
    if transcript_fasta is not None:
        return SimpleafIndex(index=index, ref=None, t2g=None)
    ref = subdir(results, f"{prefix}/ref")
    t2g = await ref.get_file("t2g_3col.tsv") or await ref.get_file("t2g.tsv")
    return SimpleafIndex(index=index, ref=ref, t2g=t2g)


@dataclass
class SimpleafQuant:
    map: Dir  # <prefix>/af_map
    quant: Dir  # <prefix>/af_quant (count matrix in alevin/)


def file_name(f: File) -> str:
    return f.path.rstrip("/").rsplit("/", 1)[-1]


async def simpleaf_quant(
    reads_1: list[File],
    reads_2: list[File],
    index: Dir,
    chemistry: str,
    prefix: str,
    t2g: File | None = None,
    resolution: str = "cr-like",
    cell_filter: str = "unfiltered-pl",
    cell_filter_number: int | None = None,
    permit_list: File | None = None,
    args: str = "",
) -> SimpleafQuant:
    """Map and quantify one sample (its lanes' R1 and R2 FASTQs, paired in order) with simpleaf.

    ``cell_filter`` is upstream's: ``unfiltered-pl`` with the protocol's
    barcode ``permit_list`` (simpleaf's own for 10x chemistries when None),
    ``explicit-pl`` with a ``permit_list``, ``knee``, or ``expect-cells`` /
    ``forced-cells`` with ``cell_filter_number``.
    """
    method = cell_filter.replace("_", "-")
    if method not in CELL_FILTERS:
        raise ValueError(f"cell_filter must be one of {CELL_FILTERS}, got {cell_filter!r}")
    if len(reads_1) != len(reads_2) or not reads_1:
        raise ValueError("simpleaf_quant needs one R2 per R1")
    filter_args = f"--{method}"
    if method in ("expect-cells", "forced-cells"):
        if cell_filter_number is None:
            raise ValueError(f"cell_filter {method!r} needs cell_filter_number")
        filter_args += f" {cell_filter_number}"
    if method == "explicit-pl" and permit_list is None:
        raise ValueError("cell_filter 'explicit-pl' needs a permit_list")
    names_1, names_2 = [file_name(f) for f in reads_1], [file_name(f) for f in reads_2]
    if len(set(names_1)) != len(names_1) or len(set(names_2)) != len(names_2):
        raise ValueError("a sample's FASTQs must have distinct file names")
    results = await simpleaf_quant_cmd(
        reads_1=reads_1,
        reads_2=reads_2,
        order_1=",".join(names_1),
        order_2=",".join(names_2),
        index=index,
        t2g=t2g,
        chemistry=chemistry,
        resolution=resolution,
        cell_filter=filter_args,
        permit_list=permit_list if method in ("unfiltered-pl", "explicit-pl") else None,
        prefix=prefix,
        args=args,
    )
    return SimpleafQuant(map=subdir(results, f"{prefix}/af_map"), quant=subdir(results, f"{prefix}/af_quant"))
