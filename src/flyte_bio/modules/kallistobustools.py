"""kallistobustools — kb-python: kallisto | bustools single-cell references and counts.

Exposes:

- :func:`kallistobustools_ref`, which builds a kallisto index and its
  transcript-to-gene map from a genome FASTA and a GTF (``kb ref``). Outputs
  keep upstream's names: ``kb_ref_out.idx``, ``t2g.txt`` and ``cdna.fa``,
  plus, for workflows other than ``standard`` (``nac``, ``lamanno``),
  ``intron.fa``, ``cdna_t2c.txt`` and ``intron_t2c.txt``.
- :func:`kallistobustools_count`, which pseudoaligns one sample's reads and
  counts cells by genes (``kb count``) with such an index, for a technology
  (``-x``: ``10XV3``, ...) and workflow. As upstream, the inputs are staged
  under their own names and passed by name, so ``kb_info.json`` records the
  same commands. The output is kb's ``<prefix>.count`` folder
  (``counts_unfiltered/cells_x_genes.mtx`` and its barcodes and genes,
  ``output.bus``, ``matrix.ec``, ``transcripts.txt``, ``inspect.json``,
  ``kb_info.json``, ``run_info.json``, the technology's barcode whitelist,
  ...).

``kb count`` gets the task's CPUs (``-t``) and its memory less 1 GB
(``-m``), as upstream.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

KB_IMAGE = "quay.io/biocontainers/kb-python:0.28.2--pyhdfd78af_2"
# Upstream: process_medium (6 CPUs, 36 GB). A human genome's index wants
# about 32 GiB: raise with <cmd>.override(resources=...).
REF_RESOURCES = flyte.Resources(cpu=4, memory="16Gi", disk="50Gi")
COUNT_RESOURCES = flyte.Resources(cpu=4, memory="16Gi", disk="50Gi")

kallistobustools_ref_cmd = shell.create(
    name="kallistobustools_ref",
    image=KB_IMAGE,
    resources=REF_RESOURCES,
    inputs={"fasta": File, "gtf": File, "workflow_mode": str},
    outputs={"results": Dir},
    script=r"""
        FA=({inputs.fasta}); GTF=({inputs.gtf})
        MODE={inputs.workflow_mode}
        OUT={outputs.results}
        W=$(mktemp -d); cd "$W"
        ln -s "${FA[0]}" "${GTF[0]}" .
        FASTA=$(basename "${FA[0]}"); GENES=$(basename "${GTF[0]}")
        if [ "$MODE" = standard ]; then
            kb \
                ref \
                -i kb_ref_out.idx \
                -g t2g.txt \
                -f1 cdna.fa \
                --workflow $MODE \
                "$FASTA" \
                "$GENES" || exit 1
            mv kb_ref_out.idx t2g.txt cdna.fa "$OUT"/
        else
            kb \
                ref \
                -i kb_ref_out.idx \
                -g t2g.txt \
                -f1 cdna.fa \
                -f2 intron.fa \
                -c1 cdna_t2c.txt \
                -c2 intron_t2c.txt \
                --workflow $MODE \
                "$FASTA" \
                "$GENES" || exit 1
            mv kb_ref_out.idx t2g.txt cdna.fa intron.fa cdna_t2c.txt intron_t2c.txt "$OUT"/
        fi
    """,
)

kallistobustools_count_cmd = shell.create(
    name="kallistobustools_count",
    image=KB_IMAGE,
    resources=COUNT_RESOURCES,
    inputs={
        "reads": list[File],
        "index": File,
        "t2g": File,
        "t1c": File | None,
        "t2c": File | None,
        "technology": str,
        "workflow_mode": str,
        "prefix": str,
        "args": str,
    },
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        shopt -s nullglob
        READS=({{inputs.reads}}); INDEX=({{inputs.index}}); T2G=({{inputs.t2g}})
        T1C=({{inputs.t1c}}); T2C=({{inputs.t2c}})
        shopt -u nullglob
        TECH={{inputs.technology}}
        MODE={{inputs.workflow_mode}}
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"

        # Every input under its own name in the working directory, passed by name as upstream.
        for F in "${{READS[@]}}" "${{INDEX[@]}}" "${{T2G[@]}}" "${{T1C[@]}}" "${{T2C[@]}}"; do ln -s "$F" .; done
        names() {{ for F in "$@"; do basename "$F"; done; }}
        CDNA=(); INTRON=()
        if [ ${{#T1C[@]}} -gt 0 ]; then CDNA=(-c1 "$(basename "${{T1C[0]}}")"); fi
        if [ ${{#T2C[@]}} -gt 0 ]; then INTRON=(-c2 "$(basename "${{T2C[0]}}")"); fi
        mapfile -t READ_NAMES < <(names "${{READS[@]}}")

        kb \
            count \
            -t $TASK_CPUS \
            -i "$(basename "${{INDEX[0]}}")" \
            -g "$(basename "${{T2G[0]}}")" \
            "${{CDNA[@]}}" \
            "${{INTRON[@]}}" \
            -x $TECH \
            --workflow $MODE \
            $ARGS \
            -o "$P.count" \
            -m $(( TASK_MEMORY_MB / 1024 - 1 ))G \
            "${{READ_NAMES[@]}}" || exit 1
        mv "$P.count" "$OUT"/
    """,
)

env = flyte.TaskEnvironment.from_task(
    "kallistobustools", kallistobustools_ref_cmd.as_task(), kallistobustools_count_cmd.as_task()
)


@dataclass
class KbReference:
    index: File  # kb_ref_out.idx
    t2g: File  # t2g.txt
    cdna: File  # cdna.fa
    intron: File | None  # intron.fa (not with the standard workflow)
    cdna_t2c: File | None  # cdna_t2c.txt (not with the standard workflow)
    intron_t2c: File | None  # intron_t2c.txt (not with the standard workflow)


async def kallistobustools_ref(fasta: File, gtf: File, workflow_mode: str = "standard") -> KbReference:
    """A kb index and transcript-to-gene map from ``fasta`` and ``gtf`` for ``workflow_mode``."""
    results = await kallistobustools_ref_cmd(fasta=fasta, gtf=gtf, workflow_mode=workflow_mode)

    async def pick(name: str) -> File:
        found = await results.get_file(name)
        if found is None:
            raise FileNotFoundError(f"kb ref wrote no {name}")
        return found

    return KbReference(
        index=await pick("kb_ref_out.idx"),
        t2g=await pick("t2g.txt"),
        cdna=await pick("cdna.fa"),
        intron=await results.get_file("intron.fa"),
        cdna_t2c=await results.get_file("cdna_t2c.txt"),
        intron_t2c=await results.get_file("intron_t2c.txt"),
    )


@dataclass
class KbCountResult:
    count: Dir  # <prefix>.count
    matrices: list[File]  # <prefix>.count/*/*.mtx


async def kallistobustools_count(
    reads: list[File],
    index: File,
    t2g: File,
    technology: str,
    workflow_mode: str = "standard",
    t1c: File | None = None,
    t2c: File | None = None,
    prefix: str = "sample",
    args: str = "",
) -> KbCountResult:
    """``kb count`` of one sample's ``reads`` with ``index`` and ``t2g`` for ``technology`` (``-x``)."""
    if not reads:
        raise ValueError("kallistobustools_count needs the sample's reads")
    results = await kallistobustools_count_cmd(
        reads=reads,
        index=index,
        t2g=t2g,
        t1c=t1c,
        t2c=t2c,
        technology=technology,
        workflow_mode=workflow_mode,
        prefix=prefix,
        args=args,
    )
    return await read_count_results(results, prefix)


async def read_count_results(results: Dir, prefix: str) -> KbCountResult:
    """Pick the outputs out of ``kallistobustools_count_cmd``'s results Dir (e.g. from an overridden task)."""
    count = Dir.from_existing_remote(results.path.rstrip("/") + f"/{prefix}.count")
    root = count.path.rstrip("/") + "/"
    matrices = []
    async for file in count.walk():
        relative = file.path.removeprefix(root)
        if relative.count("/") == 1 and relative.endswith(".mtx"):
            matrices.append(file)
    if not matrices:
        raise FileNotFoundError(f"kb count wrote no {prefix}.count/*/*.mtx")
    return KbCountResult(count=count, matrices=sorted(matrices, key=lambda f: f.path))
