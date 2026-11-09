"""centrifuge — Centrifuge taxonomic classification indexes.

Exposes :func:`centrifuge_build`, which builds a Centrifuge index from
reference FASTA files (``centrifuge-build``), given the sequence-to-taxon
conversion table (``seqid2taxid.map``), the taxonomy tree (``nodes.dmp``),
the name table (``names.dmp``) and optionally a genome size table. The index
is a Dir of ``<prefix>.{1,2,3,4}.cf``, as upstream's ``<prefix>/``.

Several FASTA files are passed to ``centrifuge-build`` as one comma-separated
list. ``-p`` is the task's CPUs; options go in ``args``.
"""

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

CENTRIFUGE_IMAGE = "community.wave.seqera.io/library/centrifuge:1.0.4.2--d96245be9f9eb1f0"
# Upstream labels it process_high; real databases need far more memory.
DEFAULT_RESOURCES = flyte.Resources(cpu=4, memory="16Gi")

centrifuge_build_cmd = shell.create(
    name="centrifuge_build",
    image=CENTRIFUGE_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "fasta": list[File],
        "conversion_table": File,
        "taxonomy_tree": File,
        "name_table": File,
        "size_table": File | None,
        "prefix": str,
        "args": str,
    },
    defaults={"args": ""},
    outputs={"index": Dir},
    script=rf"""
        {TASK_LIMITS}
        FASTA=({{inputs.fasta}}); CONVERSION=({{inputs.conversion_table}})
        TREE=({{inputs.taxonomy_tree}}); NAMES=({{inputs.name_table}})
        shopt -s nullglob; SIZES=({{inputs.size_table}}); shopt -u nullglob
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.index}}
        W=$(mktemp -d); cd "$W"
        REFERENCES=$(IFS=,; echo "${{FASTA[*]}}")
        SIZE_TABLE=()
        if [ ${{#SIZES[@]}} -gt 0 ]; then SIZE_TABLE=(--size_table "${{SIZES[0]}}"); fi
        mkdir "$P"
        centrifuge-build \
            -p $TASK_CPUS \
            "$REFERENCES" \
            "$P/$P" \
            --conversion-table "${{CONVERSION[0]}}" \
            --taxonomy-tree "${{TREE[0]}}" \
            --name-table "${{NAMES[0]}}" \
            "${{SIZE_TABLE[@]}}" \
            $ARGS \
            || exit 1
        mv "$P"/* "$OUT"/
    """,
)

env = flyte.TaskEnvironment.from_task("centrifuge", centrifuge_build_cmd.as_task())


async def centrifuge_build(
    fasta: File | list[File],
    conversion_table: File,
    taxonomy_tree: File,
    name_table: File,
    size_table: File | None = None,
    prefix: str = "sample",
    args: str = "",
) -> Dir:
    """Build a Centrifuge index of ``fasta``: a Dir of ``<prefix>.{1,2,3,4}.cf``."""
    fastas = fasta if isinstance(fasta, list) else [fasta]
    return await centrifuge_build_cmd(
        fasta=fastas,
        conversion_table=conversion_table,
        taxonomy_tree=taxonomy_tree,
        name_table=name_table,
        size_table=size_table,
        prefix=prefix,
        args=args,
    )
