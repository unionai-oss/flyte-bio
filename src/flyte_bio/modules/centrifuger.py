"""centrifuger — Centrifuger taxonomic classification indexes.

Exposes :func:`centrifuger_build`, which builds a Centrifuger index from
reference FASTA files (``centrifuger-build``), given the taxonomy tree
(``nodes.dmp``), the name table (``names.dmp``) and the sequence-to-taxon
conversion table (``seqid2taxid.map``), which Centrifuger always needs. The
index is a Dir of ``<prefix>.{1,2,3,4}.cfr``, as upstream's ``<prefix>/``.

As upstream, the references are linked into ``genomes/`` and listed with
``find`` for ``-l``. ``-t`` is the task's CPUs; options go in ``args``.
"""

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

CENTRIFUGER_IMAGE = "quay.io/biocontainers/centrifuger:1.1.2--h3be2455_0"
# Upstream labels it process_high; real databases need far more memory.
DEFAULT_RESOURCES = flyte.Resources(cpu=4, memory="16Gi")

centrifuger_build_cmd = shell.create(
    name="centrifuger_build",
    image=CENTRIFUGER_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "references": list[File],
        "taxonomy_nodes": File,
        "taxonomy_names": File,
        "conversion_table": File,
        "prefix": str,
        "args": str,
    },
    defaults={"args": ""},
    outputs={"index": Dir},
    script=rf"""
        {TASK_LIMITS}
        REFERENCES=({{inputs.references}}); NODES=({{inputs.taxonomy_nodes}})
        NAMES=({{inputs.taxonomy_names}}); CONVERSION=({{inputs.conversion_table}})
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.index}}
        W=$(mktemp -d); cd "$W"
        mkdir genomes
        for F in "${{REFERENCES[@]}}"; do ln -s "$F" genomes/"$(basename "$F")"; done
        find -L genomes/ -type f > reference_list.txt
        mkdir -p "$P"
        centrifuger-build \
            -l reference_list.txt \
            --taxonomy-tree "${{NODES[0]}}" \
            --name-table "${{NAMES[0]}}" \
            --conversion-table "${{CONVERSION[0]}}" \
            -t $TASK_CPUS \
            -o "$P/$P" \
            $ARGS \
            || exit 1
        mv "$P"/* "$OUT"/
    """,
)

env = flyte.TaskEnvironment.from_task("centrifuger", centrifuger_build_cmd.as_task())


async def centrifuger_build(
    references: File | list[File],
    taxonomy_nodes: File,
    taxonomy_names: File,
    conversion_table: File,
    prefix: str = "sample",
    args: str = "",
) -> Dir:
    """Build a Centrifuger index of ``references``: a Dir of ``<prefix>.{1,2,3,4}.cfr``."""
    files = references if isinstance(references, list) else [references]
    return await centrifuger_build_cmd(
        references=files,
        taxonomy_nodes=taxonomy_nodes,
        taxonomy_names=taxonomy_names,
        conversion_table=conversion_table,
        prefix=prefix,
        args=args,
    )
