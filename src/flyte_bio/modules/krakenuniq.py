"""krakenuniq — KrakenUniq taxonomic classification databases.

Exposes :func:`krakenuniq_build`, which builds a KrakenUniq database
(``krakenuniq-build``) from library FASTA files, the taxonomy files
(``names.dmp`` and ``nodes.dmp``, as ``krakenuniq-download taxonomy``
fetches) and a ``seqid2taxid.map``. As upstream, they're copied into
``<prefix>/library/``, ``<prefix>/taxonomy/`` and
``<prefix>/seqid2taxid.map`` for the build and removed after it. Without
``keep_intermediate`` only the files KrakenUniq classifies with remain:
``database*.kdb``, ``database.idx``, ``database.kdb.counts`` and ``taxDB``.
The database comes back as a Dir of them.

``--threads`` is the task's CPUs. As krakenuniq-build otherwise looks for
Jellyfish next to itself, ``--jellyfish-bin`` points it at the image's
``jellyfish`` (upstream's test passes the same in ``ext.args``); other
options go in ``args``.
"""

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

KRAKENUNIQ_IMAGE = "quay.io/biocontainers/krakenuniq:1.0.4--pl5321h6dccd9a_2"
# Upstream labels it process_medium; real databases need far more memory.
DEFAULT_RESOURCES = flyte.Resources(cpu=4, memory="8Gi")

krakenuniq_build_cmd = shell.create(
    name="krakenuniq_build",
    image=KRAKENUNIQ_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "library": list[File],
        "taxonomy": list[File],
        "seqid2taxid": File,
        "keep_intermediate": bool,
        "prefix": str,
        "args": str,
    },
    defaults={"keep_intermediate": False, "args": ""},
    outputs={"db": Dir},
    script=rf"""
        {TASK_LIMITS}
        LIBRARY=({{inputs.library}}); TAXONOMY=({{inputs.taxonomy}}); MAP=({{inputs.seqid2taxid}})
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.db}}
        W=$(mktemp -d); cd "$W"
        mkdir -p "$P/library" "$P/taxonomy"
        for F in "${{LIBRARY[@]}}"; do cp "$F" "$P/library/$(basename "$F")"; done
        for F in "${{TAXONOMY[@]}}"; do cp "$F" "$P/taxonomy/$(basename "$F")"; done
        cp "${{MAP[0]}}" "$P/seqid2taxid.map"
        krakenuniq-build \
            --jellyfish-bin "$(command -v jellyfish)" \
            $ARGS \
            --threads $TASK_CPUS \
            --db "$P" \
            || exit 1
        rm -r "$P/library" "$P/taxonomy" "$P/seqid2taxid.map"
        if [ {{inputs.keep_intermediate}} != true ]; then
            find -L "$P" -type f \
                -not -name "*.kdb" \
                -not -name "*idx" \
                -not -name "taxDB" \
                -not -name "*.counts" \
                -delete
            find "$P" -type d -empty -delete
        fi
        mv "$P"/* "$OUT"/
    """,
)

env = flyte.TaskEnvironment.from_task("krakenuniq", krakenuniq_build_cmd.as_task())


async def krakenuniq_build(
    library: File | list[File],
    taxonomy: list[File],
    seqid2taxid: File,
    keep_intermediate: bool = False,
    prefix: str = "sample",
    args: str = "",
) -> Dir:
    """Build a KrakenUniq database of ``library`` over ``taxonomy`` (``names.dmp``, ``nodes.dmp``)."""
    files = library if isinstance(library, list) else [library]
    return await krakenuniq_build_cmd(
        library=files,
        taxonomy=taxonomy,
        seqid2taxid=seqid2taxid,
        keep_intermediate=keep_intermediate,
        prefix=prefix,
        args=args,
    )
