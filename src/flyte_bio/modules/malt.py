"""malt — MALT (MEGAN alignment tool) indexes for metagenomic alignment.

Exposes :func:`malt_build`, which builds a MALT index (``malt-build``) from
reference FASTA files, linked into ``fa_folder/`` as upstream, optionally
with a GFF annotation (``-igff``), and a MEGAN mapping file that assigns the
references to taxa: ``mapping_db`` passed as ``-<map_type>``, e.g. a MEGAN
mapping database (``map_type="mdb"``, ``megan-*.db``) or an
accession-to-taxon binary (``map_type="a2t"``, ``*.abin``). It returns the
index, upstream's ``malt_index/`` (``index0.idx``, ``ref.{db,idx,inf}``,
``table0.{db,idx}``, ``taxonomy.{idx,map,tre}``, ...), and the build log
``malt-build.log``.

MALT is Java: the JVM heap is 80% of the task's memory (``-J-Xmx``), and a
``-J-Xmx`` in ``args`` comes later and wins. ``-t`` is the task's CPUs; the
sequence type (``--sequenceType DNA``) and other options go in ``args``.
Mapping files are several GiB, so the default resources include disk for
them.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import JAVA_HEAP_MB, TASK_LIMITS

MALT_IMAGE = "quay.io/biocontainers/malt:0.62--hdfd78af_0"
# Upstream labels it process_high; real databases need far more memory.
DEFAULT_RESOURCES = flyte.Resources(cpu=4, memory="16Gi", disk="50Gi")

malt_build_cmd = shell.create(
    name="malt_build",
    image=MALT_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"fasta": list[File], "gff": File | None, "mapping_db": File, "map_type": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        FASTA=({{inputs.fasta}}); MAPPING=({{inputs.mapping_db}})
        shopt -s nullglob; GFF=({{inputs.gff}}); shopt -u nullglob
        MAP_TYPE={{inputs.map_type}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        mkdir fa_folder
        for F in "${{FASTA[@]}}"; do ln -s "$F" fa_folder/"$(basename "$F")"; done
        IGFF=()
        if [ ${{#GFF[@]}} -gt 0 ]; then IGFF=(-igff "${{GFF[0]}}"); fi
        set -o pipefail
        malt-build \
            -J-Xmx{JAVA_HEAP_MB}m \
            $ARGS \
            -v \
            --input fa_folder/ \
            "${{IGFF[@]}}" \
            -d malt_index/ \
            -t $TASK_CPUS \
            -"$MAP_TYPE" "${{MAPPING[0]}}" \
            |& tee malt-build.log \
            || exit 1
        mv malt_index malt-build.log "$OUT"/
    """,
)

env = flyte.TaskEnvironment.from_task("malt", malt_build_cmd.as_task())


@dataclass
class MaltBuildResult:
    index: Dir  # malt_index/
    log: File  # malt-build.log


async def malt_build(
    fasta: File | list[File],
    mapping_db: File,
    map_type: str,
    gff: File | None = None,
    args: str = "",
) -> MaltBuildResult:
    """Build a MALT index of ``fasta`` with the taxon mapping ``mapping_db`` (``-<map_type>``)."""
    fastas = fasta if isinstance(fasta, list) else [fasta]
    results = await malt_build_cmd(fasta=fastas, gff=gff, mapping_db=mapping_db, map_type=map_type, args=args)
    log = await results.get_file("malt-build.log")
    if log is None:
        raise FileNotFoundError("malt-build wrote no malt-build.log")
    return MaltBuildResult(index=Dir.from_existing_remote(results.path.rstrip("/") + "/malt_index"), log=log)
