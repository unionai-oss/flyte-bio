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

:func:`malt_run` aligns reads (FASTQ or FASTA files) to a MALT index
(``malt-run``) and writes, per input file, MEGAN's ``<input name>.rma6``
and, with ``-a .`` in ``args``, the alignments in the format ``-f`` names
(e.g. ``<input name>.blastn.sam.gz``), plus the log
``<prefix>-malt-run.log``. The JVM heap and ``-t`` are set as for
:func:`malt_build`; the alignment mode (``-m BlastN``) and other options go
in ``args``.
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

malt_run_cmd = shell.create(
    name="malt_run",
    image=MALT_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"reads": list[File], "reads_order": str, "index": Dir, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        READS=({{inputs.reads}})
        # Staged list inputs expand in name order; restore the order they were given in.
        read -ra READS_NAMES <<< {{inputs.reads_order}}; READS_DIR=$(dirname "${{READS[0]}}")
        READS=(); for N in "${{READS_NAMES[@]}}"; do READS+=("$READS_DIR/$N"); done
        INDEX={{inputs.index}}
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        # The SAM header records the command line: give it upstream's relative paths.
        ln -s "$INDEX" malt_index
        for F in "${{READS[@]}}"; do ln -s "$F" "$(basename "$F")"; done
        set -o pipefail
        malt-run \
            -J-Xmx{JAVA_HEAP_MB}m \
            -t $TASK_CPUS \
            -v \
            -o . \
            $ARGS \
            --inFile "${{READS_NAMES[@]}}" \
            --index malt_index/ \
            |& tee "$P-malt-run.log" \
            || exit 1
        shopt -s nullglob
        mv ./*.rma6 ./*.tab ./*.text ./*.sam ./*.tab.gz ./*.text.gz ./*.sam.gz ./*-malt-run.log "$OUT"/
    """,
)

env = flyte.TaskEnvironment.from_task("malt", malt_build_cmd.as_task(), malt_run_cmd.as_task())


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


@dataclass
class MaltRunResult:
    rma6: list[File]  # <input name>.rma6, one per input file
    alignments: list[File]  # with -a: <input name>.<mode>.{tab,text,sam}[.gz]
    log: File  # <prefix>-malt-run.log
    results: Dir  # everything above


async def malt_run(reads: File | list[File], index: Dir, prefix: str = "sample", args: str = "") -> MaltRunResult:
    """Align ``reads`` to the MALT ``index`` (:func:`malt_build`'s)."""
    files = reads if isinstance(reads, list) else [reads]
    results = await malt_run_cmd(
        reads=files, reads_order=" ".join(f.name for f in files), index=index, prefix=prefix, args=args
    )
    return await malt_run_result(results, prefix)


async def malt_run_result(results: Dir, prefix: str) -> MaltRunResult:
    """Pick malt-run's outputs out of ``results``."""
    log = await results.get_file(f"{prefix}-malt-run.log")
    if log is None:
        raise FileNotFoundError(f"malt-run wrote no {prefix}-malt-run.log")
    everything = await results.list_files()
    alignment_suffixes = (".tab", ".text", ".sam", ".tab.gz", ".text.gz", ".sam.gz")
    return MaltRunResult(
        rma6=sorted((f for f in everything if f.name.endswith(".rma6")), key=lambda f: f.name),
        alignments=sorted((f for f in everything if f.name.endswith(alignment_suffixes)), key=lambda f: f.name),
        log=log,
        results=results,
    )
