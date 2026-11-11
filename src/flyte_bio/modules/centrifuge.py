"""centrifuge — Centrifuge taxonomic classification indexes.

Exposes :func:`centrifuge_build`, which builds a Centrifuge index from
reference FASTA files (``centrifuge-build``), given the sequence-to-taxon
conversion table (``seqid2taxid.map``), the taxonomy tree (``nodes.dmp``),
the name table (``names.dmp``) and optionally a genome size table. The index
is a Dir of ``<prefix>.{1,2,3,4}.cf``, as upstream's ``<prefix>/``.

Several FASTA files are passed to ``centrifuge-build`` as one comma-separated
list. ``-p`` is the task's CPUs; options go in ``args``.

:func:`centrifuge_centrifuge` classifies one sample's reads (single- or
paired-end) against an index Dir (found by its ``*.1.cf``). Outputs keep
upstream's names:

- ``<prefix>.report.txt``: the per-taxon report (``--report-file``);
- ``<prefix>.results.txt``: the per-read classifications (``-S``);
- with ``save_aligned`` / ``save_unaligned``: the classified / unclassified
  reads, ``<prefix>.mapped.fastq.gz`` / ``<prefix>.unmapped.fastq.gz``
  (single-end) or ``<prefix>.{mapped,unmapped}.fastq.{1,2}.gz`` (paired-end).

:func:`centrifuge_kreport` turns the per-read classifications into a
Kraken-style report (``centrifuge-kreport``), ``<prefix>.txt``.
"""

from dataclasses import dataclass

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

# Upstream labels it process_high.
CLASSIFY_RESOURCES = flyte.Resources(cpu=4, memory="16Gi")

centrifuge_centrifuge_cmd = shell.create(
    name="centrifuge_centrifuge",
    image=CENTRIFUGE_IMAGE,
    resources=CLASSIFY_RESOURCES,
    inputs={
        "reads_1": File,
        "reads_2": File | None,
        "db": Dir,
        "save_unaligned": bool,
        "save_aligned": bool,
        "prefix": str,
        "args": str,
    },
    defaults={"save_unaligned": False, "save_aligned": False, "args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        R1=({{inputs.reads_1}})
        shopt -s nullglob; R2=({{inputs.reads_2}}); shopt -u nullglob
        DB={{inputs.db}}
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        DB_NAME=$(find -L "$DB" -name "*.1.cf" -not -name "._*" | sed 's/\.1.cf$//')
        if [ -z "$DB_NAME" ]; then echo "no Centrifuge index (*.1.cf) in the database" >&2; exit 1; fi
        W=$(mktemp -d); cd "$W"
        EXTRA=()
        if [ ${{#R2[@]}} -gt 0 ]; then
            INPUT=(-1 "${{R1[0]}}" -2 "${{R2[0]}}")
            if [ {{inputs.save_unaligned}} = true ]; then EXTRA+=(--un-conc-gz "$P.unmapped.fastq.gz"); fi
            if [ {{inputs.save_aligned}} = true ]; then EXTRA+=(--al-conc-gz "$P.mapped.fastq.gz"); fi
        else
            INPUT=(-U "${{R1[0]}}")
            if [ {{inputs.save_unaligned}} = true ]; then EXTRA+=(--un-gz "$P.unmapped.fastq.gz"); fi
            if [ {{inputs.save_aligned}} = true ]; then EXTRA+=(--al-gz "$P.mapped.fastq.gz"); fi
        fi
        # A directory for the pipe files, as upstream, so parallel runs don't clash in /tmp.
        mkdir ./temp
        centrifuge \
            -x "$DB_NAME" \
            --temp-directory ./temp \
            -p $TASK_CPUS \
            "${{INPUT[@]}}" \
            --report-file "$P.report.txt" \
            -S "$P.results.txt" \
            "${{EXTRA[@]}}" \
            $ARGS \
            || exit 1
        shopt -s nullglob
        mv ./*report.txt ./*results.txt ./*.sam ./*.tab ./*.mapped.fastq*.gz ./*.unmapped.fastq*.gz "$OUT"/
    """,
)

# Upstream: process_single.
KREPORT_RESOURCES = flyte.Resources(cpu=1, memory="4Gi")

centrifuge_kreport_cmd = shell.create(
    name="centrifuge_kreport",
    image=CENTRIFUGE_IMAGE,
    resources=KREPORT_RESOURCES,
    inputs={"results": File, "db": Dir, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"report": Dir},
    script=r"""
        RESULTS=({inputs.results})
        DB={inputs.db}
        P={inputs.prefix}
        ARGS={inputs.args}
        OUT={outputs.report}
        DB_NAME=$(find -L "$DB" -name "*.1.cf" -not -name "._*" | sed 's/\.1.cf$//')
        if [ -z "$DB_NAME" ]; then echo "no Centrifuge index (*.1.cf) in the database" >&2; exit 1; fi
        centrifuge-kreport \
            $ARGS \
            -x "$DB_NAME" \
            "${RESULTS[0]}" > "$OUT/$P.txt" \
            || exit 1
    """,
)

env = flyte.TaskEnvironment.from_task(
    "centrifuge",
    centrifuge_build_cmd.as_task(),
    centrifuge_centrifuge_cmd.as_task(),
    centrifuge_kreport_cmd.as_task(),
)


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


@dataclass
class CentrifugeResult:
    report: File  # <prefix>.report.txt
    results: File  # <prefix>.results.txt
    mapped: list[File]  # with save_aligned: one file, or mates 1 and 2
    unmapped: list[File]  # with save_unaligned: one file, or mates 1 and 2
    all_results: Dir  # everything centrifuge wrote


async def centrifuge_centrifuge(
    reads_1: File,
    db: Dir,
    reads_2: File | None = None,
    save_unaligned: bool = False,
    save_aligned: bool = False,
    prefix: str = "sample",
    args: str = "",
) -> CentrifugeResult:
    """Classify one sample's reads (``reads_2`` None for single-end) against the Centrifuge index ``db``."""
    results = await centrifuge_centrifuge_cmd(
        reads_1=reads_1,
        reads_2=reads_2,
        db=db,
        save_unaligned=save_unaligned,
        save_aligned=save_aligned,
        prefix=prefix,
        args=args,
    )
    report = await results.get_file(f"{prefix}.report.txt")
    per_read = await results.get_file(f"{prefix}.results.txt")
    if report is None or per_read is None:
        raise FileNotFoundError(f"centrifuge wrote no {prefix}.report.txt or {prefix}.results.txt")

    async def reads(kind: str) -> list[File]:
        names = [f"{prefix}.{kind}.fastq.gz", f"{prefix}.{kind}.fastq.1.gz", f"{prefix}.{kind}.fastq.2.gz"]
        found = [await results.get_file(name) for name in names]
        return [f for f in found if f is not None]

    return CentrifugeResult(
        report=report,
        results=per_read,
        mapped=await reads("mapped"),
        unmapped=await reads("unmapped"),
        all_results=results,
    )


async def centrifuge_kreport(results: File, db: Dir, prefix: str = "sample", args: str = "") -> File:
    """A Kraken-style report of centrifuge's per-read ``results`` against ``db``: ``<prefix>.txt``."""
    out = await centrifuge_kreport_cmd(results=results, db=db, prefix=prefix, args=args)
    report = await out.get_file(f"{prefix}.txt")
    if report is None:
        raise FileNotFoundError(f"centrifuge-kreport wrote no {prefix}.txt")
    return report
