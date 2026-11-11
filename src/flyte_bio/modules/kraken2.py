"""kraken2 — taxonomic classification of reads against a Kraken2 database.

Exposes :func:`kraken2`, which classifies one sample's reads (single-end, or
paired-end with ``reads_2``) against ``db``, a Dir holding the database's
``hash.k2d``, ``opts.k2d`` and ``taxo.k2d``. Outputs keep upstream's names:

- the report: ``<prefix>.kraken2.report.txt`` (always);
- with ``save_output_fastqs``: the classified and unclassified reads,
  ``<prefix>.classified.fastq.gz`` / ``<prefix>.unclassified.fastq.gz``
  (single-end) or ``<prefix>.classified_{1,2}.fastq.gz`` /
  ``<prefix>.unclassified_{1,2}.fastq.gz`` (paired-end), compressed with pigz;
- with ``save_reads_assignment``: the per-read classifications,
  ``<prefix>.kraken2.classifiedreads.txt``.

kraken2 loads the whole database into memory, so a real database (the
standard one is tens of GiB) needs ``memory`` raised with
``kraken2_cmd.override(resources=...)``; the default fits small
databases like the test one.

Building a custom database takes two steps, as upstream:

- :func:`kraken2_add` adds FASTA files to a new database's library
  (``kraken2-build --add-to-library``, once per file), next to the taxonomy
  (``names.dmp``, ``nodes.dmp`` and any ``*.accession2taxid`` files) and an
  optional custom ``seqid2taxid.map``. A library holds a few files per FASTA,
  so the database directory ``<prefix>/`` (``library/added/``,
  ``taxonomy/``, ``seqid2taxid.map``) comes back as ``<prefix>.tar.gz``.
- :func:`kraken2_build` builds the database from that archive
  (``kraken2-build --build``). It returns the database's top-level files,
  ``hash.k2d``, ``opts.k2d``, ``taxo.k2d`` (plus ``seqid2taxid.map`` and
  ``unmapped.txt`` where written): a Dir :func:`kraken2` reads. With
  ``clean`` (upstream's ``cleaning``) it then removes the library and
  taxonomy (``kraken2-build --clean``); without, the whole database
  directory also comes back as ``kraken2-database.tar.gz``, which
  :func:`flyte_bio.modules.bracken.bracken_build` needs.

The thread count is the task's CPUs.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

KRAKEN2_IMAGE = "community.wave.seqera.io/library/kraken2_coreutils_pigz:920ecc6b96e2ba71"
DEFAULT_RESOURCES = flyte.Resources(cpu=4, memory="8Gi")

kraken2_cmd = shell.create(
    name="kraken2",
    image=KRAKEN2_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "reads_1": File,
        "reads_2": File | None,
        "db": Dir,
        "prefix": str,
        "save_output_fastqs": bool,
        "save_reads_assignment": bool,
        "args": str,
    },
    defaults={"save_output_fastqs": False, "save_reads_assignment": False, "args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        R1=({{inputs.reads_1}})
        shopt -s nullglob; R2=({{inputs.reads_2}}); shopt -u nullglob
        DB={{inputs.db}}
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"

        READS=("${{R1[0]}}")
        EXTRA=()
        if [ ${{#R2[@]}} -gt 0 ]; then
            READS+=("${{R2[0]}}")
            EXTRA+=(--paired)
            CLASSIFIED="$P.classified#.fastq"
            UNCLASSIFIED="$P.unclassified#.fastq"
        else
            CLASSIFIED="$P.classified.fastq"
            UNCLASSIFIED="$P.unclassified.fastq"
        fi
        SAVE_FASTQS={{inputs.save_output_fastqs}}
        if [ $SAVE_FASTQS = true ]; then
            EXTRA+=(--classified-out "$CLASSIFIED" --unclassified-out "$UNCLASSIFIED")
        fi
        if [ {{inputs.save_reads_assignment}} = true ]; then
            EXTRA+=(--output "$P.kraken2.classifiedreads.txt")
        else
            EXTRA+=(--output /dev/null)
        fi

        kraken2 \
            --db "$DB" \
            --threads $TASK_CPUS \
            --report "$P.kraken2.report.txt" \
            --gzip-compressed \
            "${{EXTRA[@]}}" \
            $ARGS \
            "${{READS[@]}}"

        if [ $SAVE_FASTQS = true ]; then pigz -p $TASK_CPUS ./*.fastq; fi
        shopt -s nullglob
        mv ./*.kraken2.report.txt ./*.classifiedreads.txt ./*.fastq.gz "$OUT"/
    """,
)

# Upstream labels both builds process_medium.
BUILD_RESOURCES = flyte.Resources(cpu=4, memory="8Gi")

kraken2_add_cmd = shell.create(
    name="kraken2_add",
    image=KRAKEN2_IMAGE,
    resources=BUILD_RESOURCES,
    inputs={
        "fasta": list[File],
        "taxonomy_names": File,
        "taxonomy_nodes": File,
        "accession2taxid": list[File],
        "seqid2taxid": File | None,
        "prefix": str,
        "args": str,
    },
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        FASTA=({{inputs.fasta}}); NAMES=({{inputs.taxonomy_names}}); NODES=({{inputs.taxonomy_nodes}})
        shopt -s nullglob; ACCESSIONS=({{inputs.accession2taxid}}); MAP=({{inputs.seqid2taxid}}); shopt -u nullglob
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"

        mkdir -p "$P/taxonomy"
        ln -s "${{NAMES[0]}}" "$P/taxonomy/names.dmp"
        ln -s "${{NODES[0]}}" "$P/taxonomy/nodes.dmp"
        for F in "${{ACCESSIONS[@]}}"; do ln -s "$F" "$P/taxonomy/$(basename "$F")"; done
        if [ ${{#MAP[@]}} -gt 0 ]; then cp "${{MAP[0]}}" "$P/seqid2taxid.map"; fi

        for F in "${{FASTA[@]}}"; do
            kraken2-build \
                --add-to-library "$F" \
                --db "$P" \
                --threads $TASK_CPUS \
                $ARGS \
                || exit 1
        done
        tar -czhf "$OUT/$P.tar.gz" "$P"
    """,
)

kraken2_build_cmd = shell.create(
    name="kraken2_build",
    image=KRAKEN2_IMAGE,
    resources=BUILD_RESOURCES,
    inputs={"library": File, "clean": bool, "args": str},
    defaults={"clean": True, "args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        LIBRARY=({{inputs.library}})
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"

        mkdir kraken2-database
        tar -xzf "${{LIBRARY[0]}}" -C kraken2-database --strip-components 1
        kraken2-build \
            --build \
            $ARGS \
            --threads $TASK_CPUS \
            --db kraken2-database/ \
            || exit 1
        if [ {{inputs.clean}} = true ]; then
            kraken2-build --clean --db kraken2-database/ || exit 1
        else
            tar -czf "$OUT/kraken2-database.tar.gz" kraken2-database
        fi
        mkdir "$OUT/kraken2-database"
        find kraken2-database -maxdepth 1 -type f -exec cp {{}} "$OUT/kraken2-database/" \;
    """,
)

env = flyte.TaskEnvironment.from_task(
    "kraken2", kraken2_cmd.as_task(), kraken2_add_cmd.as_task(), kraken2_build_cmd.as_task()
)


@dataclass
class Kraken2Result:
    report: File
    classified_reads_fastq: list[File]  # with save_output_fastqs: one file, or mates 1 and 2
    unclassified_reads_fastq: list[File]  # with save_output_fastqs: one file, or mates 1 and 2
    classified_reads_assignment: File | None  # with save_reads_assignment
    results: Dir  # everything kraken2 wrote


async def kraken2(
    reads_1: File,
    db: Dir,
    reads_2: File | None = None,
    prefix: str = "sample",
    save_output_fastqs: bool = False,
    save_reads_assignment: bool = False,
    args: str = "",
) -> Kraken2Result:
    """Classify one sample's reads against ``db``; ``reads_2`` is None for single-end reads."""
    results = await kraken2_cmd(
        reads_1=reads_1,
        reads_2=reads_2,
        db=db,
        prefix=prefix,
        save_output_fastqs=save_output_fastqs,
        save_reads_assignment=save_reads_assignment,
        args=args,
    )

    async def pick(name: str) -> File:
        found = await results.get_file(name)
        if found is None:
            raise FileNotFoundError(f"kraken2 wrote no {name}")
        return found

    mates = ["_1", "_2"] if reads_2 is not None else [""]
    classified = []
    unclassified = []
    if save_output_fastqs:
        classified = [await pick(f"{prefix}.classified{mate}.fastq.gz") for mate in mates]
        unclassified = [await pick(f"{prefix}.unclassified{mate}.fastq.gz") for mate in mates]
    assignment = None
    if save_reads_assignment:
        assignment = await pick(f"{prefix}.kraken2.classifiedreads.txt")
    return Kraken2Result(
        report=await pick(f"{prefix}.kraken2.report.txt"),
        classified_reads_fastq=classified,
        unclassified_reads_fastq=unclassified,
        classified_reads_assignment=assignment,
        results=results,
    )


async def kraken2_add(
    fasta: File | list[File],
    taxonomy_names: File,
    taxonomy_nodes: File,
    accession2taxid: File | list[File] | None = None,
    seqid2taxid: File | None = None,
    prefix: str = "sample",
    args: str = "",
) -> File:
    """Add ``fasta`` to a new database's library: ``<prefix>.tar.gz`` of the database directory ``<prefix>/``."""
    fastas = fasta if isinstance(fasta, list) else [fasta]
    if accession2taxid is None:
        accessions = []
    else:
        accessions = accession2taxid if isinstance(accession2taxid, list) else [accession2taxid]
    results = await kraken2_add_cmd(
        fasta=fastas,
        taxonomy_names=taxonomy_names,
        taxonomy_nodes=taxonomy_nodes,
        accession2taxid=accessions,
        seqid2taxid=seqid2taxid,
        prefix=prefix,
        args=args,
    )
    archive = await results.get_file(f"{prefix}.tar.gz")
    if archive is None:
        raise FileNotFoundError(f"kraken2_add wrote no {prefix}.tar.gz")
    return archive


@dataclass
class Kraken2BuildResult:
    db: Dir  # the database's top-level files: hash.k2d, opts.k2d, taxo.k2d, ...
    archive: File | None  # without clean: the whole database directory, kraken2-database.tar.gz


async def kraken2_build(library: File, clean: bool = True, args: str = "") -> Kraken2BuildResult:
    """Build a Kraken2 database from :func:`kraken2_add`'s archive; ``clean`` removes library and taxonomy."""
    results = await kraken2_build_cmd(library=library, clean=clean, args=args)
    archive = await results.get_file("kraken2-database.tar.gz")
    if not clean and archive is None:
        raise FileNotFoundError("kraken2_build wrote no kraken2-database.tar.gz")
    return Kraken2BuildResult(
        db=Dir.from_existing_remote(results.path.rstrip("/") + "/kraken2-database"),
        archive=archive,
    )
