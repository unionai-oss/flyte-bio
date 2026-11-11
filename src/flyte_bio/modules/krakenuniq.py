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

:func:`krakenuniq_preloadedkrakenuniq` classifies a batch of samples
against one database, loaded into memory once (``krakenuniq --preload``,
with ``args``) and then used for each sample (with ``args2``). With
``--preload-size`` in ``args2`` the database is instead read in chunks of
that size and not preloaded. Samples are single files or, with ``paired``,
pairs of files, each named by its entry in ``prefixes``; ``sequence_type``
(``fasta`` or ``fastq``) names the read outputs. Outputs keep upstream's
names:

- with ``report_file``: ``<prefix>.krakenuniq.report.txt``;
- with ``save_output``: the per-read classifications,
  ``<prefix>.krakenuniq.classified.txt``;
- with ``save_output_reads``: the classified and unclassified reads,
  ``<prefix>.classified.<sequence_type>.gz`` and
  ``<prefix>.unclassified.<sequence_type>.gz`` (``<prefix>.merged.*`` for
  pairs).

The database is a Dir (as :func:`krakenuniq_build` returns) or a ``.tar.gz``
of one, found by its ``database.kdb``. ``--threads`` is the task's CPUs.
"""

from dataclasses import dataclass

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

# Upstream labels it process_high; real databases need far more memory.
CLASSIFY_RESOURCES = flyte.Resources(cpu=4, memory="8Gi")

krakenuniq_preloadedkrakenuniq_cmd = shell.create(
    name="krakenuniq_preloadedkrakenuniq",
    image=KRAKENUNIQ_IMAGE,
    resources=CLASSIFY_RESOURCES,
    inputs={
        "sequences": list[File],
        "sequences_order": str,
        "prefixes": str,
        "paired": bool,
        "sequence_type": str,
        "db": Dir | None,
        "db_archive": File | None,
        "save_output_reads": bool,
        "report_file": bool,
        "save_output": bool,
        "args": str,
        "args2": str,
    },
    defaults={"args": "", "args2": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        SEQUENCES=({{inputs.sequences}})
        # Staged list inputs expand in name order; restore the order they were given in.
        read -ra SEQUENCES_NAMES <<< {{inputs.sequences_order}}; SEQUENCES_DIR=$(dirname "${{SEQUENCES[0]}}")
        SEQUENCES=(); for N in "${{SEQUENCES_NAMES[@]}}"; do SEQUENCES+=("$SEQUENCES_DIR/$N"); done
        shopt -s nullglob; ARCHIVE=({{inputs.db_archive}}); shopt -u nullglob
        PREFIX_LIST={{inputs.prefixes}}
        read -ra PREFIXES <<< "$PREFIX_LIST"
        TYPE={{inputs.sequence_type}}
        DB={{inputs.db}}
        ARGS={{inputs.args}}
        ARGS2={{inputs.args2}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        # An unset Dir input renders as a path that doesn't exist (or is empty).
        if [ ${{#ARCHIVE[@]}} -gt 0 ]; then
            mkdir database
            tar -xzf "${{ARCHIVE[0]}}" -C database || exit 1
            DB=database
        fi
        KDB=$(find -L "$DB" -name database.kdb | head -n 1)
        if [ -z "$KDB" ]; then echo "no KrakenUniq database (database.kdb) in the database" >&2; exit 1; fi
        DB=$(dirname "$KDB")

        if [[ " $ARGS2 " != *" --preload-size "* ]]; then
            krakenuniq $ARGS --db "$DB" --preload --threads $TASK_CPUS || exit 1
        fi
        STEP=1; PAIRED=()
        if [ {{inputs.paired}} = true ]; then STEP=2; PAIRED=(--paired); fi
        I=0
        for PREFIX in "${{PREFIXES[@]}}"; do
            INPUT=("${{SEQUENCES[@]:$I:$STEP}}")
            I=$(( I + STEP ))
            if [ {{inputs.paired}} = true ]; then
                CLASSIFIED="$PREFIX.merged.classified.$TYPE"; UNCLASSIFIED="$PREFIX.merged.unclassified.$TYPE"
            else
                CLASSIFIED="$PREFIX.classified.$TYPE"; UNCLASSIFIED="$PREFIX.unclassified.$TYPE"
            fi
            EXTRA=()
            if [ {{inputs.report_file}} = true ]; then EXTRA+=(--report-file "$PREFIX.krakenuniq.report.txt"); fi
            if [ {{inputs.save_output}} = true ]; then EXTRA+=(--output "$PREFIX.krakenuniq.classified.txt"); fi
            if [ {{inputs.save_output_reads}} = true ]; then
                EXTRA+=(--unclassified-out "$UNCLASSIFIED" --classified-out "$CLASSIFIED")
            fi
            krakenuniq \
                --db "$DB" \
                --threads $TASK_CPUS \
                "${{EXTRA[@]}}" \
                "${{PAIRED[@]}}" \
                $ARGS2 \
                "${{INPUT[@]}}" \
                || exit 1
        done
        shopt -s nullglob
        for F in ./*."$TYPE"; do gzip --no-name "$F"; done
        FOUND=(./*.classified."$TYPE".gz ./*.unclassified."$TYPE".gz)
        FOUND+=(./*.krakenuniq.classified.txt ./*.krakenuniq.report.txt)
        if [ ${{#FOUND[@]}} -gt 0 ]; then mv "${{FOUND[@]}}" "$OUT"/; fi
    """,
)

env = flyte.TaskEnvironment.from_task(
    "krakenuniq", krakenuniq_build_cmd.as_task(), krakenuniq_preloadedkrakenuniq_cmd.as_task()
)


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


@dataclass
class KrakenuniqResult:
    classified_reads: list[File]  # with save_output_reads, one per sample
    unclassified_reads: list[File]  # with save_output_reads, one per sample
    classified_assignment: list[File]  # with save_output: <prefix>.krakenuniq.classified.txt per sample
    report: list[File]  # with report_file: <prefix>.krakenuniq.report.txt per sample
    results: Dir  # everything above


async def krakenuniq_preloadedkrakenuniq(
    sequences: list[File],
    prefixes: list[str],
    db: Dir | File,
    sequence_type: str = "fastq",
    paired: bool = False,
    save_output_reads: bool = False,
    report_file: bool = False,
    save_output: bool = False,
    args: str = "",
    args2: str = "",
) -> KrakenuniqResult:
    """Classify a batch of samples against ``db``: ``sequences`` named by ``prefixes`` (in pairs with ``paired``)."""
    if sequence_type not in ("fasta", "fastq"):
        raise ValueError(f"sequence_type must be fasta or fastq, got {sequence_type!r}")
    if len(sequences) != len(prefixes) * (2 if paired else 1):
        raise ValueError(f"{len(sequences)} sequence files for {len(prefixes)} prefixes (paired={paired})")
    if any(not p or any(c.isspace() for c in p) for p in prefixes):
        raise ValueError(f"prefixes must be non-empty and contain no whitespace: {prefixes}")
    results = await krakenuniq_preloadedkrakenuniq_cmd(
        sequences=sequences,
        sequences_order=" ".join(f.name for f in sequences),
        prefixes=" ".join(prefixes),
        paired=paired,
        sequence_type=sequence_type,
        db=db if isinstance(db, Dir) else None,
        db_archive=db if isinstance(db, File) else None,
        save_output_reads=save_output_reads,
        report_file=report_file,
        save_output=save_output,
        args=args,
        args2=args2,
    )

    async def pick(suffix: str) -> list[File]:
        stem = ".merged" if paired and suffix.endswith(".gz") else ""
        found = [await results.get_file(f"{p}{stem}{suffix}") for p in prefixes]
        return [f for f in found if f is not None]

    return KrakenuniqResult(
        classified_reads=await pick(f".classified.{sequence_type}.gz"),
        unclassified_reads=await pick(f".unclassified.{sequence_type}.gz"),
        classified_assignment=await pick(".krakenuniq.classified.txt"),
        report=await pick(".krakenuniq.report.txt"),
        results=results,
    )
