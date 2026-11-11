"""motus — taxonomic profiling of metagenomes with mOTUs marker genes.

The mOTUs database (``db_mOTU/``, about 3 GB compressed, CC-BY-4.0) is passed
to every tool as an archive (any tar compression) and extracted in the task;
the database is found by its ``db_mOTU_versions`` file.

- :func:`motus_downloaddb` downloads the database with mOTUs' own
  ``downloadDB.py`` (passed as a File, as upstream), which checks its md5 and
  writes ``db_mOTU_versions``. It returns ``db_mOTU.tar`` (uncompressed).
- :func:`motus_preplong` converts long reads into the short-read-like
  sequences ``motus profile`` takes (``motus prep_long``): ``<prefix>.gz``.
- :func:`motus_profile` profiles one sample (``motus profile``): single-end
  reads (``-s``), a pair (``-f``/``-r``), a BAM (``-i``) or an ``.mgc`` table
  (``-m``), picked as upstream from the first file's extension and the
  number of files. It writes the profile ``<prefix>.out`` (its header records
  the command line), the log ``<prefix>.log`` and, with ``-I``/``-M`` in
  ``args``, ``<prefix>.bam`` / ``<prefix>.mgc``. The sample is named
  ``<prefix>`` (``-n``).
- :func:`motus_merge` merges several profiles into one table (``motus
  merge``): ``<prefix>.txt``, or ``<prefix>.biom`` with ``-B`` in ``args``.

``-t`` is the task's CPUs; options go in ``args``.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

MOTUS_IMAGE = "quay.io/biocontainers/motus:3.1.0--pyhdfd78af_0"
# Upstream labels downloaddb process_low; the archive needs disk.
DOWNLOAD_RESOURCES = flyte.Resources(cpu=2, memory="4Gi", disk="40Gi")
# Upstream labels prep_long and profile process_medium, merge process_single; all read the database.
DEFAULT_RESOURCES = flyte.Resources(cpu=4, memory="16Gi", disk="40Gi")
MERGE_RESOURCES = flyte.Resources(cpu=1, memory="8Gi", disk="40Gi")

# Extract the database archive into the working directory; DB is the db_mOTU directory.
EXTRACT_DB = r"""
        mkdir database
        tar -xf "${DB_ARCHIVE[0]}" -C database || exit 1
        VERSIONS=$(find database -name db_mOTU_versions | head -n 1)
        if [ -z "$VERSIONS" ]; then echo "no mOTUs database (db_mOTU_versions) in the archive" >&2; exit 1; fi
        DB=$(dirname "$VERSIONS")
""".strip()

motus_downloaddb_cmd = shell.create(
    name="motus_downloaddb",
    image=MOTUS_IMAGE,
    resources=DOWNLOAD_RESOURCES,
    inputs={"script": File, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        SCRIPT=({{inputs.script}})
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        # As upstream: the script downloads next to itself, so copy it here.
        NAME=$(basename "${{SCRIPT[0]}}"); NAME="${{NAME%%.*}}_copy.py"
        cp "${{SCRIPT[0]}}" "$NAME"
        python "$NAME" \
            $ARGS \
            -t $TASK_CPUS \
            || exit 1
        tar -cf "$OUT/db_mOTU.tar" db_mOTU
    """,
)

motus_preplong_cmd = shell.create(
    name="motus_preplong",
    image=MOTUS_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"reads": File, "db": File, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        READS=({{inputs.reads}})
        DB_ARCHIVE=({{inputs.db}})
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        {EXTRACT_DB}
        motus \
            prep_long \
            $ARGS \
            -i "${{READS[0]}}" \
            -db "$DB" \
            -t $TASK_CPUS \
            -o "$P.gz" \
            2> "$P.log" \
            || {{ tail -n 20 "$P.log" >&2; exit 1; }}
        cat "$P.log" >&2
        mv "$P.gz" "$OUT"/
    """,
)

motus_profile_cmd = shell.create(
    name="motus_profile",
    image=MOTUS_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"reads": list[File], "reads_order": str, "db": File, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        READS=({{inputs.reads}})
        # Staged list inputs expand in name order; restore the order they were given in.
        read -ra READS_NAMES <<< {{inputs.reads_order}}; READS_DIR=$(dirname "${{READS[0]}}")
        READS=(); for N in "${{READS_NAMES[@]}}"; do READS+=("$READS_DIR/$N"); done
        DB_ARCHIVE=({{inputs.db}})
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        {EXTRACT_DB}
        case "${{READS[0]}}" in
            *.bam) INPUT=(-i "${{READS[@]}}") ;;
            *.mgc) INPUT=(-m "${{READS[@]}}") ;;
            *)
                if [ ${{#READS[@]}} -gt 1 ]; then
                    INPUT=(-f "${{READS[0]}}" -r "${{READS[1]}}")
                else
                    INPUT=(-s "${{READS[0]}}")
                fi
                ;;
        esac
        motus profile \
            $ARGS \
            "${{INPUT[@]}}" \
            -db "$DB" \
            -t $TASK_CPUS \
            -n "$P" \
            -o "$P.out" \
            2> "$P.log" \
            || {{ tail -n 20 "$P.log" >&2; exit 1; }}
        cat "$P.log" >&2
        shopt -s nullglob
        mv ./*.out ./*.bam ./*.mgc ./*.log "$OUT"/
    """,
)

motus_merge_cmd = shell.create(
    name="motus_merge",
    image=MOTUS_IMAGE,
    resources=MERGE_RESOURCES,
    inputs={"profiles": list[File], "profiles_order": str, "db": File, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        PROFILES=({{inputs.profiles}})
        # Staged list inputs expand in name order; restore the order they were given in.
        read -ra PROFILES_NAMES <<< {{inputs.profiles_order}}; PROFILES_DIR=$(dirname "${{PROFILES[0]}}")
        PROFILES=(); for N in "${{PROFILES_NAMES[@]}}"; do PROFILES+=("$PROFILES_DIR/$N"); done
        DB_ARCHIVE=({{inputs.db}})
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        {EXTRACT_DB}
        SUFFIX=txt
        if [[ " $ARGS " == *" -B "* ]]; then SUFFIX=biom; fi
        INPUT=$(IFS=,; echo "${{PROFILES[*]}}")
        motus \
            merge \
            -db "$DB" \
            -i "$INPUT" \
            $ARGS \
            -o "$P.$SUFFIX" \
            || exit 1
        mv "$P.$SUFFIX" "$OUT"/
    """,
)

env = flyte.TaskEnvironment.from_task(
    "motus",
    motus_downloaddb_cmd.as_task(),
    motus_preplong_cmd.as_task(),
    motus_profile_cmd.as_task(),
    motus_merge_cmd.as_task(),
)


async def motus_downloaddb(script: File, args: str = "") -> File:
    """Download the mOTUs database with mOTUs' ``downloadDB.py``: ``db_mOTU.tar``."""
    results = await motus_downloaddb_cmd(script=script, args=args)
    archive = await results.get_file("db_mOTU.tar")
    if archive is None:
        raise FileNotFoundError("motus downloadDB.py left no db_mOTU")
    return archive


async def motus_preplong(reads: File, db: File, prefix: str = "sample", args: str = "") -> File:
    """Long ``reads`` converted for ``motus profile``: ``<prefix>.gz``."""
    results = await motus_preplong_cmd(reads=reads, db=db, prefix=prefix, args=args)
    converted = await results.get_file(f"{prefix}.gz")
    if converted is None:
        raise FileNotFoundError(f"motus prep_long wrote no {prefix}.gz")
    return converted


@dataclass
class MotusProfileResult:
    out: File  # <prefix>.out, the profile
    log: File  # <prefix>.log
    bam: File | None  # <prefix>.bam, with -I
    mgc: File | None  # <prefix>.mgc, with -M
    results: Dir  # everything above


async def motus_profile(
    reads: File | list[File], db: File, prefix: str = "sample", args: str = ""
) -> MotusProfileResult:
    """Profile one sample: single-end reads, a pair, BAM files or an ``.mgc`` table."""
    files = reads if isinstance(reads, list) else [reads]
    results = await motus_profile_cmd(
        reads=files, reads_order=" ".join(f.name for f in files), db=db, prefix=prefix, args=args
    )
    out = await results.get_file(f"{prefix}.out")
    log = await results.get_file(f"{prefix}.log")
    if out is None or log is None:
        raise FileNotFoundError(f"motus profile wrote no {prefix}.out or {prefix}.log")
    return MotusProfileResult(
        out=out,
        log=log,
        bam=await results.get_file(f"{prefix}.bam"),
        mgc=await results.get_file(f"{prefix}.mgc"),
        results=results,
    )


async def motus_merge(profiles: File | list[File], db: File, prefix: str = "sample", args: str = "") -> File:
    """Several mOTUs profiles merged into one table: ``<prefix>.txt`` (``<prefix>.biom`` with ``-B``)."""
    files = profiles if isinstance(profiles, list) else [profiles]
    results = await motus_merge_cmd(
        profiles=files, profiles_order=" ".join(f.name for f in files), db=db, prefix=prefix, args=args
    )
    suffix = "biom" if "-B" in args.split() else "txt"
    merged = await results.get_file(f"{prefix}.{suffix}")
    if merged is None:
        raise FileNotFoundError(f"motus merge wrote no {prefix}.{suffix}")
    return merged
