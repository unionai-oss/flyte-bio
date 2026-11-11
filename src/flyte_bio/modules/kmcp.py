"""kmcp — k-mer-based metagenomic profiling databases.

Exposes the two steps of building a KMCP database, as upstream:

- :func:`kmcp_compute` computes k-mers of reference sequences (``kmcp
  compute``): FASTA/FASTQ files, linked into ``genomes/``, or a directory of
  them used as ``genomes/``. It writes one ``.unik`` file per reference and a
  summary, ``_info.txt``, into ``<prefix>/``; as a database holds many
  references, that directory comes back as ``<prefix>.tar.gz``, beside
  ``_info.txt``.
- :func:`kmcp_index` indexes that archive (``kmcp index``) into
  ``<prefix>/`` (``R001/__db.yml``, ``__name_mapping.tsv``,
  ``_block*.uniki``, ...), returned as ``<prefix>.tar.gz`` with the log
  ``<prefix>.log``. Given taxonomy dump files and/or a sequence-to-taxon map,
  it copies them into ``<prefix>/taxdmp/`` for ``kmcp profile``. As
  ``_info.txt`` names the ``.unik`` files by the compute directory's name,
  the archive is extracted under that name, so the two prefixes must differ.

Then, against :func:`kmcp_index`'s archive:

- :func:`kmcp_search` searches one sample's reads (single- or paired-end)
  against the database (``kmcp search``): ``<prefix>.gz``.
- :func:`kmcp_profile` profiles search results (``kmcp profile``) at
  ``level`` ``strain``, ``assembly`` or ``species``: ``<prefix>.profile``.
  Species profiles need the taxonomy dump and sequence-to-taxon map that
  :func:`kmcp_index` copies into the database (``-X``, ``-T``); without them
  the task fails, as upstream.

``--threads`` (``-j`` for profile) is the task's CPUs; options go in
``args``.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

KMCP_IMAGE = "quay.io/biocontainers/kmcp:0.9.4--h9ee0642_0"
# Upstream labels both process_medium.
DEFAULT_RESOURCES = flyte.Resources(cpu=4, memory="8Gi")

kmcp_compute_cmd = shell.create(
    name="kmcp_compute",
    image=KMCP_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"sequences": list[File], "sequences_dir": Dir | None, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        shopt -s nullglob; SEQUENCES=({{inputs.sequences}}); shopt -u nullglob
        SEQUENCES_DIR={{inputs.sequences_dir}}
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        # An unset Dir input renders as a path that doesn't exist (or is empty).
        if [ -d "$SEQUENCES_DIR" ] && [ -n "$(ls -A "$SEQUENCES_DIR")" ]; then
            ln -s "$SEQUENCES_DIR" genomes
        else
            mkdir genomes
            for F in "${{SEQUENCES[@]}}"; do ln -s "$F" genomes/"$(basename "$F")"; done
        fi
        kmcp \
            compute \
            $ARGS \
            --threads $TASK_CPUS \
            --out-dir "$P"/ \
            --in-dir genomes/ \
            || exit 1
        cp "$P/_info.txt" "$OUT"/
        tar -czf "$OUT/$P.tar.gz" "$P"
    """,
)

kmcp_index_cmd = shell.create(
    name="kmcp_index",
    image=KMCP_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"compute": File, "taxdmp": list[File], "seq2taxidmap": list[File], "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        COMPUTE=({{inputs.compute}})
        shopt -s nullglob; TAXDMP=({{inputs.taxdmp}} {{inputs.seq2taxidmap}}); shopt -u nullglob
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        # _info.txt names the .unik files by the compute directory's own name, so keep it.
        tar -xzf "${{COMPUTE[0]}}"
        COMPUTE_DIR=$(tar -tzf "${{COMPUTE[0]}}" | awk -F/ 'NR == 1 {{ print $1 }}')
        kmcp \
            index \
            --in-dir "$COMPUTE_DIR" \
            $ARGS \
            --threads $TASK_CPUS \
            --log "$P.log" \
            --out-dir "$P" \
            || exit 1
        # Optionally copy over taxonomy files if provided for downstream module
        if [ ${{#TAXDMP[@]}} -gt 0 ]; then
            mkdir "$P/taxdmp"
            for F in "${{TAXDMP[@]}}"; do cp "$F" "$P/taxdmp/$(basename "$F")"; done
        fi
        mv "$P.log" "$OUT"/
        tar -czf "$OUT/$P.tar.gz" "$P"
    """,
)

# Upstream labels search process_medium and profile process_low.
PROFILE_RESOURCES = flyte.Resources(cpu=2, memory="4Gi")

# Extract the database archive into the working directory; DB_DIR is its top directory.
EXTRACT_DB = r"""
        tar -xzf "${DB[0]}" || exit 1
        DB_DIR=$(tar -tzf "${DB[0]}" | awk -F/ 'NR == 1 { print $1 }')
""".strip()

kmcp_search_cmd = shell.create(
    name="kmcp_search",
    image=KMCP_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"reads_1": File, "reads_2": File | None, "db": File, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        R1=({{inputs.reads_1}})
        shopt -s nullglob; R2=({{inputs.reads_2}}); shopt -u nullglob
        DB=({{inputs.db}})
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        {EXTRACT_DB}
        if [ ${{#R2[@]}} -gt 0 ]; then INPUT=(-1 "${{R1[0]}}" -2 "${{R2[0]}}"); else INPUT=("${{R1[0]}}"); fi
        kmcp \
            search \
            $ARGS \
            --threads $TASK_CPUS \
            --db-dir "$DB_DIR" \
            "${{INPUT[@]}}" \
            --out-file "$P.gz" \
            || exit 1
        mv "$P.gz" "$OUT"/
    """,
)

kmcp_profile_cmd = shell.create(
    name="kmcp_profile",
    image=KMCP_IMAGE,
    resources=PROFILE_RESOURCES,
    inputs={"search_results": File, "db": File, "level": str, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        RESULTS=({{inputs.search_results}})
        DB=({{inputs.db}})
        LEVEL={{inputs.level}}
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        {EXTRACT_DB}
        ln -s "${{RESULTS[0]}}" "$(basename "${{RESULTS[0]}}")"
        TAXONOMY=()
        if [ "$LEVEL" = species ]; then
            TAXID_MAP=$(find -L . -name '*map' -type f)
            TAXDMP=$(find -L . -name 'nodes.dmp' -type f | sed 's#nodes.dmp$##g')
            ERROR="[KMCP_PROFILE] ERROR: --level 'species' requires"
            MISSING="to be provided in the database directory. These files are not detected."
            if [ -z "$TAXID_MAP" ]; then echo "$ERROR a reference to taxid map file $MISSING" >> errors.txt; fi
            if [ -z "$TAXDMP" ]; then echo "$ERROR taxdump files $MISSING" >> errors.txt; fi
            if [ -s errors.txt ]; then
                cat errors.txt >&2
                # stderr goes to a file; the termination log carries the errors into the task's failure.
                {{ cat errors.txt > /dev/termination-log; }} 2> /dev/null
                exit 1
            fi
            TAXONOMY=(-X "$TAXDMP" -T "$TAXID_MAP")
        fi
        kmcp \
            profile \
            --level "$LEVEL" \
            -j $TASK_CPUS \
            -o "$P.profile" \
            $ARGS \
            "${{TAXONOMY[@]}}" \
            "$(basename "${{RESULTS[0]}}")" \
            || exit 1
        mv "$P.profile" "$OUT"/
    """,
)

env = flyte.TaskEnvironment.from_task(
    "kmcp",
    kmcp_compute_cmd.as_task(),
    kmcp_index_cmd.as_task(),
    kmcp_search_cmd.as_task(),
    kmcp_profile_cmd.as_task(),
)


@dataclass
class KmcpComputeResult:
    archive: File  # <prefix>.tar.gz: the directory <prefix>/ of .unik files and _info.txt
    info: File  # _info.txt


async def kmcp_compute(sequences: File | list[File] | Dir, prefix: str = "sample", args: str = "") -> KmcpComputeResult:
    """Compute k-mers of ``sequences`` (files, or a directory of them)."""
    if isinstance(sequences, Dir):
        files, directory = [], sequences
    else:
        files, directory = (sequences if isinstance(sequences, list) else [sequences]), None
    results = await kmcp_compute_cmd(sequences=files, sequences_dir=directory, prefix=prefix, args=args)
    archive = await results.get_file(f"{prefix}.tar.gz")
    info = await results.get_file("_info.txt")
    if archive is None or info is None:
        raise FileNotFoundError(f"kmcp compute wrote no {prefix}.tar.gz or _info.txt")
    return KmcpComputeResult(archive=archive, info=info)


@dataclass
class KmcpIndexResult:
    archive: File  # <prefix>.tar.gz: the database directory <prefix>/
    log: File  # <prefix>.log


async def kmcp_index(
    compute: File,
    taxdmp: list[File] | None = None,
    seq2taxidmap: File | None = None,
    prefix: str = "sample",
    args: str = "",
) -> KmcpIndexResult:
    """Index :func:`kmcp_compute`'s archive; taxonomy files are copied into ``<prefix>/taxdmp/``."""
    if compute.name == f"{prefix}.tar.gz":
        raise ValueError(f"kmcp index's prefix {prefix!r} must differ from kmcp compute's")
    results = await kmcp_index_cmd(
        compute=compute,
        taxdmp=taxdmp or [],
        seq2taxidmap=[seq2taxidmap] if seq2taxidmap is not None else [],
        prefix=prefix,
        args=args,
    )
    archive = await results.get_file(f"{prefix}.tar.gz")
    log = await results.get_file(f"{prefix}.log")
    if archive is None or log is None:
        raise FileNotFoundError(f"kmcp index wrote no {prefix}.tar.gz or {prefix}.log")
    return KmcpIndexResult(archive=archive, log=log)


async def kmcp_search(
    reads_1: File, db: File, reads_2: File | None = None, prefix: str = "sample", args: str = ""
) -> File:
    """Search one sample's reads (``reads_2`` None for single-end) against :func:`kmcp_index`'s archive.

    Returns the search results, ``<prefix>.gz``.
    """
    results = await kmcp_search_cmd(reads_1=reads_1, reads_2=reads_2, db=db, prefix=prefix, args=args)
    found = await results.get_file(f"{prefix}.gz")
    if found is None:
        raise FileNotFoundError(f"kmcp search wrote no {prefix}.gz")
    return found


PROFILE_LEVELS = ("species", "strain", "assembly")


async def kmcp_profile(search_results: File, db: File, level: str, prefix: str = "sample", args: str = "") -> File:
    """Profile :func:`kmcp_search`'s results at ``level`` (species, strain or assembly): ``<prefix>.profile``."""
    if level not in PROFILE_LEVELS:
        raise ValueError("[KMCP_PROFILE] ERROR: --level must be one of 'species', 'strain' or 'assembly'")
    results = await kmcp_profile_cmd(search_results=search_results, db=db, level=level, prefix=prefix, args=args)
    profile = await results.get_file(f"{prefix}.profile")
    if profile is None:
        raise FileNotFoundError(f"kmcp profile wrote no {prefix}.profile")
    return profile
