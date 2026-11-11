"""metacache — MetaCache taxonomic classification databases.

Exposes :func:`metacache_build`, which builds a MetaCache database
(``metacache build``) from reference sequences (FASTA or FASTQ, plain or
gzipped), linked into ``genomes/`` as upstream, optionally over a taxonomy
(``names.dmp`` and ``nodes.dmp``, plus optionally ``merged.dmp``, linked
into ``taxonomy/``) and a sequence-to-taxon map (``-taxpostmap``). The
database is ``<prefix>.meta`` and one or more ``<prefix>.cache<N>``
(``-parts`` in ``args`` splits it).

``metacache build`` is single-threaded, so the task gets one CPU.

:func:`metacache_query` maps one sample's reads (single- or paired-end,
``-pairfiles``) against a database Dir (found by its ``*.meta``): the
per-read mappings ``<prefix>.mapping.txt`` and, with ``abundances``, the
per-taxon counts ``<prefix>.abundances.txt``. ``-threads`` is the task's
CPUs; options go in ``args``.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

METACACHE_IMAGE = "quay.io/biocontainers/metacache:2.5.0--h077b44d_0"
# Upstream labels it process_single; real databases need far more memory.
DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory="8Gi")

TAXONOMY_NAMES = (["names.dmp", "nodes.dmp"], ["merged.dmp", "names.dmp", "nodes.dmp"])

metacache_build_cmd = shell.create(
    name="metacache_build",
    image=METACACHE_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"genomes": list[File], "taxonomy": list[File], "seq2taxid": File | None, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        GENOMES=({inputs.genomes})
        shopt -s nullglob; TAXONOMY=({inputs.taxonomy}); MAP=({inputs.seq2taxid}); shopt -u nullglob
        P={inputs.prefix}
        ARGS={inputs.args}
        OUT={outputs.results}
        W=$(mktemp -d); cd "$W"
        mkdir genomes
        for F in "${GENOMES[@]}"; do ln -s "$F" genomes/"$(basename "$F")"; done
        TAXONOMY_ARGS=()
        if [ ${#TAXONOMY[@]} -gt 0 ]; then
            mkdir taxonomy
            for F in "${TAXONOMY[@]}"; do ln -s "$F" taxonomy/"$(basename "$F")"; done
            TAXONOMY_ARGS=(-taxonomy taxonomy)
        fi
        MAP_ARGS=()
        if [ ${#MAP[@]} -gt 0 ]; then MAP_ARGS=(-taxpostmap "${MAP[0]}"); fi
        metacache \
            build \
            "$P" \
            genomes/ \
            "${TAXONOMY_ARGS[@]}" \
            "${MAP_ARGS[@]}" \
            $ARGS \
            || exit 1
        mv "$P".meta "$P".cache* "$OUT"/
    """,
)

# Upstream labels it process_medium.
QUERY_RESOURCES = flyte.Resources(cpu=4, memory="8Gi")

metacache_query_cmd = shell.create(
    name="metacache_query",
    image=METACACHE_IMAGE,
    resources=QUERY_RESOURCES,
    inputs={"reads_1": File, "reads_2": File | None, "db": Dir, "abundances": bool, "prefix": str, "args": str},
    defaults={"abundances": False, "args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        R1=({{inputs.reads_1}})
        shopt -s nullglob; R2=({{inputs.reads_2}}); shopt -u nullglob
        DB={{inputs.db}}
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        DB_META=$(find -L "$DB" -name "*.meta" | head -n 1)
        if [ -z "$DB_META" ]; then echo 'Database file "*.meta" not found!' >&2; exit 1; fi
        W=$(mktemp -d); cd "$W"
        INPUT=("${{R1[0]}}")
        if [ ${{#R2[@]}} -gt 0 ]; then INPUT=("${{R1[0]}}" "${{R2[0]}}" -pairfiles); fi
        ABUNDANCES=()
        if [ {{inputs.abundances}} = true ]; then ABUNDANCES=(-abundances "$P.abundances.txt"); fi
        metacache \
            query \
            "$DB_META" \
            "${{INPUT[@]}}" \
            "${{ABUNDANCES[@]}}" \
            -threads $TASK_CPUS \
            $ARGS \
            -out "$P.mapping.txt" \
            || exit 1
        shopt -s nullglob
        mv ./*mapping.txt ./*abundances.txt "$OUT"/
    """,
)

env = flyte.TaskEnvironment.from_task("metacache", metacache_build_cmd.as_task(), metacache_query_cmd.as_task())


@dataclass
class MetacacheBuildResult:
    meta: File  # <prefix>.meta
    caches: list[File]  # <prefix>.cache0, <prefix>.cache1, ...
    db: Dir  # all of the above


async def metacache_build(
    genomes: File | list[File],
    taxonomy: list[File] | None = None,
    seq2taxid: File | None = None,
    prefix: str = "sample",
    args: str = "",
) -> MetacacheBuildResult:
    """Build a MetaCache database of ``genomes``; ``taxonomy`` is ``names.dmp``, ``nodes.dmp`` (and ``merged.dmp``)."""
    files = genomes if isinstance(genomes, list) else [genomes]
    taxonomy = taxonomy or []
    if taxonomy and sorted(f.name for f in taxonomy) not in TAXONOMY_NAMES:
        raise ValueError(f"metacache taxonomy must be names.dmp and nodes.dmp (and merged.dmp), got {taxonomy}")
    db = await metacache_build_cmd(genomes=files, taxonomy=taxonomy, seq2taxid=seq2taxid, prefix=prefix, args=args)
    meta = await db.get_file(f"{prefix}.meta")
    if meta is None:
        raise FileNotFoundError(f"metacache build wrote no {prefix}.meta")
    caches = []
    while (cache := await db.get_file(f"{prefix}.cache{len(caches)}")) is not None:
        caches.append(cache)
    return MetacacheBuildResult(meta=meta, caches=caches, db=db)


@dataclass
class MetacacheQueryResult:
    mapping: File  # <prefix>.mapping.txt
    abundances: File | None  # <prefix>.abundances.txt, with abundances
    results: Dir  # both of the above


async def metacache_query(
    reads_1: File,
    db: Dir,
    reads_2: File | None = None,
    abundances: bool = False,
    prefix: str = "sample",
    args: str = "",
) -> MetacacheQueryResult:
    """Map one sample's reads (``reads_2`` None for single-end) against the MetaCache database ``db``."""
    results = await metacache_query_cmd(
        reads_1=reads_1, reads_2=reads_2, db=db, abundances=abundances, prefix=prefix, args=args
    )
    mapping = await results.get_file(f"{prefix}.mapping.txt")
    if mapping is None:
        raise FileNotFoundError(f"metacache query wrote no {prefix}.mapping.txt")
    return MetacacheQueryResult(
        mapping=mapping, abundances=await results.get_file(f"{prefix}.abundances.txt"), results=results
    )
