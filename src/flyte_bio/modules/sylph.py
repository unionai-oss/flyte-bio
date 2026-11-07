"""sylph — metagenomic profiling by k-mer containment ANI.

Exposes :func:`sylph_profile`, which profiles one sample's reads against a
sylph database (``.syldb``) or genome FASTA files (sketched on the fly) and
writes upstream's ``<prefix>.tsv``. Reads and databases are linked into the
working directory under their own names, so the table's ``Sample_file`` and
``Genome_file`` columns name them as upstream's do (sylph-tax then names its
output after the sample file).
"""

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

SYLPH_IMAGE = "quay.io/biocontainers/sylph:0.9.0--ha6fb395_0"
# Upstream labels it process_high; real databases (GTDB) want about 16 GiB.
DEFAULT_RESOURCES = flyte.Resources(cpu=4, memory="16Gi")

sylph_profile_cmd = shell.create(
    name="sylph_profile",
    image=SYLPH_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"reads_1": File, "reads_2": File | None, "database": list[File], "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        R1=({{inputs.reads_1}}); DB=({{inputs.database}})
        shopt -s nullglob; R2=({{inputs.reads_2}}); shopt -u nullglob
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        DATABASE=()
        for F in "${{DB[@]}}"; do DATABASE+=("$(basename "$F")"); ln -s "$F" "$(basename "$F")"; done
        R1_NAME=$(basename "${{R1[0]}}"); ln -s "${{R1[0]}}" "$R1_NAME"
        if [ ${{#R2[@]}} -gt 0 ]; then
            R2_NAME=$(basename "${{R2[0]}}"); ln -s "${{R2[0]}}" "$R2_NAME"
            READS=(-1 "$R1_NAME" -2 "$R2_NAME")
        else
            READS=(-r "$R1_NAME")
        fi
        sylph profile \
            -t $TASK_CPUS \
            $ARGS \
            "${{DATABASE[@]}}" \
            "${{READS[@]}}" \
            -o "$P".tsv \
            || exit 1
        mv "$P".tsv "$OUT"/
    """,
)

env = flyte.TaskEnvironment.from_task("sylph", sylph_profile_cmd.as_task())


async def sylph_profile(
    reads_1: File,
    database: File | list[File],
    reads_2: File | None = None,
    prefix: str = "sample",
    args: str = "",
) -> File:
    """Profile one sample's reads (``reads_2`` None for single-end) against ``database``: ``<prefix>.tsv``."""
    databases = database if isinstance(database, list) else [database]
    results = await sylph_profile_cmd(reads_1=reads_1, reads_2=reads_2, database=databases, prefix=prefix, args=args)
    profile = await results.get_file(f"{prefix}.tsv")
    if profile is None:
        raise FileNotFoundError(f"sylph profile wrote no {prefix}.tsv")
    return profile
