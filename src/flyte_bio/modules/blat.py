"""blat — BLAT alignment of sequences to a reference.

Exposes :func:`blat`, aligning a query FASTA (gzipped or not) to a subject
FASTA, writing ``<prefix>.psl``. As upstream, a gzipped query is decompressed
to ``<prefix>.fasta`` first.
"""

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

BLAT_IMAGE = "community.wave.seqera.io/library/ucsc-blat:472--3b87dacc9a173a42"
DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory="6Gi")

blat_cmd = shell.create(
    name="blat",
    image=BLAT_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"query": File, "subject": File, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        QUERY=({inputs.query}); SUBJECT=({inputs.subject})
        P={inputs.prefix}
        ARGS={inputs.args}
        OUT={outputs.results}
        W=$(mktemp -d); cd "$W"

        IN="${QUERY[0]}"
        if [[ $IN == *.gz ]]; then
            gunzip -cdf "$IN" > "$P.fasta"
            IN="$P.fasta"
        fi
        blat $ARGS "${SUBJECT[0]}" "$IN" "$OUT/$P.psl"
    """,
)

env = flyte.TaskEnvironment.from_task("blat", blat_cmd.as_task())


async def blat(query: File, subject: File, prefix: str = "sample", args: str = "") -> File:
    """Align ``query`` to ``subject``: ``<prefix>.psl``."""
    results = await blat_cmd(query=query, subject=subject, prefix=prefix, args=args)
    found = await results.get_file(f"{prefix}.psl")
    if found is None:
        raise FileNotFoundError(f"blat wrote no {prefix}.psl")
    return found
