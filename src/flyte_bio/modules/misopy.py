"""misopy — MISO (Mixture of Isoforms) annotation indexing.

Exposes :func:`misopy_index`, which indexes a GFF3 of alternative events or
isoforms with MISO's ``index_gff`` into the directory MISO's ``miso
--run`` reads: a pickle per gene in a directory per chromosome, the
``genes_to_filenames`` and ``compressed_ids_to_genes`` shelves and a copy of
the GFF (``genes.gff``).

As upstream, comment lines and records whose score column isn't numeric or
``.`` are dropped first. An index has a file per gene, too many to stage as a
directory, so it is returned as ``<prefix>.tar.gz``, holding the directory
``<prefix>/``. ``args`` takes index_gff's options.
"""

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

MISOPY_IMAGE = "community.wave.seqera.io/library/python_misopy:9d5a390611c447f5"
# Upstream: process_high; index_gff is single-threaded.
DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory="8Gi")

misopy_index_cmd = shell.create(
    name="misopy_index",
    image=MISOPY_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"gff3": File, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        GFF=({inputs.gff3})
        P={inputs.prefix}
        ARGS={inputs.args}
        OUT={outputs.results}
        W=$(mktemp -d); cd "$W"
        grep -v "^#" "${GFF[0]}" | \
            awk 'BEGIN{OFS="\t"} NF==9 && $6~/^[0-9.]+$/ || $6=="." {print}' > filtered.gff3
        index_gff --index filtered.gff3 "$P" $ARGS || exit 1
        # The image's Python, as it may lack tar.
        python -c 'import sys, tarfile; t = tarfile.open(sys.argv[1], "w:gz"); t.add(sys.argv[2]); t.close()' \
            "$OUT/$P.tar.gz" "$P" || exit 1
    """,
)

env = flyte.TaskEnvironment.from_task("misopy", misopy_index_cmd.as_task())


async def misopy_index(gff3: File, prefix: str, args: str = "") -> File:
    """MISO's index of ``gff3``: ``<prefix>.tar.gz``, holding the index directory ``<prefix>/``."""
    results = await misopy_index_cmd(gff3=gff3, prefix=prefix, args=args)
    found = await results.get_file(f"{prefix}.tar.gz")
    if found is None:
        raise FileNotFoundError(f"index_gff wrote no {prefix}.tar.gz")
    return found
