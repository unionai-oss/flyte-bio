"""regtools — junction extraction from RNA-seq alignments.

Exposes :func:`regtools_junctionsextract`, which runs ``regtools junctions
extract`` on an indexed BAM and writes the splice junctions it supports, in
BED12 form, to ``<prefix>.junc``. ``strand_specificity`` is ``XS`` (the
default: the aligner's XS tag), ``RF`` (first-strand) or ``FR``
(second-strand). ``args`` takes further options (``-a``, ``-m``, ``-M``, ...).
"""

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

REGTOOLS_IMAGE = "community.wave.seqera.io/library/regtools:1.0.0--461ddf16709a70cf"
# Upstream: process_low.
DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory="4Gi")

regtools_junctionsextract_cmd = shell.create(
    name="regtools_junctionsextract",
    image=REGTOOLS_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"bam": File, "bai": File, "strand_specificity": str, "prefix": str, "args": str},
    defaults={"strand_specificity": "XS", "args": ""},
    outputs={"results": Dir},
    script=r"""
        BAM=({inputs.bam}); BAI=({inputs.bai})
        STRAND={inputs.strand_specificity}
        P={inputs.prefix}
        ARGS={inputs.args}
        OUT={outputs.results}
        W=$(mktemp -d); cd "$W"
        BAM_NAME=$(basename "${BAM[0]}")
        ln -s "${BAM[0]}" "$BAM_NAME"; ln -s "${BAI[0]}" "$BAM_NAME.bai"
        regtools junctions extract $ARGS -s "$STRAND" -o "$OUT/$P.junc" "$BAM_NAME" || exit 1
    """,
)

env = flyte.TaskEnvironment.from_task("regtools", regtools_junctionsextract_cmd.as_task())


async def regtools_junctionsextract(
    bam: File, bai: File, prefix: str, strand_specificity: str = "", args: str = ""
) -> File:
    """The splice junctions of ``bam`` (indexed by ``bai``): ``<prefix>.junc``.

    An empty ``strand_specificity`` means ``XS``, as upstream.
    """
    results = await regtools_junctionsextract_cmd(
        bam=bam, bai=bai, strand_specificity=strand_specificity or "XS", prefix=prefix, args=args
    )
    found = await results.get_file(f"{prefix}.junc")
    if found is None:
        raise FileNotFoundError(f"regtools wrote no {prefix}.junc")
    return found
