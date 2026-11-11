"""bedops — conversion of genomic annotations to BED with BEDOPS.

Exposes :func:`bedops_gtf2bed`, which converts a GTF to BED with ``gtf2bed``.
As upstream, each feature's BED name (column 4) is its ``exon_id`` attribute
(``--attribute-key=exon_id``), and the output is ``<prefix>.bed``, the
prefix defaulting to the GTF's name without its extension.

:func:`bedops_convert2bed` converts any format ``convert2bed`` reads (BAM,
GFF, GTF, GVF, PSL, RepeatMasker, SAM, VCF, WIG) to ``<prefix>.bed``; as
upstream, the input format is the file's extension (``-i``).
"""

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

BEDOPS_IMAGE = "community.wave.seqera.io/library/bedops:2.4.42--1f31d470625cb9a3"
# Upstream: process_low. gtf2bed sorts its output, in memory and temporary files.
DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory="4Gi")

bedops_gtf2bed_cmd = shell.create(
    name="bedops_gtf2bed",
    image=BEDOPS_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"gtf": File, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        GTF=({inputs.gtf})
        P={inputs.prefix}
        ARGS={inputs.args}
        OUT={outputs.results}
        W=$(mktemp -d); cd "$W"
        export TMPDIR="$W"

        # convert2bed (behind gtf2bed) exits 255 on errors, which would be retried forever.
        cat "${GTF[0]}" | gtf2bed $ARGS --attribute-key=exon_id > "$OUT/$P.bed" || exit 1
    """,
)

bedops_convert2bed_cmd = shell.create(
    name="bedops_convert2bed",
    image=BEDOPS_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"file": File, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        IN=({inputs.file})
        P={inputs.prefix}
        ARGS={inputs.args}
        OUT={outputs.results}
        W=$(mktemp -d); cd "$W"
        export TMPDIR="$W"
        FORMAT="${IN[0]##*.}"

        convert2bed $ARGS -i "$FORMAT" < "${IN[0]}" > "$OUT/$P.bed" || exit 1
    """,
)

env = flyte.TaskEnvironment.from_task("bedops", bedops_gtf2bed_cmd.as_task(), bedops_convert2bed_cmd.as_task())


async def bedops_gtf2bed(gtf: File, prefix: str = "", args: str = "") -> File:
    """``gtf`` as BED: ``<prefix>.bed`` (by default the GTF's name without its extension)."""
    if not prefix:
        prefix = gtf.name.rsplit(".", 1)[0] if "." in gtf.name else gtf.name
    results = await bedops_gtf2bed_cmd(gtf=gtf, prefix=prefix, args=args)
    bed = await results.get_file(f"{prefix}.bed")
    if bed is None:
        raise FileNotFoundError(f"gtf2bed wrote no {prefix}.bed")
    return bed


async def bedops_convert2bed(file: File, prefix: str, args: str = "") -> File:
    """``file`` as BED, read in the format its extension names: ``<prefix>.bed``."""
    results = await bedops_convert2bed_cmd(file=file, prefix=prefix, args=args)
    bed = await results.get_file(f"{prefix}.bed")
    if bed is None:
        raise FileNotFoundError(f"convert2bed wrote no {prefix}.bed")
    return bed
