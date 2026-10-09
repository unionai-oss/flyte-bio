"""vcftools — VCF summaries.

Exposes :func:`vcftools`, which runs vcftools on a VCF (``.vcf``,
``.vcf.gz`` or ``.bcf``, picked from its name as upstream does) with the
given options, writing ``<prefix>.<report>`` files (e.g. ``.TsTv.count``,
``.FILTER.summary``, ``.frq``) into one output Dir.
"""

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

VCFTOOLS_IMAGE = "quay.io/biocontainers/vcftools:0.1.17--pl5321h077b44d_0"

DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory="2Gi")

BED_FLAGS = ("--bed", "--exclude-bed", "--hapcount", "--positions", "--exclude-positions")

vcftools_cmd = shell.create(
    name="vcftools",
    image=VCFTOOLS_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"variant_file": File, "bed": File | None, "bed_flag": str, "prefix": str, "args": str},
    defaults={"bed_flag": "", "args": ""},
    outputs={"results": Dir},
    script=r"""
        V=({inputs.variant_file})
        shopt -s nullglob; B=({inputs.bed}); shopt -u nullglob
        ARGS={inputs.args}
        FLAG={inputs.bed_flag}
        W=$(mktemp -d); cd "$W"
        NAME=$(basename "${V[0]}"); ln -s "${V[0]}" "$NAME"
        case "$NAME" in
            *.vcf) INPUT=(--vcf "$NAME") ;;
            *.vcf.gz) INPUT=(--gzvcf "$NAME") ;;
            *.bcf) INPUT=(--bcf "$NAME") ;;
            *) INPUT=() ;;
        esac
        BED=(); if [ -n "$FLAG" ] && [ ${#B[@]} -gt 0 ]; then BED=("$FLAG" "${B[0]}"); fi
        vcftools "${INPUT[@]}" --out {outputs.results}/{inputs.prefix} $ARGS "${BED[@]}"
    """,
)


env = flyte.TaskEnvironment.from_task(
    "vcftools",
    vcftools_cmd.as_task(),
)


async def vcftools(variant_file: File, prefix: str, args: str, bed: File | None = None) -> Dir:
    """vcftools' reports for ``variant_file``, as ``<prefix>.*``.

    As upstream, a ``bed`` is passed with whichever of ``--bed``,
    ``--exclude-bed``, ``--hapcount``, ``--positions`` or
    ``--exclude-positions`` appears (bare) in ``args``.
    """
    tokens = args.split()
    flag = next((t for t in tokens if t in BED_FLAGS), "")
    rest = " ".join(t for t in tokens if t not in BED_FLAGS)
    return await vcftools_cmd(variant_file=variant_file, bed=bed, bed_flag=flag, prefix=prefix, args=rest)
