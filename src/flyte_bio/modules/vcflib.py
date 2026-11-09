"""vcflib — VCF manipulation with vcflib.

Exposes :func:`vcflib_vcffilter` and :func:`vcflib_vcffixup`.

:func:`vcflib_vcffilter` filters a VCF's records by INFO-field
expressions (``vcffilter -f``) and/or genotypes by genotype-field expressions
(``-g``), and bgzips the result to ``<prefix>.vcf.gz`` (upstream's default
prefix is ``<id>.filter``). Each expression may
contain spaces (``DP > 4``, ``QUAL > 20 & DP > 10``) and is passed as one
``-f``/``-g`` value; vcffilter keeps records that pass every ``-f``. Other
vcffilter options (``--invert``, ``--or``, ``-t PASS``…) go in ``args``, and
bgzip's in ``bgzip_args``. At least one expression is required, as upstream.
A bgzipped VCF needs its ``tbi``, which vcffilter loads from beside it.

:func:`vcflib_vcffixup` recomputes a VCF's AC, AN and AF from its genotypes
and pipes the result into ``bcftools view``, whose ``--output-type`` (in
``args2``) sets the output, ``<prefix>.<extension>`` (upstream's default
prefix is ``<id>.fixed``); with ``--write-index`` in ``args2`` the index lands
beside it. It runs on upstream's bcftools + vcflib image; :data:`env`
aggregates both images' environments.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

VCFLIB_IMAGE = "community.wave.seqera.io/library/vcflib:1.0.14--cc8ffb2c1a080797"
# Upstream labels vcffilter process_low.
DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory="2Gi")

# Expressions come one per line, so each may contain spaces.
vcflib_vcffilter_cmd = shell.create(
    name="vcflib_vcffilter",
    image=VCFLIB_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "vcf": File,
        "tbi": File | None,
        "info_filters": str,
        "genotype_filters": str,
        "prefix": str,
        "args": str,
        "bgzip_args": str,
    },
    defaults={"info_filters": "", "genotype_filters": "", "args": "", "bgzip_args": ""},
    outputs={"results": Dir},
    script=r"""
        VCF=({inputs.vcf})
        shopt -s nullglob; TBI=({inputs.tbi}); shopt -u nullglob
        INFO_FILTERS={inputs.info_filters}
        GENOTYPE_FILTERS={inputs.genotype_filters}
        P={inputs.prefix}
        ARGS={inputs.args}
        BGZIP_ARGS={inputs.bgzip_args}
        OUT={outputs.results}
        set -f
        FILTERS=()
        while IFS= read -r EXPRESSION; do
            if [ -n "$EXPRESSION" ]; then FILTERS+=(-f "$EXPRESSION"); fi
        done <<< "$INFO_FILTERS"
        while IFS= read -r EXPRESSION; do
            if [ -n "$EXPRESSION" ]; then FILTERS+=(-g "$EXPRESSION"); fi
        done <<< "$GENOTYPE_FILTERS"
        # vcffilter looks for a bgzipped VCF's index beside it.
        W=$(mktemp -d); cd "$W"
        NAME=$(basename "${VCF[0]}")
        ln -s "${VCF[0]}" "$NAME"
        if [ ${#TBI[@]} -gt 0 ]; then ln -s "${TBI[0]}" "$NAME.tbi"; fi
        set -o pipefail
        vcffilter \
            $ARGS \
            "${FILTERS[@]}" \
            "$NAME" \
            | bgzip -c $BGZIP_ARGS > "$OUT/$P.vcf.gz" || exit 1
    """,
)

# vcffixup's image adds bcftools.
VCFFIXUP_IMAGE = "community.wave.seqera.io/library/bcftools_vcflib:0b47030679d1eff1"

vcflib_vcffixup_cmd = shell.create(
    name="vcflib_vcffixup",
    image=VCFFIXUP_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"vcf": File, "tbi": File | None, "prefix": str, "extension": str, "args": str, "args2": str},
    defaults={"args": "", "args2": ""},
    outputs={"results": Dir},
    script=r"""
        set -o pipefail
        VCF=({inputs.vcf})
        shopt -s nullglob; TBI=({inputs.tbi}); shopt -u nullglob
        P={inputs.prefix}; EXT={inputs.extension}
        ARGS={inputs.args}; ARGS2={inputs.args2}
        OUT={outputs.results}
        W=$(mktemp -d); cd "$W"
        NAME=$(basename "${VCF[0]}")
        ln -s "${VCF[0]}" "$NAME"
        if [ ${#TBI[@]} -gt 0 ]; then ln -s "${TBI[0]}" "$NAME.tbi"; fi
        vcffixup $ARGS "$NAME" | bcftools view $ARGS2 -o "$P.$EXT" || exit 1
        for F in "$P.$EXT" "$P.$EXT".csi "$P.$EXT".tbi; do
            if [ -f "$F" ]; then mv "$F" "$OUT"/; fi
        done
    """,
)

vcffilter_env = flyte.TaskEnvironment.from_task("vcflib_vcffilter", vcflib_vcffilter_cmd.as_task())
vcffixup_env = flyte.TaskEnvironment.from_task("vcflib_vcffixup", vcflib_vcffixup_cmd.as_task())
env = flyte.TaskEnvironment(name="vcflib", depends_on=[vcffilter_env, vcffixup_env])


async def vcflib_vcffilter(
    vcf: File,
    info_filters: list[str] | None = None,
    genotype_filters: list[str] | None = None,
    tbi: File | None = None,
    prefix: str = "sample.filter",
    args: str = "",
    bgzip_args: str = "",
) -> File:
    """The records of ``vcf`` passing ``info_filters`` (and genotypes ``genotype_filters``): ``<prefix>.vcf.gz``."""
    info, genotype = list(info_filters or []), list(genotype_filters or [])
    if not info and not genotype:
        raise ValueError("vcflib_vcffilter requires either info_filters (-f) or genotype_filters (-g).")
    for expression in info + genotype:
        if not expression.strip() or "\n" in expression:
            raise ValueError(f"vcffilter expressions must be one non-empty line, got {expression!r}")
    if vcf.name == f"{prefix}.vcf.gz":
        raise ValueError("Input and output names are the same, set prefix to disambiguate!")
    results = await vcflib_vcffilter_cmd(
        vcf=vcf,
        tbi=tbi,
        info_filters="\n".join(info),
        genotype_filters="\n".join(genotype),
        prefix=prefix,
        args=args,
        bgzip_args=bgzip_args,
    )
    output = await results.get_file(f"{prefix}.vcf.gz")
    if output is None:
        raise FileNotFoundError(f"vcffilter wrote no {prefix}.vcf.gz")
    return output


def vcffixup_extension(args2: str) -> str:
    """The output's extension, read from bcftools view's ``--output-type`` in ``args2`` as upstream does."""
    for options, extension in (
        (("--output-type b", "-Ob"), "bcf.gz"),
        (("--output-type u", "-Ou"), "bcf"),
        (("--output-type z", "-Oz"), "vcf.gz"),
    ):
        if any(option in args2 for option in options):
            return extension
    return "vcf"


@dataclass
class VcffixupResult:
    vcf: File  # <prefix>.vcf, .vcf.gz, .bcf or .bcf.gz
    index: File | None  # <prefix>.<extension>.csi or .tbi, with --write-index in args2
    results: Dir


async def vcflib_vcffixup(
    vcf: File, tbi: File | None = None, prefix: str = "sample.fixed", args: str = "", args2: str = ""
) -> VcffixupResult:
    """``vcf`` with AC, AN and AF recomputed from its genotypes: ``<prefix>.<extension>``.

    ``args`` go to vcffixup, ``args2`` to the ``bcftools view`` that writes
    the output (its ``--output-type`` sets the format, plain VCF by default).
    """
    extension = vcffixup_extension(args2)
    if vcf.name == f"{prefix}.{extension}":
        raise ValueError("Input and output names are the same, set prefix to disambiguate!")
    results = await vcflib_vcffixup_cmd(vcf=vcf, tbi=tbi, prefix=prefix, extension=extension, args=args, args2=args2)
    output = await results.get_file(f"{prefix}.{extension}")
    if output is None:
        raise FileNotFoundError(f"vcffixup wrote no {prefix}.{extension}")
    index = await results.get_file(f"{prefix}.{extension}.csi") or await results.get_file(f"{prefix}.{extension}.tbi")
    return VcffixupResult(vcf=output, index=index, results=results)
