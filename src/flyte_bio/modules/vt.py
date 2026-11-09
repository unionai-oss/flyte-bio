"""vt — variant tool set: decomposition and normalisation of VCFs.

Exposes:

- :func:`vt_decompose` — split multiallelic variants into biallelic records,
  optionally within intervals, as ``<prefix>.vcf.gz``.
- :func:`vt_normalize` — left-align and trim variants against the reference,
  optionally within intervals, as ``<prefix>.vcf.gz``; it also returns the
  reference's ``.fai`` (written by vt when none is given).

As upstream, ``args`` goes to vt and ``args2`` to the ``bgzip`` (``--threads``
of the task's CPUs) that compresses vt's output. Inputs are linked into the
working directory under their own names, so vt records them as upstream does.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

VT_IMAGE = "community.wave.seqera.io/library/htslib_vt:933758f4e7c3b9f1"
# Upstream's process_low.
DEFAULT_RESOURCES = flyte.Resources(cpu=2, memory="4Gi")

vt_decompose_cmd = shell.create(
    name="vt_decompose",
    image=VT_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"vcf": File, "intervals": File | None, "prefix": str, "args": str, "args2": str},
    defaults={"args": "", "args2": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        V=({{inputs.vcf}})
        shopt -s nullglob; BED=({{inputs.intervals}}); shopt -u nullglob
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        ARGS2={{inputs.args2}}
        W=$(mktemp -d); cd "$W"
        ln -s "${{V[0]}}" "$(basename "${{V[0]}}")"
        REGIONS=()
        if [ ${{#BED[@]}} -gt 0 ]; then
            BED_LINK=$(basename "${{BED[0]}}"); ln -s "${{BED[0]}}" "$BED_LINK"; REGIONS=(-i "$BED_LINK")
        fi
        vt decompose -o "$P.vcf" "${{REGIONS[@]}}" $ARGS "$(basename "${{V[0]}}")" || exit 1
        bgzip $ARGS2 --threads $TASK_CPUS "$P.vcf" || exit 1
        mv "$P.vcf.gz" {{outputs.results}}/
    """,
)

vt_normalize_cmd = shell.create(
    name="vt_normalize",
    image=VT_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "vcf": File,
        "tbi": File | None,
        "intervals": File | None,
        "fasta": File,
        "fai": File | None,
        "prefix": str,
        "args": str,
        "args2": str,
    },
    defaults={"args": "", "args2": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        V=({{inputs.vcf}}); FA=({{inputs.fasta}})
        shopt -s nullglob; TBI=({{inputs.tbi}}); BED=({{inputs.intervals}}); FAI=({{inputs.fai}}); shopt -u nullglob
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        ARGS2={{inputs.args2}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        VCF=$(basename "${{V[0]}}"); ln -s "${{V[0]}}" "$VCF"
        if [ ${{#TBI[@]}} -gt 0 ]; then ln -s "${{TBI[0]}}" "$VCF.tbi"; fi
        REF=$(basename "${{FA[0]}}"); ln -s "${{FA[0]}}" "$REF"
        if [ ${{#FAI[@]}} -gt 0 ]; then ln -s "${{FAI[0]}}" "$REF.fai"; fi
        REGIONS=()
        if [ ${{#BED[@]}} -gt 0 ]; then
            BED_LINK=$(basename "${{BED[0]}}"); ln -s "${{BED[0]}}" "$BED_LINK"; REGIONS=(-i "$BED_LINK")
        fi
        vt normalize -o "$P.vcf" -r "$REF" "${{REGIONS[@]}}" $ARGS "$VCF" || exit 1
        bgzip $ARGS2 --threads $TASK_CPUS "$P.vcf" || exit 1
        mv "$P.vcf.gz" "$OUT/"
        if [ -e "$REF.fai" ]; then cp -L "$REF.fai" "$OUT/"; fi
    """,
)

env = flyte.TaskEnvironment.from_task("vt", vt_decompose_cmd.as_task(), vt_normalize_cmd.as_task())


def check_prefix(vcf: File, prefix: str) -> None:
    if vcf.name in (f"{prefix}.vcf", f"{prefix}.vcf.gz"):
        raise ValueError("Input and output names are the same, set prefix to disambiguate!")


async def pick(results: Dir, name: str, tool: str) -> File:
    found = await results.get_file(name)
    if found is None:
        raise FileNotFoundError(f"vt {tool} wrote no {name}")
    return found


async def vt_decompose(vcf: File, prefix: str, intervals: File | None = None, args: str = "", args2: str = "") -> File:
    """``vcf`` with multiallelic variants split into biallelic records: ``<prefix>.vcf.gz``.

    Upstream's default prefix is ``<id>_decompose``.
    """
    check_prefix(vcf, prefix)
    results = await vt_decompose_cmd(vcf=vcf, intervals=intervals, prefix=prefix, args=args, args2=args2)
    return await pick(results, f"{prefix}.vcf.gz", "decompose")


@dataclass
class VtNormalizeResult:
    vcf: File  # <prefix>.vcf.gz
    fai: File | None  # <fasta name>.fai: the one given, or the one vt wrote


async def vt_normalize(
    vcf: File,
    fasta: File,
    prefix: str,
    tbi: File | None = None,
    intervals: File | None = None,
    fai: File | None = None,
    args: str = "",
    args2: str = "",
) -> VtNormalizeResult:
    """``vcf`` left-aligned and trimmed against ``fasta``: ``<prefix>.vcf.gz``.

    ``intervals`` (vt's ``-i``) needs the VCF's ``tbi``. Upstream's default
    prefix is ``<id>_normalize``.
    """
    check_prefix(vcf, prefix)
    results = await vt_normalize_cmd(
        vcf=vcf, tbi=tbi, intervals=intervals, fasta=fasta, fai=fai, prefix=prefix, args=args, args2=args2
    )
    return VtNormalizeResult(
        vcf=await pick(results, f"{prefix}.vcf.gz", "normalize"), fai=await results.get_file(f"{fasta.name}.fai")
    )
