"""svtk — svtk, structural variant toolkit.

Exposes:

- :func:`svtk_standardize` — convert one SV caller's VCF (``caller``:
  ``delly``, ``lumpy``, ``manta``, ``wham``, ``melt``, …) into svtk's
  standard VCF, ``<prefix>.vcf.gz``, optionally restricted to the contigs of
  a FASTA index (``--contigs``).
- :func:`svtk_vcf2bed` — a VCF's SVs as BED, ``<prefix>.bed``.
"""

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

# Upstream's main.nf names the bare biocontainer (Docker Hub); this is its quay.io mirror.
SVTK_IMAGE = "quay.io/biocontainers/svtk:0.0.20190615--py37h73a75cf_2"
# Upstream's process_low and process_single, scaled to the tool: svtk is single-threaded.
DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory="4Gi")

svtk_standardize_cmd = shell.create(
    name="svtk_standardize",
    image=SVTK_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"vcf": File, "fai": File | None, "caller": str, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        VCF=({inputs.vcf})
        shopt -s nullglob; FAI=({inputs.fai}); shopt -u nullglob
        ARGS={inputs.args}
        CONTIGS=()
        if [ ${#FAI[@]} -gt 0 ]; then CONTIGS=(--contigs "${FAI[0]}"); fi
        W=$(mktemp -d); cd "$W"
        svtk standardize $ARGS "${CONTIGS[@]}" "${VCF[0]}" {outputs.results}/{inputs.prefix}.vcf.gz \
            {inputs.caller} || exit 1
    """,
)

svtk_vcf2bed_cmd = shell.create(
    name="svtk_vcf2bed",
    image=SVTK_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"vcf": File, "tbi": File | None, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        VCF=({inputs.vcf})
        shopt -s nullglob; TBI=({inputs.tbi}); shopt -u nullglob
        ARGS={inputs.args}
        W=$(mktemp -d); cd "$W"
        for f in "${VCF[0]}" "${TBI[@]}"; do ln -s "$f" "$(basename "$f")"; done
        svtk vcf2bed "$(basename "${VCF[0]}")" {outputs.results}/{inputs.prefix}.bed $ARGS || exit 1
    """,
)

env = flyte.TaskEnvironment.from_task("svtk", svtk_standardize_cmd.as_task(), svtk_vcf2bed_cmd.as_task())


async def pick(results: Dir, name: str, tool: str) -> File:
    found = await results.get_file(name)
    if found is None:
        raise FileNotFoundError(f"svtk {tool} wrote no {name}")
    return found


async def svtk_standardize(vcf: File, caller: str, prefix: str, fai: File | None = None, args: str = "") -> File:
    """``vcf`` (from SV caller ``caller``) in svtk's standard form: ``<prefix>.vcf.gz``."""
    if vcf.name == f"{prefix}.vcf.gz":
        raise ValueError("input and output names are the same; choose another prefix")
    results = await svtk_standardize_cmd(vcf=vcf, fai=fai, caller=caller, prefix=prefix, args=args)
    return await pick(results, f"{prefix}.vcf.gz", "standardize")


async def svtk_vcf2bed(vcf: File, prefix: str, tbi: File | None = None, args: str = "") -> File:
    """``vcf``'s SVs as BED: ``<prefix>.bed``."""
    results = await svtk_vcf2bed_cmd(vcf=vcf, tbi=tbi, prefix=prefix, args=args)
    return await pick(results, f"{prefix}.bed", "vcf2bed")
