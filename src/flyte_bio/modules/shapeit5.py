"""shapeit5 — SHAPEIT5 haplotype phasing.

Exposes two tools of SHAPEIT 5.1.1:

- :func:`shapeit5_phasecommon` — phase the common variants of a region of a
  VCF/BCF, optionally with a reference panel (``--reference``), a scaffold of
  already-phased haplotypes (``--scaffold``), a pedigree (``--pedigree``)
  and a genetic map (``--map``), into ``<prefix>.<extension>``. As upstream,
  the extension follows ``--output-format`` in ``args`` (``bcf``, ``graph``
  or ``bh``; default ``bcf``).
- :func:`shapeit5_ligate` — ligate phased chunks (each with its index) into
  ``<prefix>.<suffix>`` (``vcf``, ``vcf.gz``, ``bcf`` or ``bcf.gz``). As
  upstream, the chunks are listed in name order.

Inputs are linked into a working directory under their own names, indexes
beside their files. ``--thread`` follows the task's CPUs.
"""

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

SHAPEIT5_IMAGE = "community.wave.seqera.io/library/shapeit5:5.1.1--09a6cb254ece8f6e"
# Upstream's process_medium (phase_common) and process_low (ligate).
MEDIUM_RESOURCES = flyte.Resources(cpu=4, memory="8Gi")
LOW_RESOURCES = flyte.Resources(cpu=2, memory="4Gi")

LINK_HELPERS = r"""
        link() { ln -s "$1" "$(basename "$1")"; basename "$1"; }
        link_index() { ln -s "$2" "$1.${2##*.}"; }
        W=$(mktemp -d); cd "$W"
"""

shapeit5_phasecommon_cmd = shell.create(
    name="shapeit5_phasecommon",
    image=SHAPEIT5_IMAGE,
    resources=MEDIUM_RESOURCES,
    inputs={
        "vcf": File,
        "index": File,
        "pedigree": File | None,
        "region": str,
        "reference": File | None,
        "reference_index": File | None,
        "scaffold": File | None,
        "scaffold_index": File | None,
        "genetic_map": File | None,
        "prefix": str,
        "extension": str,
        "args": str,
    },
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        V=({{inputs.vcf}}); I=({{inputs.index}})
        shopt -s nullglob
        PED=({{inputs.pedigree}}); REF=({{inputs.reference}}); REF_I=({{inputs.reference_index}})
        SC=({{inputs.scaffold}}); SC_I=({{inputs.scaffold_index}}); M=({{inputs.genetic_map}})
        shopt -u nullglob
        REGION={{inputs.region}}; P={{inputs.prefix}}; EXT={{inputs.extension}}
        ARGS={{inputs.args}}
        {LINK_HELPERS}
        NAME=$(link "${{V[0]}}"); link_index "$NAME" "${{I[0]}}"
        EXTRA=()
        if [ ${{#M[@]}} -gt 0 ]; then EXTRA+=(--map "$(link "${{M[0]}}")"); fi
        if [ ${{#REF[@]}} -gt 0 ]; then
            REF_NAME=$(link "${{REF[0]}}"); EXTRA+=(--reference "$REF_NAME")
            if [ ${{#REF_I[@]}} -gt 0 ]; then link_index "$REF_NAME" "${{REF_I[0]}}"; fi
        fi
        if [ ${{#SC[@]}} -gt 0 ]; then
            SC_NAME=$(link "${{SC[0]}}"); EXTRA+=(--scaffold "$SC_NAME")
            if [ ${{#SC_I[@]}} -gt 0 ]; then link_index "$SC_NAME" "${{SC_I[0]}}"; fi
        fi
        if [ ${{#PED[@]}} -gt 0 ]; then EXTRA+=(--pedigree "$(link "${{PED[0]}}")"); fi
        SHAPEIT5_phase_common $ARGS --input "$NAME" "${{EXTRA[@]}}" --region "$REGION" \
            --thread $TASK_CPUS --output "$P.$EXT" || exit 1
        mv "$P.$EXT" {{outputs.results}}/
    """,
)

# As upstream, the chunks are listed (all_files.txt) in name order.
shapeit5_ligate_cmd = shell.create(
    name="shapeit5_ligate",
    image=SHAPEIT5_IMAGE,
    resources=LOW_RESOURCES,
    inputs={"vcfs": list[File], "indexes": list[File], "prefix": str, "suffix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        shopt -s nullglob; VCFS=({{inputs.vcfs}}); IDX=({{inputs.indexes}}); shopt -u nullglob
        P={{inputs.prefix}}; SUFFIX={{inputs.suffix}}
        ARGS={{inputs.args}}
        {LINK_HELPERS}
        for F in "${{VCFS[@]}}" "${{IDX[@]}}"; do link "$F" > /dev/null; done
        for F in "${{VCFS[@]}}"; do basename "$F"; done | LC_ALL=C sort > all_files.txt
        SHAPEIT5_ligate $ARGS --input all_files.txt --thread $TASK_CPUS --output "$P.$SUFFIX" || exit 1
        mv "$P.$SUFFIX" {{outputs.results}}/
    """,
)

env = flyte.TaskEnvironment.from_task("shapeit5", shapeit5_phasecommon_cmd.as_task(), shapeit5_ligate_cmd.as_task())


def phasecommon_extension(args: str) -> str:
    """The output's extension, read from ``--output-format`` in ``args`` as upstream does."""
    for extension in ("bcf", "graph", "bh"):
        if f"--output-format {extension}" in args:
            return extension
    return "bcf"


async def shapeit5_phasecommon(
    vcf: File,
    index: File,
    region: str,
    prefix: str,
    pedigree: File | None = None,
    reference: File | None = None,
    reference_index: File | None = None,
    scaffold: File | None = None,
    scaffold_index: File | None = None,
    genetic_map: File | None = None,
    args: str = "",
) -> File:
    """Phase the common variants of ``vcf`` in ``region``: ``<prefix>.<extension>`` (BCF by default).

    ``reference``, ``scaffold`` (each with its index), ``pedigree`` and
    ``genetic_map`` are passed as ``--reference``, ``--scaffold``,
    ``--pedigree`` and ``--map``.
    """
    extension = phasecommon_extension(args)
    if vcf.name == f"{prefix}.{extension}":
        raise ValueError("Input and output names are the same, set prefix to disambiguate!")
    results = await shapeit5_phasecommon_cmd(
        vcf=vcf,
        index=index,
        pedigree=pedigree,
        region=region,
        reference=reference,
        reference_index=reference_index,
        scaffold=scaffold,
        scaffold_index=scaffold_index,
        genetic_map=genetic_map,
        prefix=prefix,
        extension=extension,
        args=args,
    )
    found = await results.get_file(f"{prefix}.{extension}")
    if found is None:
        raise FileNotFoundError(f"SHAPEIT5_phase_common wrote no {prefix}.{extension}")
    return found


async def shapeit5_ligate(
    vcfs: list[File], indexes: list[File], prefix: str, suffix: str = "vcf.gz", args: str = ""
) -> File:
    """Ligate phased chunks ``vcfs`` (with their ``indexes``, named after them) into ``<prefix>.<suffix>``."""
    if suffix not in ("vcf", "vcf.gz", "bcf", "bcf.gz"):
        raise ValueError(f"suffix must be one of vcf, vcf.gz, bcf, bcf.gz, got {suffix}")
    names = [f.name for f in vcfs]
    if len(set(names)) != len(names):
        raise ValueError(f"SHAPEIT5_ligate inputs need distinct names, got {names}")
    results = await shapeit5_ligate_cmd(vcfs=vcfs, indexes=indexes, prefix=prefix, suffix=suffix, args=args)
    found = await results.get_file(f"{prefix}.{suffix}")
    if found is None:
        raise FileNotFoundError(f"SHAPEIT5_ligate wrote no {prefix}.{suffix}")
    return found
