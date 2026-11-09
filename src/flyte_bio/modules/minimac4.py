"""minimac4 — genotype imputation with Minimac4.

Exposes the two steps of Minimac4 4.1.6:

- :func:`minimac4_compressref` — compress a phased reference panel (VCF/BCF,
  whose index Minimac4 finds beside it) into ``<prefix>.msav``.
- :func:`minimac4_impute` — impute a target VCF from a compressed panel,
  optionally with a genetic map (``--map``) and within a region
  (``--region``), into ``<prefix>.<extension>``; with ``write_sites`` it also
  writes the imputed sites' summary, ``<prefix>.sites.vcf.gz``. As upstream,
  the extension follows ``--output-format``/``-O`` in ``args`` (``bcf``,
  ``sav``, ``vcf.gz``, ``vcf``, ``ubcf``, ``usav``; default ``vcf.gz``).

Inputs are linked into a working directory under their own names, indexes
beside their files. ``--threads`` follows the task's CPUs.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

MINIMAC4_IMAGE = "quay.io/biocontainers/minimac4:4.1.6--hcb620b3_1"
# Upstream's process_low (compressref) and process_medium (impute), at the
# CPUs its tests run with.
LOW_RESOURCES = flyte.Resources(cpu=2, memory="4Gi")
MEDIUM_RESOURCES = flyte.Resources(cpu=4, memory="8Gi")

LINK_HELPERS = r"""
        link() { ln -s "$1" "$(basename "$1")"; basename "$1"; }
        link_index() { ln -s "$2" "$1.${2##*.}"; }
        W=$(mktemp -d); cd "$W"
"""

minimac4_compressref_cmd = shell.create(
    name="minimac4_compressref",
    image=MINIMAC4_IMAGE,
    resources=LOW_RESOURCES,
    inputs={"reference": File, "reference_index": File | None, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        REF=({{inputs.reference}})
        shopt -s nullglob; REF_I=({{inputs.reference_index}}); shopt -u nullglob
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        {LINK_HELPERS}
        NAME=$(link "${{REF[0]}}")
        if [ ${{#REF_I[@]}} -gt 0 ]; then link_index "$NAME" "${{REF_I[0]}}"; fi
        minimac4 --compress-reference "$NAME" $ARGS --threads $TASK_CPUS -o "$P.msav" || exit 1
        mv "$P.msav" {{outputs.results}}/
    """,
)

minimac4_impute_cmd = shell.create(
    name="minimac4_impute",
    image=MINIMAC4_IMAGE,
    resources=MEDIUM_RESOURCES,
    inputs={
        "target": File,
        "target_index": File | None,
        "reference": File,
        "genetic_map": File | None,
        "region": str,
        "write_sites": bool,
        "prefix": str,
        "extension": str,
        "args": str,
    },
    defaults={"region": "", "write_sites": False, "args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        T=({{inputs.target}}); REF=({{inputs.reference}})
        shopt -s nullglob; TI=({{inputs.target_index}}); M=({{inputs.genetic_map}}); shopt -u nullglob
        REGION={{inputs.region}}; P={{inputs.prefix}}; EXT={{inputs.extension}}
        ARGS={{inputs.args}}
        {LINK_HELPERS}
        NAME=$(link "${{T[0]}}")
        if [ ${{#TI[@]}} -gt 0 ]; then link_index "$NAME" "${{TI[0]}}"; fi
        REF_NAME=$(link "${{REF[0]}}")
        EXTRA=()
        if [ ${{#M[@]}} -gt 0 ]; then EXTRA+=(--map "$(link "${{M[0]}}")"); fi
        if [ -n "$REGION" ]; then EXTRA+=(--region "$REGION"); fi
        if [ {{inputs.write_sites}} = true ]; then EXTRA+=(--sites "$P.sites.vcf.gz"); fi
        minimac4 "$REF_NAME" "$NAME" $ARGS "${{EXTRA[@]}}" --threads $TASK_CPUS -o "$P.$EXT" || exit 1
        mv "$P.$EXT" {{outputs.results}}/
        if [ -f "$P".sites.vcf.gz ]; then mv "$P".sites.vcf.gz {{outputs.results}}/; fi
    """,
)

env = flyte.TaskEnvironment.from_task("minimac4", minimac4_compressref_cmd.as_task(), minimac4_impute_cmd.as_task())


async def minimac4_compressref(
    reference: File, prefix: str, reference_index: File | None = None, args: str = ""
) -> File:
    """``<prefix>.msav``: the phased panel ``reference`` in Minimac4's compressed format."""
    results = await minimac4_compressref_cmd(
        reference=reference, reference_index=reference_index, prefix=prefix, args=args
    )
    found = await results.get_file(f"{prefix}.msav")
    if found is None:
        raise FileNotFoundError(f"minimac4 wrote no {prefix}.msav")
    return found


def impute_extension(args: str) -> str:
    """The imputed output's extension, read from ``--output-format``/``-O`` in ``args`` as upstream does."""
    for extension in ("bcf", "sav", "vcf.gz", "vcf", "ubcf", "usav"):
        if f"--output-format {extension}" in args or f"-O {extension}" in args:
            return extension
    return "vcf.gz"


@dataclass
class Minimac4Result:
    vcf: File  # <prefix>.<extension>: the imputed genotypes
    sites: File | None  # <prefix>.sites.vcf.gz, with write_sites
    results: Dir


async def minimac4_impute(
    target: File,
    reference: File,
    prefix: str,
    target_index: File | None = None,
    genetic_map: File | None = None,
    region: str = "",
    write_sites: bool = False,
    args: str = "",
) -> Minimac4Result:
    """Impute ``target`` from the compressed panel ``reference`` (an ``.msav``): ``<prefix>.<extension>``.

    ``genetic_map`` is passed as ``--map`` and ``region`` as ``--region``;
    ``write_sites`` adds ``<prefix>.sites.vcf.gz`` (``--sites``).
    """
    extension = impute_extension(args)
    results = await minimac4_impute_cmd(
        target=target,
        target_index=target_index,
        reference=reference,
        genetic_map=genetic_map,
        region=region,
        write_sites=write_sites,
        prefix=prefix,
        extension=extension,
        args=args,
    )
    imputed = await results.get_file(f"{prefix}.{extension}")
    if imputed is None:
        raise FileNotFoundError(f"minimac4 wrote no {prefix}.{extension}")
    sites = await results.get_file(f"{prefix}.sites.vcf.gz") if write_sites else None
    if write_sites and sites is None:
        raise FileNotFoundError(f"minimac4 wrote no {prefix}.sites.vcf.gz")
    return Minimac4Result(vcf=imputed, sites=sites, results=results)
