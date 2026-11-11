"""stranger — annotate repeat-expansion calls with their pathogenicity.

Exposes :func:`stranger`, which annotates a VCF of short tandem repeat
genotypes (from ExpansionHunter or TRGT) with each locus's normal and
pathologic repeat-count thresholds and a status (``normal``,
``pre_mutation``, ``full_mutation``), read from a repeat catalogue
(``--repeats-file``; stranger's bundled catalogue without one).

As upstream, the annotated VCF is bgzipped and tabix-indexed:
``<prefix>.vcf.gz`` and ``<prefix>.vcf.gz.tbi``, with ``prefix`` defaulting
to ``<sample>_stranger``. ``args`` go to stranger, ``args2`` to ``bgzip`` and
``args3`` to ``tabix``; both get the task's CPUs as threads.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

STRANGER_IMAGE = "community.wave.seqera.io/library/htslib_pip_python_stranger:beca035d9e9d6787"
DEFAULT_RESOURCES = flyte.Resources(cpu=2, memory="4Gi")

stranger_cmd = shell.create(
    name="stranger",
    image=STRANGER_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"vcf": File, "variant_catalog": File | None, "prefix": str, "args": str, "args2": str, "args3": str},
    defaults={"args": "", "args2": "", "args3": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        set -o pipefail
        VCF=({{inputs.vcf}})
        shopt -s nullglob; CATALOG=({{inputs.variant_catalog}}); shopt -u nullglob
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        ARGS2={{inputs.args2}}
        ARGS3={{inputs.args3}}
        W=$(mktemp -d); cd "$W"
        for F in "${{VCF[0]}}" "${{CATALOG[@]}}"; do ln -s "$F" "$(basename "$F")"; done
        REPEATS=()
        if [ ${{#CATALOG[@]}} -gt 0 ]; then REPEATS=(--repeats-file "$(basename "${{CATALOG[0]}}")"); fi
        stranger $ARGS "$(basename "${{VCF[0]}}")" "${{REPEATS[@]}}" \
            | bgzip $ARGS2 -c --threads $TASK_CPUS > "$P.vcf.gz" || exit 1
        tabix $ARGS3 --threads $TASK_CPUS "$P.vcf.gz" || exit 1
        mv "$P.vcf.gz" "$P.vcf.gz.tbi" {{outputs.results}}/
    """,
)

env = flyte.TaskEnvironment.from_task("stranger", stranger_cmd.as_task())


@dataclass
class StrangerResult:
    vcf: File  # <prefix>.vcf.gz
    tbi: File  # <prefix>.vcf.gz.tbi
    results: Dir


async def stranger(
    vcf: File,
    variant_catalog: File | None = None,
    prefix: str = "sample_stranger",
    args: str = "",
    args2: str = "",
    args3: str = "",
) -> StrangerResult:
    """Annotate the repeat-expansion calls in ``vcf`` against ``variant_catalog``."""
    name = vcf.name or vcf.path.rsplit("/", 1)[-1]
    if name == f"{prefix}.vcf.gz":
        raise ValueError("stranger: input and output names are the same; pass another prefix")
    results = await stranger_cmd(
        vcf=vcf, variant_catalog=variant_catalog, prefix=prefix, args=args, args2=args2, args3=args3
    )
    out = {}
    for key, file_name in (("vcf", f"{prefix}.vcf.gz"), ("tbi", f"{prefix}.vcf.gz.tbi")):
        found = await results.get_file(file_name)
        if found is None:
            raise FileNotFoundError(f"stranger wrote no {file_name}")
        out[key] = found
    return StrangerResult(vcf=out["vcf"], tbi=out["tbi"], results=results)
