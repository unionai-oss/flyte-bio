"""beagle5 — Beagle 5 genotype phasing and imputation.

Exposes :func:`beagle5_beagle`, which phases a VCF's genotypes (``gt=``) and,
given a reference panel (``ref=``), imputes the panel's ungenotyped markers.
Optional inputs are a PLINK genetic map (``map=``), a chromosome or region
(``chrom=``) and files of samples and markers to exclude
(``excludesamples=``, ``excludemarkers=``). Other Beagle parameters go in
``args`` as ``key=value`` pairs.

Outputs keep upstream's names: ``<prefix>.vcf.gz`` and Beagle's log,
``<prefix>.log`` (upstream's default prefix is ``<id>.bglout``).

As upstream, Java's heap (``-Xmx``) is 80% of the task's memory. Beagle's
``nthreads`` follows the task's CPUs (upstream leaves it to Beagle, which
uses every core the machine has).
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import JAVA_HEAP_MB, TASK_LIMITS

BEAGLE_IMAGE = "community.wave.seqera.io/library/beagle:5.5_27Feb25.75f--30bdd31e53923776"
# Upstream labels Beagle process_high; its tests ran capped at 4 CPUs.
DEFAULT_RESOURCES = flyte.Resources(cpu=4, memory="8Gi")

beagle5_beagle_cmd = shell.create(
    name="beagle5_beagle",
    image=BEAGLE_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "vcf": File,
        "vcf_index": File | None,
        "reference": File | None,
        "reference_index": File | None,
        "genetic_map": File | None,
        "exclude_samples": File | None,
        "exclude_markers": File | None,
        "region": str,
        "prefix": str,
        "args": str,
    },
    defaults={"region": "", "args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        V=({{inputs.vcf}})
        shopt -s nullglob
        VI=({{inputs.vcf_index}}); REF=({{inputs.reference}}); REF_I=({{inputs.reference_index}})
        M=({{inputs.genetic_map}}); XS=({{inputs.exclude_samples}}); XM=({{inputs.exclude_markers}})
        shopt -u nullglob
        REGION={{inputs.region}}; P={{inputs.prefix}}
        ARGS={{inputs.args}}
        W=$(mktemp -d); cd "$W"
        link() {{ ln -s "$1" "$(basename "$1")"; basename "$1"; }}
        NAME=$(link "${{V[0]}}")
        for F in "${{VI[@]}}" "${{REF_I[@]}}"; do link "$F" > /dev/null; done
        EXTRA=()
        if [ ${{#REF[@]}} -gt 0 ]; then EXTRA+=("ref=$(link "${{REF[0]}}")"); fi
        if [ ${{#M[@]}} -gt 0 ]; then EXTRA+=("map=$(link "${{M[0]}}")"); fi
        if [ -n "$REGION" ]; then EXTRA+=("chrom=$REGION"); fi
        if [ ${{#XS[@]}} -gt 0 ]; then EXTRA+=("excludesamples=$(link "${{XS[0]}}")"); fi
        if [ ${{#XM[@]}} -gt 0 ]; then EXTRA+=("excludemarkers=$(link "${{XM[0]}}")"); fi
        beagle -Xmx{JAVA_HEAP_MB}M gt="$NAME" out="$P" nthreads=$TASK_CPUS $ARGS "${{EXTRA[@]}}" || exit 1
        mv "$P.vcf.gz" "$P.log" {{outputs.results}}/
    """,
)

env = flyte.TaskEnvironment.from_task("beagle5", beagle5_beagle_cmd.as_task())


@dataclass
class BeagleResult:
    vcf: File  # <prefix>.vcf.gz: phased (and imputed) genotypes
    log: File  # <prefix>.log
    results: Dir


async def beagle5_beagle(
    vcf: File,
    prefix: str,
    vcf_index: File | None = None,
    reference: File | None = None,
    reference_index: File | None = None,
    genetic_map: File | None = None,
    exclude_samples: File | None = None,
    exclude_markers: File | None = None,
    region: str = "",
    args: str = "",
) -> BeagleResult:
    """Phase ``vcf`` (and impute it from the panel ``reference``): ``<prefix>.vcf.gz`` and ``<prefix>.log``.

    ``region`` (``chrom=``, e.g. ``chr22`` or ``chr22:1-1000000``) restricts
    the analysis, e.g. to the chromosome a single-chromosome panel covers.
    """
    results = await beagle5_beagle_cmd(
        vcf=vcf,
        vcf_index=vcf_index,
        reference=reference,
        reference_index=reference_index,
        genetic_map=genetic_map,
        exclude_samples=exclude_samples,
        exclude_markers=exclude_markers,
        region=region,
        prefix=prefix,
        args=args,
    )
    output, log = await results.get_file(f"{prefix}.vcf.gz"), await results.get_file(f"{prefix}.log")
    if output is None or log is None:
        raise FileNotFoundError(f"beagle wrote no {prefix}.vcf.gz and {prefix}.log")
    return BeagleResult(vcf=output, log=log, results=results)
