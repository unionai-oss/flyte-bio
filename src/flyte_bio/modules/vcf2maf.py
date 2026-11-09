"""vcf2maf — conversion of a VCF to a Mutation Annotation Format (MAF) table.

Exposes :func:`vcf2maf`, which runs ``vcf2maf.pl`` on an uncompressed VCF
against ``fasta`` (for the flanking reference context) and writes
``<prefix>.maf``. With ``vep_cache`` (a directory holding
``<species>/<release>_<assembly>/``) vcf2maf first annotates the variants
with Ensembl VEP from the image; without one, pass ``--inhibit-vep`` in
``args`` to convert the VCF as it is. ``args`` takes vcf2maf's options, e.g.
``--ncbi-build GRCh38 --species homo_sapiens``.

As upstream, VEP is found next to ``vep`` on the ``PATH`` when a cache is
given. vcf2maf writes its intermediate ``.vep.vcf`` and the FASTA's ``.fai``
next to its inputs, so both are symlinked into a working directory first.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

VCF2MAF_IMAGE = "community.wave.seqera.io/library/ensembl-vep_vcf2maf:2d40b60b4834af73"
# Upstream's process_low; vcf2maf runs VEP with its default 4 forks.
DEFAULT_RESOURCES = flyte.Resources(cpu=4, memory="8Gi")

vcf2maf_cmd = shell.create(
    name="vcf2maf",
    image=VCF2MAF_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"vcf": File, "fasta": File, "vep_cache": Dir | None, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        VCF=({inputs.vcf}); FA=({inputs.fasta})
        CACHE={inputs.vep_cache}
        P={inputs.prefix}
        ARGS={inputs.args}
        W=$(mktemp -d); cd "$W"
        ln -s "${VCF[0]}" "${FA[0]}" .
        VEP_CMD=()
        VEP_CACHE_CMD=()
        if [ -d "$CACHE" ]; then
            VEP_CMD=(--vep-path "$(dirname "$(type -p vep)")")
            VEP_CACHE_CMD=(--vep-data "$CACHE")
        fi
        vcf2maf.pl \
            $ARGS \
            "${VEP_CMD[@]}" \
            "${VEP_CACHE_CMD[@]}" \
            --ref-fasta "$(basename "${FA[0]}")" \
            --input-vcf "$(basename "${VCF[0]}")" \
            --output-maf "$P.maf" || exit 1
        mv "$P.maf" {outputs.results}/
    """,
)

env = flyte.TaskEnvironment.from_task("vcf2maf", vcf2maf_cmd.as_task())


@dataclass
class Vcf2mafResult:
    maf: File  # <prefix>.maf
    results: Dir


async def vcf2maf(
    vcf: File, fasta: File, vep_cache: Dir | None = None, prefix: str = "sample", args: str = ""
) -> Vcf2mafResult:
    """Convert the uncompressed ``vcf`` to ``<prefix>.maf``, annotating with VEP when ``vep_cache`` is given."""
    results = await vcf2maf_cmd(vcf=vcf, fasta=fasta, vep_cache=vep_cache, prefix=prefix, args=args)
    maf = await results.get_file(f"{prefix}.maf")
    if maf is None:
        raise FileNotFoundError(f"vcf2maf wrote no {prefix}.maf")
    return Vcf2mafResult(maf=maf, results=results)
