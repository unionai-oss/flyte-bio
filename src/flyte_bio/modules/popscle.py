"""popscle — genotype-based demultiplexing of pooled single-cell data.

Exposes:

- :func:`popscle_dscpileup` — pileup of a droplet-barcoded BAM at a VCF's
  variants (``popscle dsc-pileup``): a Dir of ``<prefix>.{cel,plp,var,umi}.gz``.
- :func:`popscle_demuxlet` — assign droplets to the donors of a genotyped VCF
  (``popscle demuxlet``), from a BAM or from a dsc-pileup: ``<prefix>.best``.
- :func:`popscle_freemuxlet` — cluster droplets into ``n_sample`` donors
  without genotypes (``popscle freemuxlet``) from a dsc-pileup:
  ``<prefix>.clust1.samples.gz``, ``<prefix>.clust1.vcf.gz``,
  ``<prefix>.lmix`` and, with ``--aux-files`` in ``args``,
  ``<prefix>.clust0.samples.gz`` and ``<prefix>.clust0.vcf.gz``.

A pileup is passed between tasks as the Dir dsc-pileup writes; demuxlet and
freemuxlet find its prefix from its ``.plp.gz`` file. ``args`` go to popscle
(e.g. ``--field GT`` for demuxlet).
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

POPSCLE_IMAGE = "quay.io/biocontainers/popscle:0.1beta--h2c78cec_0"
# Upstream: process_medium (dsc-pileup, demuxlet), process_high (freemuxlet); popscle is single-threaded.
DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory="8Gi")
HIGH_RESOURCES = flyte.Resources(cpu=1, memory="16Gi")

# Sets PLP to the pileup Dir's prefix path, found from its .plp.gz file.
FIND_PILEUP = r"""
        PLP_FILES=("$PILEUP"/*.plp.gz)
        if [ ! -e "${PLP_FILES[0]}" ]; then echo "no .plp.gz file in the pileup" >&2; exit 1; fi
        PLP="${PLP_FILES[0]%.plp.gz}"
"""

popscle_dscpileup_cmd = shell.create(
    name="popscle_dscpileup",
    image=POPSCLE_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"bam": File, "vcf": File, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"pileup": Dir},
    script=r"""
        BAM=({inputs.bam}); VCF=({inputs.vcf})
        ARGS={inputs.args}
        popscle dsc-pileup --sam "${BAM[0]}" --vcf "${VCF[0]}" --out {outputs.pileup}/{inputs.prefix} $ARGS || exit 1
    """,
)

popscle_demuxlet_cmd = shell.create(
    name="popscle_demuxlet",
    image=POPSCLE_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"bam": File | None, "pileup": Dir | None, "donor_genotype": File, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        shopt -s nullglob; BAM=({inputs.bam}); shopt -u nullglob
        VCF=({inputs.donor_genotype})
        PILEUP={inputs.pileup}
        ARGS={inputs.args}
        if [ -d "$PILEUP" ]; then
FIND_PILEUP
            INPUT=(--plp "$PLP")
        else
            INPUT=(--sam "${BAM[0]}")
        fi
        popscle demuxlet "${INPUT[@]}" --vcf "${VCF[0]}" --out {outputs.results}/{inputs.prefix} $ARGS || exit 1
    """.replace("FIND_PILEUP", FIND_PILEUP.strip("\n")),
)

popscle_freemuxlet_cmd = shell.create(
    name="popscle_freemuxlet",
    image=POPSCLE_IMAGE,
    resources=HIGH_RESOURCES,
    inputs={"pileup": Dir, "n_sample": int, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        PILEUP={inputs.pileup}
        ARGS={inputs.args}
FIND_PILEUP
        popscle freemuxlet --plp "$PLP" --out {outputs.results}/{inputs.prefix} --nsample {inputs.n_sample} $ARGS \
            || exit 1
    """.replace("FIND_PILEUP", FIND_PILEUP.strip("\n")),
)

env = flyte.TaskEnvironment.from_task(
    "popscle",
    popscle_dscpileup_cmd.as_task(),
    popscle_demuxlet_cmd.as_task(),
    popscle_freemuxlet_cmd.as_task(),
)


async def pick(results: Dir, name: str, tool: str) -> File:
    found = await results.get_file(name)
    if found is None:
        raise FileNotFoundError(f"popscle {tool} wrote no {name}")
    return found


async def popscle_dscpileup(bam: File, vcf: File, prefix: str, args: str = "") -> Dir:
    """``popscle dsc-pileup`` of ``bam`` at ``vcf``'s variants: a Dir of ``<prefix>.{cel,plp,var,umi}.gz``."""
    return await popscle_dscpileup_cmd(bam=bam, vcf=vcf, prefix=prefix, args=args)


async def popscle_demuxlet(
    donor_genotype: File,
    prefix: str,
    bam: File | None = None,
    pileup: Dir | None = None,
    args: str = "",
) -> File:
    """``popscle demuxlet`` against ``donor_genotype``: ``<prefix>.best``.

    Reads come from a dsc-pileup Dir (``pileup``) or from ``bam``; pass one.
    """
    if (bam is None) == (pileup is None):
        raise ValueError("popscle_demuxlet needs exactly one of bam and pileup")
    results = await popscle_demuxlet_cmd(
        bam=bam, pileup=pileup, donor_genotype=donor_genotype, prefix=prefix, args=args
    )
    return await pick(results, f"{prefix}.best", "demuxlet")


@dataclass
class FreemuxletResult:
    samples: File  # <prefix>.clust1.samples.gz: each droplet's cluster
    vcf: File  # <prefix>.clust1.vcf.gz: the clusters' genotypes
    lmix: File  # <prefix>.lmix
    singlet_samples: File | None  # <prefix>.clust0.samples.gz, with --aux-files
    singlet_vcf: File | None  # <prefix>.clust0.vcf.gz, with --aux-files
    results: Dir


async def popscle_freemuxlet(pileup: Dir, n_sample: int, prefix: str, args: str = "") -> FreemuxletResult:
    """``popscle freemuxlet``: cluster a dsc-pileup's droplets into ``n_sample`` donors."""
    results = await popscle_freemuxlet_cmd(pileup=pileup, n_sample=n_sample, prefix=prefix, args=args)
    return FreemuxletResult(
        samples=await pick(results, f"{prefix}.clust1.samples.gz", "freemuxlet"),
        vcf=await pick(results, f"{prefix}.clust1.vcf.gz", "freemuxlet"),
        lmix=await pick(results, f"{prefix}.lmix", "freemuxlet"),
        singlet_samples=await results.get_file(f"{prefix}.clust0.samples.gz"),
        singlet_vcf=await results.get_file(f"{prefix}.clust0.vcf.gz"),
        results=results,
    )
