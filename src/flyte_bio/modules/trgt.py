"""trgt — tandem-repeat genotyping and visualisation from PacBio HiFi reads with TRGT.

Exposes:

- :func:`trgt_genotype` — genotype the repeats in a BED of repeat
  definitions from one sample's aligned BAM (with its index), against the
  reference and its ``.fai``: ``<prefix>.vcf.gz`` and the reads spanning
  each repeat, ``<prefix>.spanning.bam`` (unsorted). ``karyotype`` (``XX``,
  ``XY``, or a file of per-contig ploidies) sets the sex-chromosome
  ploidy; empty leaves TRGT's default.
- :func:`trgt_plot` — plot one repeat (``repeat_id``) from TRGT's
  genotypes (a sorted, indexed VCF) and its sorted, indexed spanning reads:
  ``<prefix>.png``, unless ``args`` sets the image itself (``--image`` or
  ``-o``), in which case the plot is whatever ``.png``, ``.pdf`` or ``.svg``
  it names.

``--threads`` for genotyping follows the task's CPUs. ``args`` go to the
tool and aren't shell-quoted. The inputs are linked into the working
directory under their own names, indexes beside their files.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

TRGT_IMAGE = "community.wave.seqera.io/library/trgt:5.1.0--66e3b26326e566e7"
# Upstream labels genotype process_medium and plot process_low; its tests run with 4 CPUs.
GENOTYPE_RESOURCES = flyte.Resources(cpu=4, memory="8Gi")
PLOT_RESOURCES = flyte.Resources(cpu=1, memory="4Gi")

# Links FA as $FA_NAME with FAI beside it.
LINK_REFERENCE = r"""
        FA_NAME=$(basename "${FA[0]}"); ln -s "${FA[0]}" "$FA_NAME"
        ln -s "${FAI[0]}" "$FA_NAME.fai"
"""

trgt_genotype_cmd = shell.create(
    name="trgt_genotype",
    image=TRGT_IMAGE,
    resources=GENOTYPE_RESOURCES,
    inputs={
        "bam": File,
        "bai": File,
        "karyotype": str,
        "fasta": File,
        "fai": File,
        "repeats": File,
        "prefix": str,
        "args": str,
    },
    defaults={"karyotype": "", "args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        BAM=({{inputs.bam}}); BAI=({{inputs.bai}}); FA=({{inputs.fasta}}); FAI=({{inputs.fai}})
        REPEATS=({{inputs.repeats}})
        KARYOTYPE={{inputs.karyotype}}
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        W=$(mktemp -d); cd "$W"
        BAM_NAME=$(basename "${{BAM[0]}}"); ln -s "${{BAM[0]}}" "$BAM_NAME"
        ln -s "${{BAI[0]}}" "$BAM_NAME.${{BAI[0]##*.}}"
        {LINK_REFERENCE}
        REPEATS_NAME=$(basename "${{REPEATS[0]}}"); ln -s "${{REPEATS[0]}}" "$REPEATS_NAME"
        EXTRA=()
        if [ -n "$KARYOTYPE" ]; then EXTRA+=(--karyotype "$KARYOTYPE"); fi

        trgt genotype \
            $ARGS \
            --genome "$FA_NAME" \
            --reads "$BAM_NAME" \
            --repeats "$REPEATS_NAME" \
            "${{EXTRA[@]}}" \
            --threads $TASK_CPUS \
            --output-prefix "$P" || exit 1

        shopt -s nullglob
        mv "$P".vcf.gz "$P".spanning.ba[m] {{outputs.results}}/
    """,
)

trgt_plot_cmd = shell.create(
    name="trgt_plot",
    image=TRGT_IMAGE,
    resources=PLOT_RESOURCES,
    inputs={
        "bam": File,
        "bai": File,
        "vcf": File,
        "tbi": File,
        "repeat_id": str,
        "fasta": File,
        "fai": File,
        "repeats": File,
        "prefix": str,
        "args": str,
    },
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        BAM=({{inputs.bam}}); BAI=({{inputs.bai}}); VCF=({{inputs.vcf}}); TBI=({{inputs.tbi}})
        FA=({{inputs.fasta}}); FAI=({{inputs.fai}}); REPEATS=({{inputs.repeats}})
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        W=$(mktemp -d); mkdir "$W/out"; cd "$W"
        ln -s "${{BAM[0]}}" spanning.bam; ln -s "${{BAI[0]}}" spanning.bam.bai
        ln -s "${{VCF[0]}}" genotypes.vcf.gz; ln -s "${{TBI[0]}}" genotypes.vcf.gz.${{TBI[0]##*.}}
        {LINK_REFERENCE}
        REPEATS_NAME=$(basename "${{REPEATS[0]}}"); ln -s "${{REPEATS[0]}}" "$REPEATS_NAME"
        # As upstream, plot to <prefix>.png unless args name the image.
        IMAGE=(--image "out/$P.png")
        if [[ " $ARGS " == *--image* || " $ARGS " == *-o* ]]; then IMAGE=(); fi

        trgt plot \
            $ARGS \
            --genome "$FA_NAME" \
            --repeats "$REPEATS_NAME" \
            --spanning-reads spanning.bam \
            --vcf genotypes.vcf.gz \
            --repeat-id {{inputs.repeat_id}} \
            "${{IMAGE[@]}}" || exit 1

        shopt -s nullglob
        mv ./*.png ./*.pdf ./*.svg out/*.png {{outputs.results}}/
    """,
)

env = flyte.TaskEnvironment.from_task("trgt", trgt_genotype_cmd.as_task(), trgt_plot_cmd.as_task())


@dataclass
class TrgtGenotypeResult:
    vcf: File  # <prefix>.vcf.gz
    spanning_bam: File | None  # <prefix>.spanning.bam
    results: Dir  # everything above


async def trgt_genotype(
    bam: File,
    bai: File,
    fasta: File,
    fai: File,
    repeats: File,
    karyotype: str = "",
    prefix: str = "sample",
    args: str = "",
) -> TrgtGenotypeResult:
    """Genotype the tandem repeats in ``repeats`` (a TRGT repeat BED) from ``bam`` against ``fasta``."""
    results = await trgt_genotype_cmd(
        bam=bam, bai=bai, karyotype=karyotype, fasta=fasta, fai=fai, repeats=repeats, prefix=prefix, args=args
    )
    vcf = await results.get_file(f"{prefix}.vcf.gz")
    if vcf is None:
        raise FileNotFoundError(f"TRGT wrote no {prefix}.vcf.gz")
    return TrgtGenotypeResult(vcf=vcf, spanning_bam=await results.get_file(f"{prefix}.spanning.bam"), results=results)


async def trgt_plot(
    bam: File,
    bai: File,
    vcf: File,
    tbi: File,
    repeat_id: str,
    fasta: File,
    fai: File,
    repeats: File,
    prefix: str = "sample",
    args: str = "",
) -> File:
    """Plot repeat ``repeat_id`` from sorted, indexed genotypes and spanning reads: ``<prefix>.png`` by default."""
    results = await trgt_plot_cmd(
        bam=bam,
        bai=bai,
        vcf=vcf,
        tbi=tbi,
        repeat_id=repeat_id,
        fasta=fasta,
        fai=fai,
        repeats=repeats,
        prefix=prefix,
        args=args,
    )
    plots = [f async for f in results.walk()]
    if len(plots) != 1:
        raise FileNotFoundError(f"TRGT plot wrote {len(plots)} images, expected one")
    return plots[0]
