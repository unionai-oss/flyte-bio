"""graphtyper — population-scale genotyping with GraphTyper.

Exposes:

- :func:`graphtyper_genotype`, which genotypes small variants in one or more
  BAMs (each with its index) against a reference FASTA (with its ``.fai``),
  in the regions of ``region_file`` or of ``--region`` in ``args`` (one is
  required). As upstream, the BAMs are listed in name order in
  ``bam_list.txt``, and GraphTyper writes one bgzipped VCF per region chunk
  with its ``.tbi``: ``results/<chromosome>/<start>-<end>.vcf.gz``. It uses
  the task's CPUs.
- :func:`graphtyper_vcfconcatenate`, which concatenates VCFs of consecutive
  regions (in the given order) into ``<prefix>.vcf.gz`` with its ``.tbi``. As
  upstream, the VCFs are concatenated in batches of 1000, and the batches
  then into the output.

Upstream pins GraphTyper 2.7.2 for genotyping and 2.7.7 for concatenation, so
each has its own environment and ``env`` aggregates both.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

GENOTYPE_IMAGE = "quay.io/biocontainers/graphtyper:2.7.2--h7d7f7ad_0"
VCFCONCATENATE_IMAGE = "quay.io/biocontainers/graphtyper:2.7.7--h7594796_1"
# Upstream labels genotyping process_medium and concatenation process_low.
GENOTYPE_RESOURCES = flyte.Resources(cpu=4, memory="8Gi")
VCFCONCATENATE_RESOURCES = flyte.Resources(cpu=2, memory="4Gi")
BATCH_SIZE = 1000

# The BAMs and their indexes are linked side by side, under their own names,
# and listed in name order (C collation, as Groovy sorts paths); the reference
# sits next to its .fai.
graphtyper_genotype_cmd = shell.create(
    name="graphtyper_genotype",
    image=GENOTYPE_IMAGE,
    resources=GENOTYPE_RESOURCES,
    inputs={
        "bam": list[File],
        "bai": list[File],
        "ref": File,
        "ref_fai": File,
        "region_file": File | None,
        "args": str,
    },
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        BAM=({{inputs.bam}})
        BAI=({{inputs.bai}})
        REF=({{inputs.ref}})
        FAI=({{inputs.ref_fai}})
        shopt -s nullglob; REGIONS=({{inputs.region_file}}); shopt -u nullglob
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        mkdir inputs
        for F in "${{BAM[@]}}" "${{BAI[@]}}"; do ln -s "$F" inputs/; done
        REF_NAME=$(basename "${{REF[0]}}")
        ln -s "${{REF[0]}}" "$REF_NAME"
        ln -s "${{FAI[0]}}" "$REF_NAME.fai"
        for F in "${{BAM[@]}}"; do echo "inputs/$(basename "$F")"; done | LC_ALL=C sort > bam_list.txt
        REGION=()
        if [ ${{#REGIONS[@]}} -gt 0 ]; then REGION=(--region_file "${{REGIONS[0]}}"); fi
        graphtyper \
            genotype \
            "$REF_NAME" \
            $ARGS \
            --sams bam_list.txt \
            --threads $TASK_CPUS \
            "${{REGION[@]}}" || exit 1
        # As upstream, only results/*/*.vcf.gz (not results/input_sites/...) and their indexes.
        shopt -s nullglob
        for V in results/*/*.vcf.gz; do
            mkdir -p "$OUT/$(dirname "$V")"
            mv "$V" "$V.tbi" "$OUT/$(dirname "$V")"/
        done
    """,
)

# ``names`` carries the caller's order (one name per line): a list[File] is
# staged under the files' names and globbed alphabetically.
graphtyper_vcfconcatenate_cmd = shell.create(
    name="graphtyper_vcfconcatenate",
    image=VCFCONCATENATE_IMAGE,
    resources=VCFCONCATENATE_RESOURCES,
    inputs={"vcf": list[File], "names": str, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        VCF=({{inputs.vcf}})
        NAMES={{inputs.names}}
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        D=$(dirname "${{VCF[0]}}")
        W=$(mktemp -d); cd "$W"
        FILES=()
        while IFS= read -r NAME; do
            ln -s "$D/$NAME" "$NAME"
            FILES+=("$NAME")
        done <<< "$NAMES"
        # Batching, as upstream: too many VCFs overflow the command line.
        INDEX=0
        for (( START=0; START < ${{#FILES[@]}}; START+={BATCH_SIZE} )); do
            graphtyper vcf_concatenate "${{FILES[@]:START:{BATCH_SIZE}}}" $ARGS \
                --output="${{P}}_subset_$INDEX.vcf.gz" || exit 1
            INDEX=$(( INDEX + 1 ))
        done
        graphtyper vcf_concatenate \
            "${{P}}"_subset_*.vcf.gz \
            $ARGS \
            --write_tbi \
            --output="$P.vcf.gz" || exit 1
        rm "${{P}}"_subset_*.vcf.gz
        mv "$P.vcf.gz" "$P.vcf.gz.tbi" "$OUT"/
    """,
)

# One image per environment, so `env` aggregates two.
graphtyper_genotype_env = flyte.TaskEnvironment.from_task("graphtyper_genotype", graphtyper_genotype_cmd.as_task())
graphtyper_vcfconcatenate_env = flyte.TaskEnvironment.from_task(
    "graphtyper_vcfconcatenate", graphtyper_vcfconcatenate_cmd.as_task()
)
env = flyte.TaskEnvironment(name="graphtyper", depends_on=[graphtyper_genotype_env, graphtyper_vcfconcatenate_env])


@dataclass
class GraphtyperGenotypeResult:
    vcf: list[File]  # results/<chromosome>/<start>-<end>.vcf.gz, sorted by path
    tbi: list[File]  # their .tbi indexes, in the same order
    results: Dir  # results/


async def graphtyper_genotype(
    bam: File | list[File],
    bai: File | list[File],
    ref: File,
    ref_fai: File,
    region_file: File | None = None,
    args: str = "",
) -> GraphtyperGenotypeResult:
    """Genotype ``bam`` (one or more, with their ``bai``) in ``region_file``'s regions or ``--region`` in ``args``."""
    bams = [bam] if isinstance(bam, File) else list(bam)
    bais = [bai] if isinstance(bai, File) else list(bai)
    if region_file is None and "region" not in args:
        raise ValueError(
            "graphtyper_genotype requires either a region file or a region specified using '--region' in args"
        )
    names = [f.name for f in bams + bais]
    if len(set(names)) != len(names):
        raise ValueError(f"graphtyper_genotype needs distinct BAM and index names, got {names}")
    out = await graphtyper_genotype_cmd(
        bam=bams, bai=bais, ref=ref, ref_fai=ref_fai, region_file=region_file, args=args
    )
    base = out.path.rstrip("/")
    paths = sorted([f.path[len(base) :].lstrip("/") async for f in out.walk()])
    return GraphtyperGenotypeResult(
        vcf=[File.from_existing_remote(f"{base}/{p}") for p in paths if p.endswith(".vcf.gz")],
        tbi=[File.from_existing_remote(f"{base}/{p}") for p in paths if p.endswith(".vcf.gz.tbi")],
        results=Dir.from_existing_remote(f"{base}/results"),
    )


@dataclass
class GraphtyperVcfConcatenateResult:
    vcf: File  # <prefix>.vcf.gz
    tbi: File  # <prefix>.vcf.gz.tbi


async def graphtyper_vcfconcatenate(vcf: list[File], prefix: str, args: str = "") -> GraphtyperVcfConcatenateResult:
    """Concatenate ``vcf`` (consecutive regions, in order) into ``<prefix>.vcf.gz`` with its ``.tbi``."""
    names = [f.name for f in vcf]
    if not names:
        raise ValueError("graphtyper_vcfconcatenate needs at least one VCF")
    if len(set(names)) != len(names):
        raise ValueError(f"graphtyper_vcfconcatenate needs distinct VCF names, got {names}")
    if f"{prefix}.vcf.gz" in names:
        raise ValueError("Input and output names are the same, set prefix to disambiguate!")
    results = await graphtyper_vcfconcatenate_cmd(vcf=vcf, names="\n".join(names), prefix=prefix, args=args)
    output = await results.get_file(f"{prefix}.vcf.gz")
    index = await results.get_file(f"{prefix}.vcf.gz.tbi")
    if output is None or index is None:
        raise FileNotFoundError(f"graphtyper vcf_concatenate wrote no {prefix}.vcf.gz or its .tbi")
    return GraphtyperVcfConcatenateResult(vcf=output, tbi=index)
