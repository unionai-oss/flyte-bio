"""jasminesv — Jasmine merging of structural-variant calls across samples or callers.

Exposes :func:`jasminesv`, which merges several VCFs of structural variants
(plain or gzipped) into one, ``<prefix>.vcf.gz``. Optionally:

- ``bams`` (with ``bais``), ``fasta`` and ``--run_iris`` in ``args`` refine
  insertion sequences and positions with Iris (``iris_args`` from ``args2``,
  a comma-separated list such as ``threads=1,--pacbio``);
- ``sample_dists`` gives per-sample merging distances;
- ``chr_norm`` with ``--normalize_chrs`` maps chromosome names.

As upstream:

- Gzipped VCFs are decompressed first (``bgzip -d``, with ``args4``), and
  Jasmine reads every VCF in the working directory in name order, as does its
  BAM list.
- ``threads`` follows the task's CPUs; other options (``min_support=2``,
  ``max_dist=…``, ``--output_genotypes``) go in ``args``.
- The merged VCF is bgzipped (with ``args3``). When nothing merges, Jasmine
  writes an empty file, and the output is the first VCF's header instead.
"""

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

JASMINESV_IMAGE = "quay.io/biocontainers/jasminesv:1.1.5--hdfd78af_0"

# Upstream labels it process_single.
DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory="6Gi")

jasminesv_cmd = shell.create(
    name="jasminesv",
    image=JASMINESV_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "vcfs": list[File],
        "bams": list[File],
        "bais": list[File],
        "sample_dists": File | None,
        "fasta": File | None,
        "fai": File | None,
        "chr_norm": File | None,
        "prefix": str,
        "args": str,
        "args2": str,
        "args3": str,
        "args4": str,
    },
    defaults={"args": "", "args2": "", "args3": "", "args4": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        set -o pipefail
        VCFS=({{inputs.vcfs}})
        shopt -s nullglob
        BAMS=({{inputs.bams}}); BAIS=({{inputs.bais}}); DISTS=({{inputs.sample_dists}})
        FA=({{inputs.fasta}}); FAI=({{inputs.fai}}); NORM=({{inputs.chr_norm}})
        shopt -u nullglob
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        ARGS2={{inputs.args2}}
        ARGS3={{inputs.args3}}
        ARGS4={{inputs.args4}}
        W=$(mktemp -d); cd "$W"
        for F in "${{BAMS[@]}}" "${{BAIS[@]}}" "${{DISTS[@]}}" "${{FA[@]}}" "${{FAI[@]}}" "${{NORM[@]}}"; do
            ln -s "$F" "$(basename "$F")"
        done
        for F in "${{VCFS[@]}}"; do
            NAME=$(basename "$F")
            if [[ "$NAME" == *.vcf.gz ]]; then
                bgzip -d -c --threads $TASK_CPUS $ARGS4 "$F" > "${{NAME%.gz}}" || exit 1
            else
                ln -s "$F" "$NAME"
            fi
        done
        FIRST_VCF=$(basename "${{VCFS[0]}}"); FIRST_VCF=${{FIRST_VCF%.gz}}

        ls *.vcf > vcfs.txt
        EXTRA=()
        if [ ${{#BAMS[@]}} -gt 0 ]; then ls *.bam > bams.txt; EXTRA+=(bam_list=bams.txt); fi
        if [ -n "$ARGS2" ]; then EXTRA+=("iris_args=$ARGS2"); fi
        if [ ${{#DISTS[@]}} -gt 0 ]; then EXTRA+=("sample_dists=$(basename "${{DISTS[0]}}")"); fi
        if [ ${{#NORM[@]}} -gt 0 ]; then EXTRA+=("chr_norm_file=$(basename "${{NORM[0]}}")"); fi
        GENOME=""
        if [ ${{#FA[@]}} -gt 0 ]; then GENOME=$(basename "${{FA[0]}}"); fi

        jasmine \
            file_list=vcfs.txt \
            out_file="$P.vcf" \
            threads=$TASK_CPUS \
            genome_file="$GENOME" \
            "${{EXTRA[@]}}" \
            $ARGS || exit 1

        if [ -s "$P.vcf" ]; then
            echo "The file is not empty"
            bgzip --threads $TASK_CPUS $ARGS3 "$P.vcf" || exit 1
        else
            echo "The file is empty, using the header of the first VCF as output file"
            grep "#" "$FIRST_VCF" | bgzip --threads $TASK_CPUS $ARGS3 > "$P.vcf.gz" || exit 1
        fi
        mv "$P.vcf.gz" {{outputs.results}}/
    """,
)

env = flyte.TaskEnvironment.from_task("jasminesv", jasminesv_cmd.as_task())


async def jasminesv(
    vcfs: list[File],
    prefix: str,
    bams: list[File] | None = None,
    bais: list[File] | None = None,
    sample_dists: File | None = None,
    fasta: File | None = None,
    fai: File | None = None,
    chr_norm: File | None = None,
    args: str = "",
    args2: str = "",
    args3: str = "",
    args4: str = "",
) -> File:
    """Merge the structural variants in ``vcfs`` into ``<prefix>.vcf.gz``."""
    if not vcfs:
        raise ValueError("jasminesv: needs at least one VCF")
    for vcf in vcfs:
        name = vcf.name or vcf.path.rsplit("/", 1)[-1]
        if name.startswith(f"{prefix}.vcf"):
            raise ValueError("jasminesv: input and output names are the same; pass another prefix")
    results = await jasminesv_cmd(
        vcfs=vcfs,
        bams=bams or [],
        bais=bais or [],
        sample_dists=sample_dists,
        fasta=fasta,
        fai=fai,
        chr_norm=chr_norm,
        prefix=prefix,
        args=args,
        args2=args2,
        args3=args3,
        args4=args4,
    )
    vcf = await results.get_file(f"{prefix}.vcf.gz")
    if vcf is None:
        raise FileNotFoundError(f"jasmine wrote no {prefix}.vcf.gz")
    return vcf
