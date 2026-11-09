"""manta — structural variant and indel calling with Manta.

Exposes :func:`manta_somatic`, which calls somatic and germline structural
variants (and candidate small indels, for Strelka) from a tumour/normal pair
of BAMs or CRAMs, optionally within the regions of a bgzipped, indexed BED
(``--callRegions``) and with a Manta configuration file (``--config``).

As upstream, Manta's four VCFs (each with its ``.tbi``) are renamed:

- ``<prefix>.candidate_small_indels.vcf.gz`` — candidate small indels, for
  Strelka's ``--indelCandidates``;
- ``<prefix>.candidate_sv.vcf.gz`` — every candidate SV;
- ``<prefix>.diploid_sv.vcf.gz`` — germline SVs of the normal;
- ``<prefix>.somatic_sv.vcf.gz`` — somatic SVs of the tumour.

Alignments, the reference and the BED are linked beside their indices under
their own names. Manta runs locally with ``-j`` of the task's CPUs.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

MANTA_IMAGE = "community.wave.seqera.io/library/manta_python:0eb71149179b3920"
# Upstream's process_medium asks for 6 CPUs and 36 GB; the test data needs far less.
DEFAULT_RESOURCES = flyte.Resources(cpu=4, memory="8Gi")

# Manta's output VCFs (under manta/results/variants) and the names upstream gives them.
OUTPUTS = {
    "candidateSmallIndels": "candidate_small_indels",
    "candidateSV": "candidate_sv",
    "diploidSV": "diploid_sv",
    "somaticSV": "somatic_sv",
}
RENAMES = "\n".join(
    f'        mv manta/results/variants/{name}.vcf{ext} "$OUT/$P.{renamed}.vcf{ext}"'
    for name, renamed in OUTPUTS.items()
    for ext in (".gz", ".gz.tbi")
)

manta_somatic_cmd = shell.create(
    name="manta_somatic",
    image=MANTA_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "normal": File,
        "normal_index": File,
        "tumor": File,
        "tumor_index": File,
        "target_bed": File | None,
        "target_bed_tbi": File | None,
        "fasta": File,
        "fai": File,
        "config": File | None,
        "prefix": str,
        "args": str,
    },
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        N=({{inputs.normal}}); NI=({{inputs.normal_index}}); T=({{inputs.tumor}}); TI=({{inputs.tumor_index}})
        FA=({{inputs.fasta}}); FAI=({{inputs.fai}})
        shopt -s nullglob; BED=({{inputs.target_bed}}); BED_TBI=({{inputs.target_bed_tbi}}); CONFIG=({{inputs.config}})
        shopt -u nullglob
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        for f in "${{N[0]}}" "${{NI[0]}}" "${{T[0]}}" "${{TI[0]}}" "${{FA[0]}}"; do ln -s "$f" "$(basename "$f")"; done
        ln -s "${{FAI[0]}}" "$(basename "${{FA[0]}}").fai"
        OPTIONAL=()
        if [ ${{#CONFIG[@]}} -gt 0 ]; then OPTIONAL+=(--config "${{CONFIG[0]}}"); fi
        if [ ${{#BED[@]}} -gt 0 ]; then
            ln -s "${{BED[0]}}" "$(basename "${{BED[0]}}")"
            ln -s "${{BED_TBI[0]}}" "$(basename "${{BED[0]}}").tbi"
            OPTIONAL+=(--callRegions "$(basename "${{BED[0]}}")")
        fi
        configManta.py \
            --tumorBam "$(basename "${{T[0]}}")" \
            --normalBam "$(basename "${{N[0]}}")" \
            --reference "$(basename "${{FA[0]}}")" \
            --runDir manta \
            "${{OPTIONAL[@]}}" \
            $ARGS || exit 1
        python manta/runWorkflow.py -m local -j $TASK_CPUS || exit 1
{RENAMES}
    """,
)

env = flyte.TaskEnvironment.from_task("manta", manta_somatic_cmd.as_task())


@dataclass
class MantaSomaticResult:
    candidate_small_indels_vcf: File  # <prefix>.candidate_small_indels.vcf.gz
    candidate_small_indels_vcf_tbi: File
    candidate_sv_vcf: File  # <prefix>.candidate_sv.vcf.gz
    candidate_sv_vcf_tbi: File
    diploid_sv_vcf: File  # <prefix>.diploid_sv.vcf.gz
    diploid_sv_vcf_tbi: File
    somatic_sv_vcf: File  # <prefix>.somatic_sv.vcf.gz
    somatic_sv_vcf_tbi: File


async def manta_somatic(
    normal: File,
    normal_index: File,
    tumor: File,
    tumor_index: File,
    fasta: File,
    fai: File,
    prefix: str,
    target_bed: File | None = None,
    target_bed_tbi: File | None = None,
    config: File | None = None,
    args: str = "",
) -> MantaSomaticResult:
    """Somatic and germline SVs of a tumour/normal pair (BAM or CRAM, with their indices).

    ``target_bed`` is a bgzipped BED of the regions to call, with its
    ``target_bed_tbi``; ``config`` is a Manta configuration (``.ini``) file.
    """
    if (target_bed is None) != (target_bed_tbi is None):
        raise ValueError("target_bed and target_bed_tbi go together")
    results = await manta_somatic_cmd(
        normal=normal,
        normal_index=normal_index,
        tumor=tumor,
        tumor_index=tumor_index,
        target_bed=target_bed,
        target_bed_tbi=target_bed_tbi,
        fasta=fasta,
        fai=fai,
        config=config,
        prefix=prefix,
        args=args,
    )

    async def pick(name: str) -> File:
        found = await results.get_file(name)
        if found is None:
            raise FileNotFoundError(f"Manta wrote no {name}")
        return found

    files = {}
    for renamed in OUTPUTS.values():
        files[f"{renamed}_vcf"] = await pick(f"{prefix}.{renamed}.vcf.gz")
        files[f"{renamed}_vcf_tbi"] = await pick(f"{prefix}.{renamed}.vcf.gz.tbi")
    return MantaSomaticResult(**files)
