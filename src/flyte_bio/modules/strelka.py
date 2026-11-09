"""strelka — small somatic variant calling with Strelka2.

Exposes :func:`strelka_somatic`, which calls somatic SNVs and indels from a
tumour/normal pair of BAMs or CRAMs, optionally within the regions of a
bgzipped, indexed BED (``--callRegions``) and with Manta's candidate small
indels (``--indelCandidates``, from :func:`flyte_bio.modules.manta.manta_somatic`),
as upstream's best practice recommends.

As upstream, the calls are renamed ``<prefix>.somatic_indels.vcf.gz`` and
``<prefix>.somatic_snvs.vcf.gz`` (each with its ``.tbi``), the workflow's
e-mail notification is switched off, and Strelka runs locally with ``-j`` of
the task's CPUs. Alignments, the reference, the candidates and the BED are
linked beside their indices under their own names.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

STRELKA_IMAGE = "quay.io/biocontainers/strelka:2.9.10--h9ee0642_1"
# Upstream's process_medium asks for 6 CPUs and 36 GB; the test data needs far less.
DEFAULT_RESOURCES = flyte.Resources(cpu=4, memory="8Gi")

strelka_somatic_cmd = shell.create(
    name="strelka_somatic",
    image=STRELKA_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "normal": File,
        "normal_index": File,
        "tumor": File,
        "tumor_index": File,
        "indel_candidates": File | None,
        "indel_candidates_tbi": File | None,
        "target_bed": File | None,
        "target_bed_tbi": File | None,
        "fasta": File,
        "fai": File,
        "prefix": str,
        "args": str,
    },
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        N=({{inputs.normal}}); NI=({{inputs.normal_index}}); T=({{inputs.tumor}}); TI=({{inputs.tumor_index}})
        FA=({{inputs.fasta}}); FAI=({{inputs.fai}})
        shopt -s nullglob
        CANDIDATES=({{inputs.indel_candidates}}); CANDIDATES_TBI=({{inputs.indel_candidates_tbi}})
        BED=({{inputs.target_bed}}); BED_TBI=({{inputs.target_bed_tbi}})
        shopt -u nullglob
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        for f in "${{N[0]}}" "${{NI[0]}}" "${{T[0]}}" "${{TI[0]}}" "${{FA[0]}}"; do ln -s "$f" "$(basename "$f")"; done
        ln -s "${{FAI[0]}}" "$(basename "${{FA[0]}}").fai"
        OPTIONAL=()
        if [ ${{#BED[@]}} -gt 0 ]; then
            ln -s "${{BED[0]}}" "$(basename "${{BED[0]}}")"
            ln -s "${{BED_TBI[0]}}" "$(basename "${{BED[0]}}").tbi"
            OPTIONAL+=(--callRegions "$(basename "${{BED[0]}}")")
        fi
        if [ ${{#CANDIDATES[@]}} -gt 0 ]; then
            ln -s "${{CANDIDATES[0]}}" "$(basename "${{CANDIDATES[0]}}")"
            ln -s "${{CANDIDATES_TBI[0]}}" "$(basename "${{CANDIDATES[0]}}").tbi"
            OPTIONAL+=(--indelCandidates "$(basename "${{CANDIDATES[0]}}")")
        fi
        configureStrelkaSomaticWorkflow.py \
            --tumor "$(basename "${{T[0]}}")" \
            --normal "$(basename "${{N[0]}}")" \
            --referenceFasta "$(basename "${{FA[0]}}")" \
            "${{OPTIONAL[@]}}" \
            $ARGS \
            --runDir strelka || exit 1
        sed -i s/"isEmail = isLocalSmtp()"/"isEmail = False"/g strelka/runWorkflow.py
        python strelka/runWorkflow.py -m local -j $TASK_CPUS || exit 1
        mv strelka/results/variants/somatic.indels.vcf.gz "$OUT/$P.somatic_indels.vcf.gz"
        mv strelka/results/variants/somatic.indels.vcf.gz.tbi "$OUT/$P.somatic_indels.vcf.gz.tbi"
        mv strelka/results/variants/somatic.snvs.vcf.gz "$OUT/$P.somatic_snvs.vcf.gz"
        mv strelka/results/variants/somatic.snvs.vcf.gz.tbi "$OUT/$P.somatic_snvs.vcf.gz.tbi"
    """,
)

env = flyte.TaskEnvironment.from_task("strelka", strelka_somatic_cmd.as_task())


@dataclass
class StrelkaSomaticResult:
    vcf_indels: File  # <prefix>.somatic_indels.vcf.gz
    vcf_indels_tbi: File
    vcf_snvs: File  # <prefix>.somatic_snvs.vcf.gz
    vcf_snvs_tbi: File


async def strelka_somatic(
    normal: File,
    normal_index: File,
    tumor: File,
    tumor_index: File,
    fasta: File,
    fai: File,
    prefix: str,
    indel_candidates: File | None = None,
    indel_candidates_tbi: File | None = None,
    target_bed: File | None = None,
    target_bed_tbi: File | None = None,
    args: str = "",
) -> StrelkaSomaticResult:
    """Somatic SNVs and indels of a tumour/normal pair (BAM or CRAM, with their indices).

    ``indel_candidates`` is a bgzipped VCF of candidate small indels (Manta's
    ``candidate_small_indels_vcf``), ``target_bed`` a bgzipped BED of the
    regions to call; each comes with its ``.tbi``. Pass ``--exome`` (or
    ``--targeted``) in ``args`` for targeted data.
    """
    if (indel_candidates is None) != (indel_candidates_tbi is None):
        raise ValueError("indel_candidates and indel_candidates_tbi go together")
    if (target_bed is None) != (target_bed_tbi is None):
        raise ValueError("target_bed and target_bed_tbi go together")
    results = await strelka_somatic_cmd(
        normal=normal,
        normal_index=normal_index,
        tumor=tumor,
        tumor_index=tumor_index,
        indel_candidates=indel_candidates,
        indel_candidates_tbi=indel_candidates_tbi,
        target_bed=target_bed,
        target_bed_tbi=target_bed_tbi,
        fasta=fasta,
        fai=fai,
        prefix=prefix,
        args=args,
    )

    async def pick(name: str) -> File:
        found = await results.get_file(name)
        if found is None:
            raise FileNotFoundError(f"Strelka wrote no {name}")
        return found

    return StrelkaSomaticResult(
        vcf_indels=await pick(f"{prefix}.somatic_indels.vcf.gz"),
        vcf_indels_tbi=await pick(f"{prefix}.somatic_indels.vcf.gz.tbi"),
        vcf_snvs=await pick(f"{prefix}.somatic_snvs.vcf.gz"),
        vcf_snvs_tbi=await pick(f"{prefix}.somatic_snvs.vcf.gz.tbi"),
    )
