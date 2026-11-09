"""pbcpgtools — CpG methylation scores from PacBio HiFi reads with pb-CpG-tools.

Exposes :func:`pbcpgtools_alignedbamtocpgscores`, which turns one sample's
aligned BAM (with its index), whose reads carry 5mC base modifications
(``MM``/``ML`` tags), into per-site CpG methylation scores:

- ``<prefix>.combined.bed.gz`` (with ``.tbi``) and ``<prefix>.combined.bw``
  — scores from all reads;
- ``<prefix>.hap1.*`` and ``<prefix>.hap2.*`` — the same per haplotype, when
  the reads are haplotagged (``HP`` tags).

``--threads`` follows the task's CPUs; the BED headers record the command
line, thread count included. ``args`` go to ``aligned_bam_to_cpg_scores``
and aren't shell-quoted.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

PBCPGTOOLS_IMAGE = "quay.io/biocontainers/pb-cpg-tools:3.0.0--h9ee0642_0"
# Upstream's label is process_medium; its tests run with 4 CPUs.
DEFAULT_RESOURCES = flyte.Resources(cpu=4, memory="8Gi")

pbcpgtools_alignedbamtocpgscores_cmd = shell.create(
    name="pbcpgtools_alignedbamtocpgscores",
    image=PBCPGTOOLS_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"bam": File, "bai": File, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        BAM=({{inputs.bam}}); BAI=({{inputs.bai}})
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        W=$(mktemp -d); cd "$W"
        BAM_NAME=$(basename "${{BAM[0]}}"); ln -s "${{BAM[0]}}" "$BAM_NAME"
        ln -s "${{BAI[0]}}" "$BAM_NAME.${{BAI[0]##*.}}"

        aligned_bam_to_cpg_scores \
            --bam "$BAM_NAME" \
            --output-prefix "$P" \
            --threads $TASK_CPUS \
            $ARGS || exit 1

        shopt -s nullglob
        mv "$P".combined.* "$P".hap1.* "$P".hap2.* {{outputs.results}}/
    """,
)

env = flyte.TaskEnvironment.from_task("pbcpgtools", pbcpgtools_alignedbamtocpgscores_cmd.as_task())


@dataclass
class CpgScores:
    bed: File  # <prefix>.<track>.bed.gz
    bed_index: File  # its .tbi
    bigwig: File  # <prefix>.<track>.bw


@dataclass
class PbcpgtoolsResult:
    combined: CpgScores  # from all reads
    hap1: CpgScores | None  # from haplotype 1's reads, when the reads are haplotagged
    hap2: CpgScores | None  # from haplotype 2's reads
    results: Dir  # everything above


async def pbcpgtools_alignedbamtocpgscores(
    bam: File, bai: File, prefix: str = "sample", args: str = ""
) -> PbcpgtoolsResult:
    """CpG methylation scores from the 5mC-tagged reads in ``bam`` (with its ``bai``)."""
    results = await pbcpgtools_alignedbamtocpgscores_cmd(bam=bam, bai=bai, prefix=prefix, args=args)
    return await pbcpgtools_scores(results, prefix)


async def pbcpgtools_scores(results: Dir, prefix: str) -> PbcpgtoolsResult:
    """Sort the files of an ``aligned_bam_to_cpg_scores`` run into a :class:`PbcpgtoolsResult`."""

    async def track(name: str) -> CpgScores | None:
        files = [await results.get_file(f"{prefix}.{name}.{ext}") for ext in ("bed.gz", "bed.gz.tbi", "bw")]
        if all(f is None for f in files):
            return None
        bed, bed_index, bigwig = files
        if bed is None or bed_index is None or bigwig is None:
            raise FileNotFoundError(f"aligned_bam_to_cpg_scores wrote only some of {prefix}.{name}.*")
        return CpgScores(bed=bed, bed_index=bed_index, bigwig=bigwig)

    combined = await track("combined")
    if combined is None:
        raise FileNotFoundError(f"aligned_bam_to_cpg_scores wrote no {prefix}.combined.*")
    return PbcpgtoolsResult(combined=combined, hap1=await track("hap1"), hap2=await track("hap2"), results=results)
