"""hificnv — copy-number variant calling from PacBio HiFi reads with HiFiCNV.

Exposes :func:`hificnv`, which calls CNVs from one sample's aligned HiFi BAM
(with its index) against the reference, optionally with the sample's small
variants (``maf``, a VCF whose allele frequencies HiFiCNV tracks), regions to
exclude (``exclude``, a BED) and expected copy numbers (``expected_cn``, a
BED). HiFiCNV names its outputs ``<prefix>.<sample>.*``, the sample coming
from the BAM's read group:

- ``<prefix>.<sample>.copynum.bedgraph`` — the copy-number track;
- ``<prefix>.<sample>.depth.bw`` — the depth track;
- ``<prefix>.<sample>.maf.bw`` — the allele-frequency track, with ``maf``;
- ``<prefix>.<sample>.vcf.gz`` — the CNV calls.

``--threads`` follows the task's CPUs. ``args`` go to HiFiCNV (e.g.
``--cov-regex .`` to count depth on every contig); they aren't shell-quoted.
The inputs are linked into the working directory under their own names, so
HiFiCNV can index the reference beside it.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

HIFICNV_IMAGE = "community.wave.seqera.io/library/hificnv:1.0.1--b7e433ac6789e2d2"
# Upstream's label is process_medium; its tests run with 4 CPUs.
DEFAULT_RESOURCES = flyte.Resources(cpu=4, memory="8Gi")

hificnv_cmd = shell.create(
    name="hificnv",
    image=HIFICNV_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "bam": File,
        "bai": File,
        "maf": File | None,
        "fasta": File,
        "exclude": File | None,
        "expected_cn": File | None,
        "prefix": str,
        "args": str,
    },
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        BAM=({{inputs.bam}}); BAI=({{inputs.bai}}); FA=({{inputs.fasta}})
        shopt -s nullglob; MAF=({{inputs.maf}}); EXCLUDE=({{inputs.exclude}}); EXPECTED=({{inputs.expected_cn}})
        shopt -u nullglob
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        W=$(mktemp -d); cd "$W"
        BAM_NAME=$(basename "${{BAM[0]}}"); ln -s "${{BAM[0]}}" "$BAM_NAME"
        ln -s "${{BAI[0]}}" "$BAM_NAME.${{BAI[0]##*.}}"
        FA_NAME=$(basename "${{FA[0]}}"); ln -s "${{FA[0]}}" "$FA_NAME"
        EXTRA=()
        for F in "${{MAF[@]}}" "${{EXCLUDE[@]}}" "${{EXPECTED[@]}}"; do ln -s "$F" "$(basename "$F")"; done
        if [ ${{#MAF[@]}} -gt 0 ]; then EXTRA+=(--maf "$(basename "${{MAF[0]}}")"); fi
        if [ ${{#EXCLUDE[@]}} -gt 0 ]; then EXTRA+=(--exclude "$(basename "${{EXCLUDE[0]}}")"); fi
        if [ ${{#EXPECTED[@]}} -gt 0 ]; then EXTRA+=(--expected-cn "$(basename "${{EXPECTED[0]}}")"); fi

        hificnv \
            --bam "$BAM_NAME" \
            --ref "$FA_NAME" \
            "${{EXTRA[@]}}" \
            --threads $TASK_CPUS \
            --output-prefix "$P" \
            $ARGS || exit 1

        shopt -s nullglob
        mv "$P".*.copynum.bedgraph "$P".*.depth.bw "$P".*.maf.bw "$P".*.vcf.gz {{outputs.results}}/
    """,
)

env = flyte.TaskEnvironment.from_task("hificnv", hificnv_cmd.as_task())


@dataclass
class HificnvResult:
    copynum: File  # <prefix>.<sample>.copynum.bedgraph
    depth: File  # <prefix>.<sample>.depth.bw
    maf: File | None  # <prefix>.<sample>.maf.bw, with a maf VCF
    vcf: File  # <prefix>.<sample>.vcf.gz
    results: Dir  # everything above


async def hificnv(
    bam: File,
    bai: File,
    fasta: File,
    maf: File | None = None,
    exclude: File | None = None,
    expected_cn: File | None = None,
    prefix: str = "sample",
    args: str = "",
) -> HificnvResult:
    """Call copy-number variants in ``bam`` (with its ``bai``) against ``fasta``."""
    results = await hificnv_cmd(
        bam=bam, bai=bai, maf=maf, fasta=fasta, exclude=exclude, expected_cn=expected_cn, prefix=prefix, args=args
    )
    files = [f async for f in results.walk()]

    def ending(suffix: str, required: bool = True) -> File | None:
        found = [f for f in files if f.name.startswith(f"{prefix}.") and f.name.endswith(suffix)]
        if len(found) > 1 or (required and not found):
            raise FileNotFoundError(f"HiFiCNV wrote {len(found)} {prefix}.*{suffix} files, expected one")
        return found[0] if found else None

    copynum, depth, vcf = ending(".copynum.bedgraph"), ending(".depth.bw"), ending(".vcf.gz")
    assert copynum is not None and depth is not None and vcf is not None
    return HificnvResult(copynum=copynum, depth=depth, maf=ending(".maf.bw", required=False), vcf=vcf, results=results)
