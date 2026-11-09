"""sawfish — structural-variant and copy-number calling from PacBio HiFi reads with sawfish.

Exposes:

- :func:`sawfish_discover` — find one sample's SV candidates (and, from its
  depth, copy-number segments) from its aligned BAM (with its index)
  against the reference, optionally with expected copy numbers
  (``expected_cn_bed``), the sample's small variants for allele
  frequencies (``maf_vcf``) and regions to leave out of CNV calling
  (``cnv_exclude_regions_bed``). The result is sawfish's discover
  directory: ``assembly.regions.bed``, ``candidate.sv.bcf`` (+ ``.csi``),
  ``contig.alignment.bam`` (+ ``.csi``), ``depth.mpack``,
  ``discover.settings.json``, ``max.depth.bed``, ``run.stats.json``,
  ``sawfish.log``, the ``debug.*`` files and, when written,
  ``copynum.{bedgraph,mpack}``, ``genome.gclevels.mpack``,
  ``sample.gcbias.mpack``, ``maf.mpack`` and ``expected.copy.number.bed``.
- :func:`sawfish_jointcall` — genotype SVs and CNVs jointly across samples'
  discover directories: ``<prefix>_genotyped.sv.vcf.gz`` (+ ``.tbi``),
  ``contig.alignment.bam`` (+ ``.csi``), ``run.stats.json``,
  ``sawfish.log`` and per sample ``samples/<sample>/depth.bw`` and, when
  written, ``copynum.bedgraph``, ``gc_bias_corrected_depth.bw``,
  ``copynum.summary.json`` and ``maf.bw``.

As upstream, discover runs with ``--disable-path-canonicalization``, so its
directory records the BAM and reference by the relative names they're
linked under; joint calling reads the BAMs again by those names, so it
takes the same BAMs. Discover's output directory is named after its
``prefix``, which its contig BAM's header records; each discover directory
keeps that name for joint calling, which passes them (sorted by name) as ``--sample``, or,
with ``sample_csv``, reads them from that CSV instead.

``--threads`` follows the task's CPUs. ``args`` go to sawfish and aren't
shell-quoted.
"""

import asyncio
import tempfile
from dataclasses import dataclass
from pathlib import Path

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

SAWFISH_IMAGE = "community.wave.seqera.io/library/sawfish:2.3.1--1e964dca544b9188"
# Upstream labels both process_medium; its tests run with 4 CPUs, which joint calling records in its VCF header.
DEFAULT_RESOURCES = flyte.Resources(cpu=4, memory="8Gi")

sawfish_discover_cmd = shell.create(
    name="sawfish_discover",
    image=SAWFISH_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "bam": File,
        "bai": File,
        "fasta": File,
        "expected_cn_bed": File | None,
        "maf_vcf": File | None,
        "cnv_exclude_regions_bed": File | None,
        "prefix": str,
        "args": str,
    },
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        BAM=({{inputs.bam}}); BAI=({{inputs.bai}}); FA=({{inputs.fasta}})
        shopt -s nullglob
        EXPECTED=({{inputs.expected_cn_bed}}); MAF=({{inputs.maf_vcf}}); EXCLUDE=({{inputs.cnv_exclude_regions_bed}})
        shopt -u nullglob
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        W=$(mktemp -d); cd "$W"
        BAM_NAME=$(basename "${{BAM[0]}}"); ln -s "${{BAM[0]}}" "$BAM_NAME"
        ln -s "${{BAI[0]}}" "$BAM_NAME.${{BAI[0]##*.}}"
        FA_NAME=$(basename "${{FA[0]}}"); ln -s "${{FA[0]}}" "$FA_NAME"
        for F in "${{EXPECTED[@]}}" "${{MAF[@]}}" "${{EXCLUDE[@]}}"; do ln -s "$F" "$(basename "$F")"; done
        EXTRA=()
        if [ ${{#EXPECTED[@]}} -gt 0 ]; then EXTRA+=(--expected-cn "$(basename "${{EXPECTED[0]}}")"); fi
        if [ ${{#EXCLUDE[@]}} -gt 0 ]; then
            EXTRA+=(--cnv-excluded-regions "$(basename "${{EXCLUDE[0]}}")")
        fi
        if [ ${{#MAF[@]}} -gt 0 ]; then EXTRA+=(--maf "$(basename "${{MAF[0]}}")"); fi

        sawfish discover \
            --disable-path-canonicalization \
            --ref "$FA_NAME" \
            --bam "$BAM_NAME" \
            --threads $TASK_CPUS \
            $ARGS \
            "${{EXTRA[@]}}" \
            --output-dir "$P" || exit 1

        mv "$P"/* {{outputs.results}}/
    """,
)

# `samples` holds one discover directory per sample, named as joint calling should name it.
sawfish_jointcall_cmd = shell.create(
    name="sawfish_jointcall",
    image=SAWFISH_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "samples": Dir,
        "fasta": File,
        "bams": list[File],
        "bais": list[File],
        "sample_csv": File | None,
        "prefix": str,
        "args": str,
    },
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        FA=({{inputs.fasta}}); BAMS=({{inputs.bams}}); BAIS=({{inputs.bais}})
        shopt -s nullglob; CSV=({{inputs.sample_csv}}); shopt -u nullglob
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        W=$(mktemp -d); cd "$W"
        FA_NAME=$(basename "${{FA[0]}}"); ln -s "${{FA[0]}}" "$FA_NAME"
        for I in "${{!BAMS[@]}}"; do
            NAME=$(basename "${{BAMS[$I]}}"); ln -s "${{BAMS[$I]}}" "$NAME"
            ln -s "${{BAIS[$I]}}" "$NAME.${{BAIS[$I]##*.}}"
        done
        SAMPLES=()
        for D in $(cd {{inputs.samples}} && LC_ALL=C ls -d */); do
            D=${{D%/}}; ln -s {{inputs.samples}}/"$D" "$D"; SAMPLES+=(--sample "$D")
        done
        if [ ${{#CSV[@]}} -gt 0 ]; then
            CSV_NAME=$(basename "${{CSV[0]}}"); ln -s "${{CSV[0]}}" "$CSV_NAME"
            SAMPLES=(--sample-csv "$CSV_NAME")
        fi

        sawfish joint-call \
            --threads $TASK_CPUS \
            --ref "$FA_NAME" \
            "${{SAMPLES[@]}}" \
            $ARGS \
            --output-dir "$P" || exit 1

        # As upstream, the VCF is renamed to carry the prefix.
        if [ -f "$P/genotyped.sv.vcf.gz" ]; then mv "$P/genotyped.sv.vcf.gz" "$P/${{P}}_genotyped.sv.vcf.gz"; fi
        if [ -f "$P/genotyped.sv.vcf.gz.tbi" ]; then
            mv "$P/genotyped.sv.vcf.gz.tbi" "$P/${{P}}_genotyped.sv.vcf.gz.tbi"
        fi
        mv "$P"/* {{outputs.results}}/
    """,
)

env = flyte.TaskEnvironment.from_task("sawfish", sawfish_discover_cmd.as_task(), sawfish_jointcall_cmd.as_task())


@dataclass
class SawfishDiscoverResult:
    name: str  # the sample's name for joint calling: discover's prefix
    candidate_sv_bcf: File  # candidate.sv.bcf
    candidate_sv_bcf_csi: File  # candidate.sv.bcf.csi
    contig_alignment_bam: File  # contig.alignment.bam
    contig_alignment_bam_csi: File  # contig.alignment.bam.csi
    depth_mpack: File  # depth.mpack
    copynum_bedgraph: File | None  # copynum.bedgraph
    log: File  # sawfish.log
    discover_dir: Dir  # the whole discover directory


async def sawfish_discover(
    bam: File,
    bai: File,
    fasta: File,
    expected_cn_bed: File | None = None,
    maf_vcf: File | None = None,
    cnv_exclude_regions_bed: File | None = None,
    prefix: str = "sample",
    args: str = "",
) -> SawfishDiscoverResult:
    """Discover SV candidates and copy-number segments in ``bam`` (with its ``bai``) against ``fasta``."""
    results = await sawfish_discover_cmd(
        bam=bam,
        bai=bai,
        fasta=fasta,
        expected_cn_bed=expected_cn_bed,
        maf_vcf=maf_vcf,
        cnv_exclude_regions_bed=cnv_exclude_regions_bed,
        prefix=prefix,
        args=args,
    )

    async def required(name: str) -> File:
        found = await results.get_file(name)
        if found is None:
            raise FileNotFoundError(f"sawfish discover wrote no {name}")
        return found

    return SawfishDiscoverResult(
        name=prefix,
        candidate_sv_bcf=await required("candidate.sv.bcf"),
        candidate_sv_bcf_csi=await required("candidate.sv.bcf.csi"),
        contig_alignment_bam=await required("contig.alignment.bam"),
        contig_alignment_bam_csi=await required("contig.alignment.bam.csi"),
        depth_mpack=await required("depth.mpack"),
        copynum_bedgraph=await results.get_file("copynum.bedgraph"),
        log=await required("sawfish.log"),
        discover_dir=results,
    )


async def sawfish_samples(discovered: list[SawfishDiscoverResult]) -> Dir:
    """The discover directories of ``discovered`` in one Dir, each named after its sample.

    Shell tasks can't take a ``list[Dir]``, so they're copied through the caller's pod.
    """
    names = [d.name for d in discovered]
    if len(set(names)) != len(names):
        raise ValueError(f"sample names aren't unique: {names}")
    root = Path(tempfile.mkdtemp(prefix="sawfish_samples_"))
    await asyncio.gather(*(d.discover_dir.download(str(root / d.name)) for d in discovered))
    return await Dir.from_local(str(root))


@dataclass
class SawfishJointcallResult:
    vcf: File  # <prefix>_genotyped.sv.vcf.gz
    tbi: File  # its .tbi
    bam: File  # contig.alignment.bam
    bam_index: File  # contig.alignment.bam.csi
    stats: File  # run.stats.json
    log: File  # sawfish.log
    results: Dir  # everything joint calling wrote, with samples/<sample>/ tracks


async def sawfish_jointcall(
    discovered: list[SawfishDiscoverResult],
    fasta: File,
    bams: list[File],
    bais: list[File],
    sample_csv: File | None = None,
    prefix: str = "cohort",
    args: str = "",
) -> SawfishJointcallResult:
    """Jointly genotype the samples in ``discovered``, reading their ``bams`` (the BAMs discover read)."""
    samples = await sawfish_samples(discovered)
    results = await sawfish_jointcall_cmd(
        samples=samples, fasta=fasta, bams=bams, bais=bais, sample_csv=sample_csv, prefix=prefix, args=args
    )

    async def required(name: str) -> File:
        found = await results.get_file(name)
        if found is None:
            raise FileNotFoundError(f"sawfish joint-call wrote no {name}")
        return found

    return SawfishJointcallResult(
        vcf=await required(f"{prefix}_genotyped.sv.vcf.gz"),
        tbi=await required(f"{prefix}_genotyped.sv.vcf.gz.tbi"),
        bam=await required("contig.alignment.bam"),
        bam_index=await required("contig.alignment.bam.csi"),
        stats=await required("run.stats.json"),
        log=await required("sawfish.log"),
        results=results,
    )
