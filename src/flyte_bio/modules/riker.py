"""riker — alignment QC metrics from a BAM or CRAM, many collectors in one pass.

Exposes :func:`riker_multi`, which runs ``riker multi`` (the collectors
chosen with ``--tools`` in ``args``) and returns a Dir of
``<prefix>.<metric>`` files, under upstream's names:

- ``alignment-metrics.txt``, ``base-distribution-by-cycle.txt``,
  ``mean-quality-by-cycle.txt``, ``quality-score-distribution.txt``;
- ``isize-metrics.txt``, ``isize-histogram.txt``;
- ``wgs-metrics.txt``, ``wgs-coverage.txt``, ``wgs-per-base.bw``;
- ``gcbias-detail.txt``, ``gcbias-summary.txt``;
- ``hybcap-metrics.txt``, ``hybcap-per-target.txt``, ``hybcap-per-base.txt*``;
- ``error-mismatch.txt``, ``error-overlap.txt``, ``error-indel.txt``;
- ``rna-metrics.txt``, ``rna-biotype.txt``, ``rna-insert-size.txt``,
  ``rna-insert-size-histogram.txt``;
- charts as ``*.pdf``.

Optional inputs feed the collectors that use them (a reference for CRAM,
``wgs`` and ``gcbias``; baits and targets for ``hybcap``, which must be given
together; a gene model for ``rna``). As upstream, ``error_vcf`` is used only
together with its index.
"""

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

RIKER_IMAGE = "community.wave.seqera.io/library/riker:0.5.0--fc5310a356298f75"

# Upstream's process_medium, capped at the 4 CPUs its tests run with.
DEFAULT_RESOURCES = flyte.Resources(cpu=4, memory="8Gi")

# The BAM/CRAM, VCF and FASTA are linked next to their indexes in in/;
# riker writes <prefix>.* into out/, which is then moved to the output Dir.
riker_multi_cmd = shell.create(
    name="riker_multi",
    image=RIKER_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "alignment": File,
        "alignment_index": File | None,
        "fasta": File | None,
        "fai": File | None,
        "error_vcf": File | None,
        "error_vcf_index": File | None,
        "error_intervals": File | None,
        "gcbias_exclude_intervals": File | None,
        "hybcap_baits": File | None,
        "hybcap_targets": File | None,
        "rna_gene_model": File | None,
        "rna_ribosomal_intervals": File | None,
        "wgs_intervals": File | None,
        "prefix": str,
        "args": str,
    },
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        ALN=({{inputs.alignment}})
        shopt -s nullglob
        ALN_INDEX=({{inputs.alignment_index}}); FA=({{inputs.fasta}}); FAI=({{inputs.fai}})
        VCF=({{inputs.error_vcf}}); VCF_INDEX=({{inputs.error_vcf_index}})
        ERROR_INTERVALS=({{inputs.error_intervals}}); GC_EXCLUDE=({{inputs.gcbias_exclude_intervals}})
        BAITS=({{inputs.hybcap_baits}}); TARGETS=({{inputs.hybcap_targets}})
        GENE_MODEL=({{inputs.rna_gene_model}}); RIBOSOMAL=({{inputs.rna_ribosomal_intervals}})
        WGS_INTERVALS=({{inputs.wgs_intervals}})
        shopt -u nullglob
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        W=$(mktemp -d); mkdir "$W/in" "$W/out"
        for f in "${{ALN[@]}}" "${{ALN_INDEX[@]}}" "${{FA[@]}}" "${{FAI[@]}}" "${{VCF[@]}}" "${{VCF_INDEX[@]}}"; do
            ln -s "$f" "$W/in/"
        done
        OPTS=(-i "$W/in/$(basename "${{ALN[0]}}")")
        if [ ${{#FA[@]}} -gt 0 ]; then OPTS+=(--reference "$W/in/$(basename "${{FA[0]}}")"); fi
        if [ ${{#BAITS[@]}} -gt 0 ] && [ ${{#TARGETS[@]}} -gt 0 ]; then
            OPTS+=(--hybcap::baits "${{BAITS[0]}}" --hybcap::targets "${{TARGETS[0]}}")
        fi
        if [ ${{#VCF[@]}} -gt 0 ] && [ ${{#VCF_INDEX[@]}} -gt 0 ]; then
            OPTS+=(--error::vcf "$W/in/$(basename "${{VCF[0]}}")")
        fi
        if [ ${{#ERROR_INTERVALS[@]}} -gt 0 ]; then OPTS+=(--error::intervals "${{ERROR_INTERVALS[0]}}"); fi
        if [ ${{#GC_EXCLUDE[@]}} -gt 0 ]; then OPTS+=(--gcbias::exclude-intervals "${{GC_EXCLUDE[0]}}"); fi
        if [ ${{#GENE_MODEL[@]}} -gt 0 ]; then OPTS+=(--rna::gene-model "${{GENE_MODEL[0]}}"); fi
        if [ ${{#RIBOSOMAL[@]}} -gt 0 ]; then OPTS+=(--rna::ribosomal-intervals "${{RIBOSOMAL[0]}}"); fi
        if [ ${{#WGS_INTERVALS[@]}} -gt 0 ]; then OPTS+=(--wgs::intervals "${{WGS_INTERVALS[0]}}"); fi
        cd "$W/out"
        riker multi "${{OPTS[@]}}" -o "$P" --threads $TASK_CPUS $ARGS || exit 1
        mv ./* {{outputs.results}}/
    """,
)

env = flyte.TaskEnvironment.from_task("riker", riker_multi_cmd.as_task())


async def riker_multi(
    alignment: File,
    prefix: str,
    alignment_index: File | None = None,
    fasta: File | None = None,
    fai: File | None = None,
    error_vcf: File | None = None,
    error_vcf_index: File | None = None,
    error_intervals: File | None = None,
    gcbias_exclude_intervals: File | None = None,
    hybcap_baits: File | None = None,
    hybcap_targets: File | None = None,
    rna_gene_model: File | None = None,
    rna_ribosomal_intervals: File | None = None,
    wgs_intervals: File | None = None,
    args: str = "",
) -> Dir:
    """``riker multi`` QC metrics of ``alignment``: a Dir of ``<prefix>.<metric>`` files."""
    if (hybcap_baits is None) != (hybcap_targets is None):
        raise ValueError("riker multi: hybcap baits and targets must be given together, or neither")
    return await riker_multi_cmd(
        alignment=alignment,
        alignment_index=alignment_index,
        fasta=fasta,
        fai=fai,
        error_vcf=error_vcf,
        error_vcf_index=error_vcf_index,
        error_intervals=error_intervals,
        gcbias_exclude_intervals=gcbias_exclude_intervals,
        hybcap_baits=hybcap_baits,
        hybcap_targets=hybcap_targets,
        rna_gene_model=rna_gene_model,
        rna_ribosomal_intervals=rna_ribosomal_intervals,
        wgs_intervals=wgs_intervals,
        prefix=prefix,
        args=args,
    )
