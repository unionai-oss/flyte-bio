"""sentieon — Sentieon's licensed, accelerated alignment, QC and variant calling.

**Untested: these tasks need a Sentieon licence, which we don't have.** They
are ported from upstream's modules without being run; their tests are
placeholders that report as not tested.

Exposes:

- :func:`sentieon_bwaindex` — a BWA index of a FASTA, a Dir of
  ``<prefix>.{amb,ann,bwt,pac,sa}`` (prefix: the FASTA's base name).
- :func:`sentieon_bwamem` — align reads with Sentieon's BWA and sort them into
  ``<prefix>`` (``.bam`` or ``.cram``) with its ``.bai`` / ``.crai``.
- :func:`sentieon_readwriter` — merge and/or convert alignments into
  ``<prefix>`` (``.bam`` or ``.cram``) with its index.
- :func:`sentieon_dedup` — mark (or remove) duplicates: ``<prefix>`` with its
  index, ``<prefix without extension>.score``, ``<prefix>.metrics`` and
  ``<prefix>.metrics.multiqc.tsv``.
- :func:`sentieon_datametrics` — alignment QC: ``<prefix>_{mq,qd,is,gc,aln}_metrics.txt``
  and ``<prefix>_gc_summary.txt``, plus PDF plots when asked.
- :func:`sentieon_wgsmetrics` — whole-genome coverage metrics, ``<prefix>.txt``.
- :func:`sentieon_qualcal` — a base-quality recalibration table,
  ``<prefix>.table`` (or, given a table, the post-recalibration table, CSV
  and PDF report, and optionally recalibrated alignments).
- :func:`sentieon_haplotyper` — germline calling: ``<prefix>.unfiltered.vcf.gz``
  and/or ``<prefix>.g.vcf.gz``, each with its ``.tbi``.
- :func:`sentieon_gvcftyper` — joint genotyping of GVCFs into ``<prefix>.vcf.gz``
  (+ ``.tbi``).
- :func:`sentieon_dnascope` — DNAscope germline calling, outputs as Haplotyper's.
- :func:`sentieon_dnamodelapply` — apply a DNAscope model to a VCF, ``<prefix>.vcf.gz``
  (+ ``.tbi``).
- :func:`sentieon_tnscope` — somatic calling, ``<prefix>.vcf.gz`` (+ ``.tbi``).
- :func:`sentieon_varcal` — VQSR model: ``<prefix>.recal`` (+ ``.idx``),
  ``<prefix>.tranches`` and an optional ``*plots.R``.
- :func:`sentieon_applyvarcal` — apply a VQSR model, ``<prefix>.vcf.gz`` (+ ``.tbi``).
- :func:`sentieon_staralign` — Sentieon's STAR: everything STAR writes under
  ``<prefix>.`` (logs, BAMs, junctions, signal tracks), with unmapped mates
  renamed to ``<prefix>.unmapped_{1,2}.fastq.gz``.
- :func:`sentieon_rsempreparereference` — an RSEM (optionally RSEM + STAR)
  reference, a Dir holding ``genome.*`` and the FASTA, and
  ``genome.transcripts.fa``.
- :func:`sentieon_rsemcalculateexpression` — RSEM quantification with
  Sentieon's STAR: ``<prefix>.{genes,isoforms}.results``, ``<prefix>.stat/``
  and optional BAMs.

Licence: as upstream, the tools read the licence from ``SENTIEON_LICENSE``
(a licence server's ``host:port``, or the path of a licence file). The tasks
get it from Flyte secrets, mounted as environment variables:

- ``SENTIEON_LICENSE`` — the licence server address (secret ``sentieon_license``
  by default, the only one the tasks request unless told otherwise);
- ``SENTIEON_LICENSE_BASE64`` — a base64-encoded licence file, which the
  script decodes into a temporary file and points ``SENTIEON_LICENSE`` at, as
  upstream does with its ``SENTIEON_LICENSE_BASE64`` secret;
- ``SENTIEON_AUTH_MECH`` and ``SENTIEON_AUTH_DATA`` — Sentieon's
  authentication mechanism and its data, for licence servers that ask for them.

A missing secret fails the task, so request only those you created: build the
list with :func:`sentieon_secrets` and pass it as any wrapper's ``secrets``.

Threads are the task's CPUs; Sentieon's BWA gets 90% of the task's memory.
Inputs are linked into the working directory under their own names, so
alignments, VCFs and FASTAs sit next to their indexes.
"""

from dataclasses import dataclass
from fnmatch import fnmatch

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

# Upstream builds the same Sentieon release (202503.03) into two images; each
# module keeps the one upstream names.
SENTIEON_IMAGE = "community.wave.seqera.io/library/sentieon:202503.03--df1987151f8b6d33"
SENTIEON_CALLING_IMAGE = "community.wave.seqera.io/library/sentieon:202503.03--5e34aa16344b911c"
# Sentieon with findutils (bwamem looks the index up with `find`).
SENTIEON_FINDUTILS_IMAGE = "community.wave.seqera.io/library/sentieon_findutils:59df946da6a128ad"
# Sentieon with gnuplot, for `sentieon plot`.
SENTIEON_GNUPLOT_IMAGE = "community.wave.seqera.io/library/sentieon_gnuplot:84874a80b33697db"
RSEM_SENTIEON_PREPARE_IMAGE = "community.wave.seqera.io/library/rsem_sentieon:2c2a3e1e4a56ec97"
RSEM_SENTIEON_CALCULATE_IMAGE = "community.wave.seqera.io/library/rsem_sentieon:3e4315fa0b636313"

# Upstream's process_low / process_medium / process_high, capped at what the cluster schedules.
LOW_RESOURCES = flyte.Resources(cpu=2, memory="12Gi")
MEDIUM_RESOURCES = flyte.Resources(cpu=6, memory="32Gi")
HIGH_RESOURCES = flyte.Resources(cpu=8, memory="32Gi")

LICENSE_SECRET = "sentieon_license"


def sentieon_secrets(
    license: str | None = LICENSE_SECRET,
    license_base64: str | None = None,
    auth_mech: str | None = None,
    auth_data: str | None = None,
) -> list[flyte.Secret]:
    """The Flyte secrets that carry the licence, by secret key; None leaves one out.

    ``license`` holds the licence server address (``SENTIEON_LICENSE``),
    ``license_base64`` a base64-encoded licence file (``SENTIEON_LICENSE_BASE64``),
    ``auth_mech`` and ``auth_data`` the server's authentication
    (``SENTIEON_AUTH_MECH``, ``SENTIEON_AUTH_DATA``).
    """
    wanted = {
        "SENTIEON_LICENSE": license,
        "SENTIEON_LICENSE_BASE64": license_base64,
        "SENTIEON_AUTH_MECH": auth_mech,
        "SENTIEON_AUTH_DATA": auth_data,
    }
    return [flyte.Secret(key=key, as_env_var=var) for var, key in wanted.items() if key is not None]


DEFAULT_SECRETS = sentieon_secrets()

# Every script starts with this: the task's limits, upstream's licence-file
# handling, a working directory, and `link`, which links files into it under
# their own names.
PRELUDE = (
    TASK_LIMITS
    + r"""
        if [ -n "${SENTIEON_LICENSE_BASE64:-}" ]; then
            export SENTIEON_LICENSE=$(mktemp); echo -e "$SENTIEON_LICENSE_BASE64" | base64 -d > "$SENTIEON_LICENSE"
        fi
        link() { for f in "$@"; do ln -sf "$f" "$(basename "$f")"; done; }
        W=$(mktemp -d); cd "$W"
"""
)


def sentieon_task(name: str, image: str, resources: flyte.Resources, script: str, **kwargs):
    return shell.create(
        name=name,
        image=image,
        resources=resources,
        secrets=DEFAULT_SECRETS,
        outputs={"results": Dir},
        script=PRELUDE + script,
        **kwargs,
    )


sentieon_bwaindex_cmd = sentieon_task(
    "sentieon_bwaindex",
    SENTIEON_IMAGE,
    HIGH_RESOURCES,
    inputs={"fasta": File, "prefix": str, "args": str},
    defaults={"prefix": "", "args": ""},
    script=r"""
        FA=({inputs.fasta}); link "${FA[0]}"
        NAME=$(basename "${FA[0]}")
        P={inputs.prefix}
        if [ -z "$P" ]; then P="${NAME%.*}"; fi
        ARGS={inputs.args}
        sentieon bwa index $ARGS -p {outputs.results}/"$P" "$NAME"
    """,
)

sentieon_bwamem_cmd = sentieon_task(
    "sentieon_bwamem",
    SENTIEON_FINDUTILS_IMAGE,
    HIGH_RESOURCES,
    inputs={
        "reads_1": File,
        "reads_2": File | None,
        "index": Dir,
        "fasta": File,
        "fai": File | None,
        "prefix": str,
        "args": str,
    },
    defaults={"args": ""},
    script=r"""
        R1=({inputs.reads_1}); FA=({inputs.fasta})
        shopt -s nullglob; R2=({inputs.reads_2}); FAI=({inputs.fai}); shopt -u nullglob
        link "${FA[0]}" "${FAI[@]}"
        P={inputs.prefix}
        ARGS={inputs.args}
        export bwt_max_mem="$(( TASK_MEMORY_MB * 9 / 10 / 1024 ))G"
        INDEX=$(find -L {inputs.index} -name "*.amb" | sed 's/.amb//')
        sentieon bwa mem $ARGS -t $TASK_CPUS "$INDEX" "${R1[0]}" "${R2[@]}" \
            | sentieon util sort -r "$(basename "${FA[0]}")" -t $TASK_CPUS -o "$P" --sam2bam -
        if [[ "$P" == *.cram ]]; then rm -f "$P.bai"; fi
        shopt -s nullglob
        mv "$P" "$P".* {outputs.results}/
    """,
)

sentieon_readwriter_cmd = sentieon_task(
    "sentieon_readwriter",
    SENTIEON_IMAGE,
    MEDIUM_RESOURCES,
    inputs={
        "alignments": list[File],
        "indexes": list[File],
        "fasta": File | None,
        "fai": File | None,
        "prefix": str,
        "args": str,
        "args2": str,
    },
    defaults={"args": "", "args2": ""},
    script=r"""
        ALN=({inputs.alignments}); IDX=({inputs.indexes})
        shopt -s nullglob; FA=({inputs.fasta}); FAI=({inputs.fai}); shopt -u nullglob
        link "${ALN[@]}" "${IDX[@]}" "${FA[@]}" "${FAI[@]}"
        P={inputs.prefix}
        ARGS={inputs.args}; ARGS2={inputs.args2}
        REF=(); if [ ${#FA[@]} -gt 0 ]; then REF=(-r "$(basename "${FA[0]}")"); fi
        INPUTS=(); for a in "${ALN[@]}"; do INPUTS+=(-i "$(basename "$a")"); done
        sentieon driver -t $TASK_CPUS "${REF[@]}" $ARGS "${INPUTS[@]}" --algo ReadWriter $ARGS2 "$P"
        shopt -s nullglob
        mv "$P" "$P".* {outputs.results}/
    """,
)

sentieon_dedup_cmd = sentieon_task(
    "sentieon_dedup",
    SENTIEON_IMAGE,
    MEDIUM_RESOURCES,
    inputs={
        "bams": list[File],
        "bais": list[File],
        "fasta": File,
        "fai": File,
        "prefix": str,
        "args": str,
        "args2": str,
        "args3": str,
        "args4": str,
    },
    defaults={"args": "", "args2": "", "args3": "", "args4": ""},
    script=r"""
        BAMS=({inputs.bams}); BAIS=({inputs.bais}); FA=({inputs.fasta}); FAI=({inputs.fai})
        link "${BAMS[@]}" "${BAIS[@]}" "${FA[0]}" "${FAI[0]}"
        REF=$(basename "${FA[0]}")
        P={inputs.prefix}
        ARGS={inputs.args}; ARGS2={inputs.args2}; ARGS3={inputs.args3}; ARGS4={inputs.args4}
        METRICS="$P.metrics"
        BASE="${P%.*}"
        INPUTS=(); for b in "${BAMS[@]}"; do INPUTS+=(-i "$(basename "$b")"); done
        sentieon driver $ARGS -t $TASK_CPUS "${INPUTS[@]}" -r "$REF" --algo LocusCollector $ARGS2 \
            --fun score_info "$BASE.score"
        sentieon driver $ARGS3 -t $TASK_CPUS "${INPUTS[@]}" -r "$REF" --algo Dedup $ARGS4 \
            --score_info "$BASE.score" --metrics "$METRICS" "$P"
        # As upstream: the metrics' first lines as a table MultiQC reads as custom content.
        head -3 "$METRICS" > "$METRICS.multiqc.tsv"
        mv "$P" "$P".* "$BASE.score" {outputs.results}/
    """,
)

sentieon_datametrics_cmd = sentieon_task(
    "sentieon_datametrics",
    SENTIEON_GNUPLOT_IMAGE,
    MEDIUM_RESOURCES,
    inputs={
        "bams": list[File],
        "bais": list[File],
        "fasta": File,
        "fai": File,
        "plot_results": bool,
        "prefix": str,
        "args": str,
    },
    defaults={"plot_results": False, "args": ""},
    script=r"""
        BAMS=({inputs.bams}); BAIS=({inputs.bais}); FA=({inputs.fasta}); FAI=({inputs.fai})
        link "${BAMS[@]}" "${BAIS[@]}" "${FA[0]}" "${FAI[0]}"
        P={inputs.prefix}
        ARGS={inputs.args}
        INPUTS=(); for b in "${BAMS[@]}"; do INPUTS+=(-i "$(basename "$b")"); done
        sentieon driver -t $TASK_CPUS -r "$(basename "${FA[0]}")" "${INPUTS[@]}" $ARGS \
            --algo GCBias --summary "${P}_gc_summary.txt" "${P}_gc_metrics.txt" \
            --algo MeanQualityByCycle "${P}_mq_metrics.txt" \
            --algo QualDistribution "${P}_qd_metrics.txt" \
            --algo InsertSizeMetricAlgo "${P}_is_metrics.txt" \
            --algo AlignmentStat "${P}_aln_metrics.txt"
        if {inputs.plot_results}; then
            sentieon plot GCBias -o "${P}_gc_metrics.pdf" "${P}_gc_metrics.txt"
            sentieon plot MeanQualityByCycle -o "${P}_mq_metrics.pdf" "${P}_mq_metrics.txt"
            sentieon plot QualDistribution -o "${P}_qd_metrics.pdf" "${P}_qd_metrics.txt"
            sentieon plot InsertSizeMetricAlgo -o "${P}_is_metrics.pdf" "${P}_is_metrics.txt"
        fi
        shopt -s nullglob
        mv ./*_metrics.txt ./*_gc_summary.txt ./*_metrics.pdf {outputs.results}/
    """,
)

sentieon_wgsmetrics_cmd = sentieon_task(
    "sentieon_wgsmetrics",
    SENTIEON_IMAGE,
    MEDIUM_RESOURCES,
    inputs={
        "bams": list[File],
        "bais": list[File],
        "fasta": File,
        "fai": File,
        "intervals": File | None,
        "prefix": str,
        "args": str,
    },
    defaults={"args": ""},
    script=r"""
        BAMS=({inputs.bams}); BAIS=({inputs.bais}); FA=({inputs.fasta}); FAI=({inputs.fai})
        shopt -s nullglob; IV=({inputs.intervals}); shopt -u nullglob
        link "${BAMS[@]}" "${BAIS[@]}" "${FA[0]}" "${FAI[0]}" "${IV[@]}"
        P={inputs.prefix}
        ARGS={inputs.args}
        INTERVAL=(); if [ ${#IV[@]} -gt 0 ]; then INTERVAL=(--interval "$(basename "${IV[0]}")"); fi
        INPUTS=(); for b in "${BAMS[@]}"; do INPUTS+=(-i "$(basename "$b")"); done
        sentieon driver -t $TASK_CPUS -r "$(basename "${FA[0]}")" "${INPUTS[@]}" "${INTERVAL[@]}" $ARGS \
            --algo WgsMetricsAlgo {outputs.results}/"$P.txt"
    """,
)

# Without a table: QualCal writes <prefix>.table. With one (and a single
# alignment): the post-recalibration table, the before/after CSV and its PDF
# plot, and with generate_recalibrated_bams the recalibrated alignment. As
# upstream, the plot step reads <prefix>.table, so the given table should be
# named that.
sentieon_qualcal_cmd = sentieon_task(
    "sentieon_qualcal",
    SENTIEON_GNUPLOT_IMAGE,
    MEDIUM_RESOURCES,
    inputs={
        "alignments": list[File],
        "indexes": list[File],
        "fasta": File,
        "fai": File,
        "known_sites": list[File],
        "known_sites_tbi": list[File],
        "recalibration_table": File | None,
        "generate_recalibrated_bams": bool,
        "prefix": str,
        "args": str,
    },
    defaults={"known_sites": [], "known_sites_tbi": [], "generate_recalibrated_bams": False, "args": ""},
    script=r"""
        ALN=({inputs.alignments}); IDX=({inputs.indexes}); FA=({inputs.fasta}); FAI=({inputs.fai})
        shopt -s nullglob
        SITES=({inputs.known_sites}); SITES_TBI=({inputs.known_sites_tbi}); TABLE=({inputs.recalibration_table})
        shopt -u nullglob
        link "${ALN[@]}" "${IDX[@]}" "${FA[0]}" "${FAI[0]}" "${SITES[@]}" "${SITES_TBI[@]}" "${TABLE[@]}"
        REF=$(basename "${FA[0]}")
        P={inputs.prefix}
        ARGS={inputs.args}
        KNOWN=(); for s in "${SITES[@]}"; do KNOWN+=(-k "$(basename "$s")"); done
        if [ ${#TABLE[@]} -eq 0 ]; then
            INPUTS=(); for a in "${ALN[@]}"; do INPUTS+=(-i "$(basename "$a")"); done
            sentieon driver -t $TASK_CPUS -r "$REF" "${INPUTS[@]}" --algo QualCal $ARGS "${KNOWN[@]}" "$P.table"
        else
            IN=$(basename "${ALN[0]}")
            if [[ "$IN" == *.bam ]]; then SUFFIX=bam; else SUFFIX=cram; fi
            RECAL=()
            if {inputs.generate_recalibrated_bams}; then RECAL=(--algo ReadWriter "$P.recalibrated.$SUFFIX"); fi
            sentieon driver -t $TASK_CPUS -r "$REF" -i "$IN" -q "$(basename "${TABLE[0]}")" \
                --algo QualCal $ARGS "${KNOWN[@]}" "$P.table.post" "${RECAL[@]}"
            sentieon driver -t $TASK_CPUS --algo QualCal --plot --before "$P.table" --after "$P.table.post" \
                $ARGS "$P.csv"
            sentieon plot QualCal -o "$P.pdf" "$P.csv"
        fi
        shopt -s nullglob
        for f in "$P".table "$P".table.post "$P".recalibrated.* "$P".csv "$P".pdf; do
            if [ -e "$f" ] && [ ! -L "$f" ]; then mv "$f" {outputs.results}/; fi
        done
    """,
)

# Haplotyper and DNAscope can write a VCF and a GVCF in one pass over the reads.
sentieon_haplotyper_cmd = sentieon_task(
    "sentieon_haplotyper",
    SENTIEON_CALLING_IMAGE,
    MEDIUM_RESOURCES,
    inputs={
        "alignments": list[File],
        "indexes": list[File],
        "intervals": File | None,
        "recal_table": File | None,
        "fasta": File,
        "fai": File,
        "dbsnp": File | None,
        "dbsnp_tbi": File | None,
        "emit_vcf": str,
        "emit_gvcf": bool,
        "prefix": str,
        "args": str,
        "args2": str,
        "args3": str,
    },
    defaults={"args": "", "args2": "", "args3": ""},
    script=r"""
        ALN=({inputs.alignments}); IDX=({inputs.indexes}); FA=({inputs.fasta}); FAI=({inputs.fai})
        shopt -s nullglob
        IV=({inputs.intervals}); RT=({inputs.recal_table}); DBSNP=({inputs.dbsnp}); DBSNP_TBI=({inputs.dbsnp_tbi})
        shopt -u nullglob
        link "${ALN[@]}" "${IDX[@]}" "${FA[0]}" "${FAI[0]}" "${IV[@]}" "${RT[@]}" "${DBSNP[@]}" "${DBSNP_TBI[@]}"
        P={inputs.prefix}
        ARGS={inputs.args}; ARGS2={inputs.args2}; ARGS3={inputs.args3}
        EMIT_VCF={inputs.emit_vcf}
        INPUTS=(); for a in "${ALN[@]}"; do INPUTS+=(-i "$(basename "$a")"); done
        BASE=(--algo Haplotyper)
        if [ ${#DBSNP[@]} -gt 0 ]; then BASE+=(-d "$(basename "${DBSNP[0]}")"); fi
        INTERVAL=(); if [ ${#IV[@]} -gt 0 ]; then INTERVAL=(--interval "$(basename "${IV[0]}")"); fi
        RECAL=(); if [ ${#RT[@]} -gt 0 ]; then RECAL=(-q "$(basename "${RT[0]}")"); fi
        VCF=(); if [ -n "$EMIT_VCF" ]; then VCF=("${BASE[@]}" $ARGS2 --emit_mode "$EMIT_VCF" "$P.unfiltered.vcf.gz"); fi
        GVCF=(); if {inputs.emit_gvcf}; then GVCF=("${BASE[@]}" $ARGS3 --emit_mode gvcf "$P.g.vcf.gz"); fi
        sentieon driver $ARGS -r "$(basename "${FA[0]}")" -t $TASK_CPUS "${INTERVAL[@]}" "${INPUTS[@]}" \
            "${RECAL[@]}" "${VCF[@]}" "${GVCF[@]}"
        shopt -s nullglob
        mv ./*.unfiltered.vcf.gz* ./*.g.vcf.gz* {outputs.results}/
    """,
)

# As upstream, GVCFtyper runs without a thread count.
sentieon_gvcftyper_cmd = sentieon_task(
    "sentieon_gvcftyper",
    SENTIEON_CALLING_IMAGE,
    HIGH_RESOURCES,
    inputs={
        "gvcfs": list[File],
        "tbis": list[File],
        "intervals": File | None,
        "fasta": File,
        "fai": File,
        "dbsnp": File | None,
        "dbsnp_tbi": File | None,
        "prefix": str,
        "args": str,
    },
    defaults={"args": ""},
    script=r"""
        GVCFS=({inputs.gvcfs}); TBIS=({inputs.tbis}); FA=({inputs.fasta}); FAI=({inputs.fai})
        shopt -s nullglob
        IV=({inputs.intervals}); DBSNP=({inputs.dbsnp}); DBSNP_TBI=({inputs.dbsnp_tbi})
        shopt -u nullglob
        link "${GVCFS[@]}" "${TBIS[@]}" "${FA[0]}" "${FAI[0]}" "${IV[@]}" "${DBSNP[@]}" "${DBSNP_TBI[@]}"
        P={inputs.prefix}
        ARGS={inputs.args}
        INPUTS=(); for g in "${GVCFS[@]}"; do INPUTS+=(-v "$(basename "$g")"); done
        DB=(); if [ ${#DBSNP[@]} -gt 0 ]; then DB=(--dbsnp "$(basename "${DBSNP[0]}")"); fi
        INTERVAL=(); if [ ${#IV[@]} -gt 0 ]; then INTERVAL=(--interval "$(basename "${IV[0]}")"); fi
        sentieon driver -r "$(basename "${FA[0]}")" "${INTERVAL[@]}" --algo GVCFtyper $ARGS "${INPUTS[@]}" \
            "${DB[@]}" "$P.vcf.gz"
        mv "$P".vcf.gz "$P".vcf.gz.tbi {outputs.results}/
    """,
)

sentieon_dnascope_cmd = sentieon_task(
    "sentieon_dnascope",
    SENTIEON_CALLING_IMAGE,
    HIGH_RESOURCES,
    inputs={
        "bam": File,
        "bai": File,
        "intervals": File | None,
        "fasta": File,
        "fai": File,
        "dbsnp": File | None,
        "dbsnp_tbi": File | None,
        "ml_model": File | None,
        "pcr_indel_model": str,
        "emit_vcf": str,
        "emit_gvcf": bool,
        "prefix": str,
        "args": str,
        "args2": str,
        "args3": str,
    },
    defaults={"args": "", "args2": "", "args3": ""},
    script=r"""
        BAM=({inputs.bam}); BAI=({inputs.bai}); FA=({inputs.fasta}); FAI=({inputs.fai})
        shopt -s nullglob
        IV=({inputs.intervals}); DBSNP=({inputs.dbsnp}); DBSNP_TBI=({inputs.dbsnp_tbi}); MODEL=({inputs.ml_model})
        shopt -u nullglob
        link "${BAM[0]}" "${BAI[0]}" "${FA[0]}" "${FAI[0]}" "${IV[@]}" "${DBSNP[@]}" "${DBSNP_TBI[@]}" "${MODEL[@]}"
        P={inputs.prefix}
        ARGS={inputs.args}; ARGS2={inputs.args2}; ARGS3={inputs.args3}
        EMIT_VCF={inputs.emit_vcf}; PCR={inputs.pcr_indel_model}
        BASE=(--algo DNAscope)
        if [ ${#DBSNP[@]} -gt 0 ]; then BASE+=(-d "$(basename "${DBSNP[0]}")"); fi
        MODEL_OPTS=()
        if [ ${#MODEL[@]} -gt 0 ]; then MODEL_OPTS+=(--model "$(basename "${MODEL[0]}")"); fi
        if [ -n "$PCR" ]; then MODEL_OPTS+=(--pcr_indel_model "$PCR"); fi
        INTERVAL=(); if [ ${#IV[@]} -gt 0 ]; then INTERVAL=(--interval "$(basename "${IV[0]}")"); fi
        VCF=()
        if [ -n "$EMIT_VCF" ]; then
            VCF=("${BASE[@]}" $ARGS2 "${MODEL_OPTS[@]}" --emit_mode "$EMIT_VCF" "$P.unfiltered.vcf.gz")
        fi
        GVCF=()
        if {inputs.emit_gvcf}; then GVCF=("${BASE[@]}" $ARGS3 "${MODEL_OPTS[@]}" --emit_mode gvcf "$P.g.vcf.gz"); fi
        sentieon driver $ARGS -r "$(basename "${FA[0]}")" -t $TASK_CPUS -i "$(basename "${BAM[0]}")" \
            "${INTERVAL[@]}" "${VCF[@]}" "${GVCF[@]}"
        shopt -s nullglob
        mv ./*.unfiltered.vcf.gz* ./*.g.vcf.gz* {outputs.results}/
    """,
)

sentieon_dnamodelapply_cmd = sentieon_task(
    "sentieon_dnamodelapply",
    SENTIEON_CALLING_IMAGE,
    HIGH_RESOURCES,
    inputs={"vcf": File, "idx": File, "fasta": File, "fai": File, "ml_model": File, "prefix": str, "args": str},
    defaults={"args": ""},
    script=r"""
        VCF=({inputs.vcf}); IDX=({inputs.idx}); FA=({inputs.fasta}); FAI=({inputs.fai}); MODEL=({inputs.ml_model})
        link "${VCF[0]}" "${IDX[0]}" "${FA[0]}" "${FAI[0]}" "${MODEL[0]}"
        P={inputs.prefix}
        ARGS={inputs.args}
        sentieon driver -t $TASK_CPUS -r "$(basename "${FA[0]}")" $ARGS --algo DNAModelApply \
            --model "$(basename "${MODEL[0]}")" -v "$(basename "${VCF[0]}")" "$P.vcf.gz"
        mv "$P".vcf.gz "$P".vcf.gz.tbi {outputs.results}/
    """,
)

sentieon_tnscope_cmd = sentieon_task(
    "sentieon_tnscope",
    SENTIEON_CALLING_IMAGE,
    HIGH_RESOURCES,
    inputs={
        "alignments": list[File],
        "indexes": list[File],
        "intervals": File | None,
        "fasta": File,
        "fai": File,
        "dbsnp": File | None,
        "dbsnp_tbi": File | None,
        "pon": File | None,
        "pon_tbi": File | None,
        "cosmic": File | None,
        "cosmic_tbi": File | None,
        "prefix": str,
        "args": str,
        "args2": str,
    },
    defaults={"args": "", "args2": ""},
    script=r"""
        ALN=({inputs.alignments}); IDX=({inputs.indexes}); FA=({inputs.fasta}); FAI=({inputs.fai})
        shopt -s nullglob
        IV=({inputs.intervals}); DBSNP=({inputs.dbsnp}); DBSNP_TBI=({inputs.dbsnp_tbi})
        PON=({inputs.pon}); PON_TBI=({inputs.pon_tbi}); COSMIC=({inputs.cosmic}); COSMIC_TBI=({inputs.cosmic_tbi})
        shopt -u nullglob
        link "${ALN[@]}" "${IDX[@]}" "${FA[0]}" "${FAI[0]}" "${IV[@]}" "${DBSNP[@]}" "${DBSNP_TBI[@]}" \
            "${PON[@]}" "${PON_TBI[@]}" "${COSMIC[@]}" "${COSMIC_TBI[@]}"
        P={inputs.prefix}
        ARGS={inputs.args}; ARGS2={inputs.args2}
        INPUTS=(); for a in "${ALN[@]}"; do INPUTS+=(-i "$(basename "$a")"); done
        INTERVAL=(); if [ ${#IV[@]} -gt 0 ]; then INTERVAL=(--interval "$(basename "${IV[0]}")"); fi
        EXTRA=()
        if [ ${#COSMIC[@]} -gt 0 ]; then EXTRA+=(--cosmic "$(basename "${COSMIC[0]}")"); fi
        if [ ${#DBSNP[@]} -gt 0 ]; then EXTRA+=(--dbsnp "$(basename "${DBSNP[0]}")"); fi
        if [ ${#PON[@]} -gt 0 ]; then EXTRA+=(--pon "$(basename "${PON[0]}")"); fi
        sentieon driver -t $TASK_CPUS -r "$(basename "${FA[0]}")" "${INPUTS[@]}" "${INTERVAL[@]}" $ARGS \
            --algo TNscope $ARGS2 "${EXTRA[@]}" "$P.vcf.gz"
        mv "$P".vcf.gz "$P".vcf.gz.tbi {outputs.results}/
    """,
)

# `resources` is the resource options :func:`varcal_resource_options` builds
# from upstream's labels; they name the resource VCFs, which are linked into
# the working directory under their own names. As upstream, VarCal runs
# without a thread count.
sentieon_varcal_cmd = sentieon_task(
    "sentieon_varcal",
    SENTIEON_IMAGE,
    LOW_RESOURCES,
    inputs={
        "vcf": File,
        "tbi": File,
        "resource_vcfs": list[File],
        "resource_tbis": list[File],
        "resources": str,
        "fasta": File,
        "fai": File,
        "prefix": str,
        "args": str,
    },
    defaults={"args": ""},
    script=r"""
        VCF=({inputs.vcf}); TBI=({inputs.tbi}); FA=({inputs.fasta}); FAI=({inputs.fai})
        RES=({inputs.resource_vcfs}); RES_TBI=({inputs.resource_tbis})
        link "${VCF[0]}" "${TBI[0]}" "${FA[0]}" "${FAI[0]}" "${RES[@]}" "${RES_TBI[@]}"
        P={inputs.prefix}
        RESOURCES={inputs.resources}
        ARGS={inputs.args}
        sentieon driver -r "$(basename "${FA[0]}")" --algo VarCal -v "$(basename "${VCF[0]}")" \
            --tranches_file "$P.tranches" $RESOURCES $ARGS "$P.recal"
        shopt -s nullglob
        mv "$P".recal "$P".recal.idx "$P".tranches ./*plots.R {outputs.results}/
    """,
)

sentieon_applyvarcal_cmd = sentieon_task(
    "sentieon_applyvarcal",
    SENTIEON_CALLING_IMAGE,
    LOW_RESOURCES,
    inputs={
        "vcf": File,
        "vcf_tbi": File,
        "recal": File,
        "recal_index": File,
        "tranches": File,
        "fasta": File,
        "fai": File,
        "prefix": str,
        "args": str,
        "args2": str,
    },
    defaults={"args": "", "args2": ""},
    script=r"""
        VCF=({inputs.vcf}); TBI=({inputs.vcf_tbi}); RECAL=({inputs.recal}); RECAL_IDX=({inputs.recal_index})
        TRANCHES=({inputs.tranches}); FA=({inputs.fasta}); FAI=({inputs.fai})
        link "${VCF[0]}" "${TBI[0]}" "${RECAL[0]}" "${RECAL_IDX[0]}" "${TRANCHES[0]}" "${FA[0]}" "${FAI[0]}"
        P={inputs.prefix}
        ARGS={inputs.args}; ARGS2={inputs.args2}
        sentieon driver -r "$(basename "${FA[0]}")" -t $TASK_CPUS $ARGS --algo ApplyVarCal \
            -v "$(basename "${VCF[0]}")" --recal "$(basename "${RECAL[0]}")" \
            --tranches_file "$(basename "${TRANCHES[0]}")" $ARGS2 "$P.vcf.gz"
        mv "$P".vcf.gz "$P".vcf.gz.tbi {outputs.results}/
    """,
)

# Mates are separate lists (reads_2 empty for single-end); each list is staged
# under its files' names and globbed alphabetically, then joined with commas
# as STAR takes several runs. As upstream, the options default to an
# unsorted BAM and a read group of ID:<prefix> SM:<prefix> unless args set them.
sentieon_staralign_cmd = sentieon_task(
    "sentieon_staralign",
    SENTIEON_CALLING_IMAGE,
    HIGH_RESOURCES,
    inputs={
        "reads_1": list[File],
        "reads_2": list[File],
        "index": Dir,
        "gtf": File | None,
        "ignore_sjdbgtf": bool,
        "prefix": str,
        "args": str,
    },
    defaults={"reads_2": [], "ignore_sjdbgtf": False, "args": ""},
    script=r"""
        shopt -s nullglob; R1=({inputs.reads_1}); R2=({inputs.reads_2}); GTF=({inputs.gtf}); shopt -u nullglob
        P={inputs.prefix}
        ARGS={inputs.args}
        READS=("$(IFS=,; echo "${R1[*]}")")
        if [ ${#R2[@]} -gt 0 ]; then READS+=("$(IFS=,; echo "${R2[*]}")"); fi
        OPTS=()
        if [[ "$ARGS" != *--outSAMtype* ]]; then OPTS+=(--outSAMtype BAM Unsorted); fi
        if ! {inputs.ignore_sjdbgtf}; then OPTS+=(--sjdbGTFfile "${GTF[0]}"); fi
        if [[ "$ARGS" != *--outSAMattrRGline* ]]; then OPTS+=(--outSAMattrRGline "ID:$P" "SM:$P"); fi
        sentieon STAR --genomeDir {inputs.index} --readFilesIn "${READS[@]}" --runThreadN $TASK_CPUS \
            --outFileNamePrefix "$P." "${OPTS[@]}" $ARGS
        if [[ "$ARGS" == *"--outSAMtype BAM Unsorted SortedByCoordinate"* ]]; then
            mv "$P.Aligned.out.bam" "$P.Aligned.unsort.out.bam"
        fi
        if [ -f "$P.Unmapped.out.mate1" ]; then
            mv "$P.Unmapped.out.mate1" "$P.unmapped_1.fastq"; gzip "$P.unmapped_1.fastq"
        fi
        if [ -f "$P.Unmapped.out.mate2" ]; then
            mv "$P.Unmapped.out.mate2" "$P.unmapped_2.fastq"; gzip "$P.unmapped_2.fastq"
        fi
        find . -maxdepth 1 -type f -name "$P.*" -exec mv {} {outputs.results}/ \;
    """,
)

# The reference Dir is upstream's `rsem/`: the FASTA, RSEM's genome.* and, with
# `--star` in args, STAR's index, built by Sentieon's STAR (linked as `STAR`).
sentieon_rsempreparereference_cmd = sentieon_task(
    "sentieon_rsempreparereference",
    RSEM_SENTIEON_PREPARE_IMAGE,
    HIGH_RESOURCES,
    inputs={"fasta": File, "gtf": File, "args": str, "args2": str},
    defaults={"args": "", "args2": ""},
    script=r"""
        FA=({inputs.fasta}); GTF=({inputs.gtf})
        OUT={outputs.results}
        NAME=$(basename "${FA[0]}")
        cp "${FA[0]}" "$OUT/$NAME"
        ARGS_STR={inputs.args}
        ARGS=($ARGS_STR)
        ln -sf "$(which sentieon)" ./STAR
        export PATH=".:$PATH"
        # As upstream: `--star` builds STAR's index too, and every option naming --star is left out of RSEM's.
        STAR_MODE=false; RSEM_ARGS=()
        for a in "${ARGS[@]}"; do
            if [ "$a" = --star ]; then STAR_MODE=true; fi
            if [[ "$a" != *--star* ]]; then RSEM_ARGS+=("$a"); fi
        done
        if [ $STAR_MODE = true ]; then
            ARGS2={inputs.args2}
            STAR --runMode genomeGenerate --genomeDir "$OUT/" --genomeFastaFiles "$OUT/$NAME" \
                --sjdbGTFfile "${GTF[0]}" --runThreadN $TASK_CPUS \
                --limitGenomeGenerateRAM $(( TASK_MEMORY_MB * 1048576 - 100000000 )) $ARGS2
            rsem-prepare-reference --gtf "${GTF[0]}" --num-threads $TASK_CPUS "${RSEM_ARGS[@]}" \
                "$OUT/$NAME" "$OUT/genome"
        else
            rsem-prepare-reference --gtf "${GTF[0]}" --num-threads $TASK_CPUS "${ARGS[@]}" "$OUT/$NAME" "$OUT/genome"
        fi
    """,
)

# `single_end` is upstream's meta.single_end: "true", "false", or "" when
# unknown, which detects pairing from the reads (as upstream, an undetected
# pairing passes the literal "unknown" on to RSEM).
sentieon_rsemcalculateexpression_cmd = sentieon_task(
    "sentieon_rsemcalculateexpression",
    RSEM_SENTIEON_CALCULATE_IMAGE,
    HIGH_RESOURCES,
    inputs={
        "reads_1": File,
        "reads_2": File | None,
        "index": Dir,
        "strandedness": str,
        "single_end": str,
        "prefix": str,
        "args": str,
    },
    defaults={"strandedness": "", "single_end": "", "args": ""},
    script=r"""
        R1=({inputs.reads_1}); shopt -s nullglob; R2=({inputs.reads_2}); shopt -u nullglob
        READS=("${R1[0]}" "${R2[@]}")
        P={inputs.prefix}
        ARGS={inputs.args}
        if [[ "$ARGS" != *--temporary-folder* ]]; then ARGS="$ARGS --temporary-folder ./tmp/"; fi
        STRANDED={inputs.strandedness}
        STRANDEDNESS=()
        if [ "$STRANDED" = forward ] || [ "$STRANDED" = reverse ]; then STRANDEDNESS=(--strandedness "$STRANDED"); fi
        IS_BAM=false
        if [ ${#READS[@]} -eq 1 ] && [[ "${READS[0],,}" == *.bam ]]; then IS_BAM=true; fi
        ALIGNMENT=(); if [ $IS_BAM = true ]; then ALIGNMENT=(--alignments); fi
        case {inputs.single_end} in
            true) PAIRED_END="" ;;
            false) PAIRED_END="--paired-end" ;;
            *) PAIRED_END="unknown" ;;
        esac
        INDEX=$(find -L {inputs.index} -name "*.grp" | sed 's/\.grp$//')
        ln -sf "$(which sentieon)" ./STAR
        export PATH=".:$PATH"
        PAIRED_END_FLAG="$PAIRED_END"
        if [ "$PAIRED_END" == "unknown" ]; then
            if [ $IS_BAM = true ]; then
                samtools flagstat "${READS[0]}" | grep -q 'paired in sequencing' && PAIRED_END_FLAG="--paired-end"
            else
                [ ${#READS[@]} -gt 1 ] && PAIRED_END_FLAG="--paired-end"
            fi
        fi
        rsem-calculate-expression --num-threads $TASK_CPUS "${ALIGNMENT[@]}" "${STRANDEDNESS[@]}" \
            $PAIRED_END_FLAG $ARGS "${READS[@]}" "$INDEX" "$P"
        shopt -s nullglob
        mv "$P".genes.results "$P".isoforms.results "$P".stat ./*.log "$P".*.bam {outputs.results}/
    """,
)

sentieon_env = flyte.TaskEnvironment.from_task(
    "sentieon",
    sentieon_bwaindex_cmd.as_task(),
    sentieon_readwriter_cmd.as_task(),
    sentieon_dedup_cmd.as_task(),
    sentieon_wgsmetrics_cmd.as_task(),
    sentieon_varcal_cmd.as_task(),
)
sentieon_calling_env = flyte.TaskEnvironment.from_task(
    "sentieon_calling",
    sentieon_haplotyper_cmd.as_task(),
    sentieon_gvcftyper_cmd.as_task(),
    sentieon_dnascope_cmd.as_task(),
    sentieon_dnamodelapply_cmd.as_task(),
    sentieon_tnscope_cmd.as_task(),
    sentieon_applyvarcal_cmd.as_task(),
    sentieon_staralign_cmd.as_task(),
)
sentieon_findutils_env = flyte.TaskEnvironment.from_task("sentieon_findutils", sentieon_bwamem_cmd.as_task())
sentieon_gnuplot_env = flyte.TaskEnvironment.from_task(
    "sentieon_gnuplot", sentieon_datametrics_cmd.as_task(), sentieon_qualcal_cmd.as_task()
)
rsem_sentieon_prepare_env = flyte.TaskEnvironment.from_task(
    "rsem_sentieon_prepare", sentieon_rsempreparereference_cmd.as_task()
)
rsem_sentieon_calculate_env = flyte.TaskEnvironment.from_task(
    "rsem_sentieon_calculate", sentieon_rsemcalculateexpression_cmd.as_task()
)
env = flyte.TaskEnvironment(
    name="sentieon_modules",
    depends_on=[
        sentieon_env,
        sentieon_calling_env,
        sentieon_findutils_env,
        sentieon_gnuplot_env,
        rsem_sentieon_prepare_env,
        rsem_sentieon_calculate_env,
    ],
)


def licensed(cmd, secrets: list[flyte.Secret] | None):
    """``cmd``, requesting ``secrets`` instead of :data:`DEFAULT_SECRETS` when given."""
    return cmd if secrets is None else cmd.override(secrets=secrets)


async def pick(results: Dir, name: str, tool: str) -> File:
    found = await results.get_file(name)
    if found is None:
        raise FileNotFoundError(f"Sentieon {tool} wrote no {name}")
    return found


async def matching(results: Dir, pattern: str) -> list[File]:
    """The files in ``results`` whose names match the glob ``pattern``, by name."""
    found = [f async for f in results.walk() if fnmatch(f.name, pattern)]
    return sorted(found, key=lambda f: f.name)


def index_name(prefix: str) -> str:
    """An alignment's index name, as upstream: ``.bai`` for a BAM, else ``.crai``."""
    return f"{prefix}.bai" if prefix.rsplit(".", 1)[-1] == "bam" else f"{prefix}.crai"


@dataclass
class AlignmentResult:
    alignment: File  # <prefix>: .bam or .cram
    index: File  # <prefix>.bai or <prefix>.crai


async def sentieon_bwaindex(
    fasta: File, prefix: str = "", args: str = "", secrets: list[flyte.Secret] | None = None
) -> Dir:
    """A BWA index of ``fasta``: a Dir of ``<prefix>.{amb,ann,bwt,pac,sa}`` (prefix: the FASTA's base name)."""
    return await licensed(sentieon_bwaindex_cmd, secrets)(fasta=fasta, prefix=prefix, args=args)


async def sentieon_bwamem(
    reads_1: File,
    index: Dir,
    fasta: File,
    prefix: str,
    reads_2: File | None = None,
    fai: File | None = None,
    args: str = "",
    secrets: list[flyte.Secret] | None = None,
) -> AlignmentResult:
    """Align reads (``reads_2`` None for single-end) and sort them into ``prefix`` (``<id>.bam`` or ``<id>.cram``).

    The read group goes in ``args`` (``-R @RG\\tID:...``), without quotes.
    """
    results = await licensed(sentieon_bwamem_cmd, secrets)(
        reads_1=reads_1, reads_2=reads_2, index=index, fasta=fasta, fai=fai, prefix=prefix, args=args
    )
    return AlignmentResult(
        alignment=await pick(results, prefix, "bwa mem"), index=await pick(results, index_name(prefix), "bwa mem")
    )


async def sentieon_readwriter(
    alignments: list[File],
    indexes: list[File],
    prefix: str,
    fasta: File | None = None,
    fai: File | None = None,
    args: str = "",
    args2: str = "",
    secrets: list[flyte.Secret] | None = None,
) -> AlignmentResult:
    """Merge (and convert) alignments into ``prefix`` (``<id>.bam`` or ``<id>.cram``) with its index."""
    results = await licensed(sentieon_readwriter_cmd, secrets)(
        alignments=alignments, indexes=indexes, fasta=fasta, fai=fai, prefix=prefix, args=args, args2=args2
    )
    return AlignmentResult(
        alignment=await pick(results, prefix, "ReadWriter"),
        index=await pick(results, index_name(prefix), "ReadWriter"),
    )


@dataclass
class DedupResult:
    alignment: File  # <prefix>: .cram or .bam
    index: File  # <prefix>.crai or <prefix>.bai
    score: File  # <prefix without extension>.score
    metrics: File  # <prefix>.metrics
    metrics_multiqc_tsv: File  # <prefix>.metrics.multiqc.tsv


async def sentieon_dedup(
    bams: list[File],
    bais: list[File],
    fasta: File,
    fai: File,
    prefix: str,
    args: str = "",
    args2: str = "",
    args3: str = "",
    args4: str = "",
    secrets: list[flyte.Secret] | None = None,
) -> DedupResult:
    """Mark duplicates (``--rmdup`` in ``args4`` removes them) into ``prefix`` (``<id>.cram`` upstream)."""
    results = await licensed(sentieon_dedup_cmd, secrets)(
        bams=bams, bais=bais, fasta=fasta, fai=fai, prefix=prefix, args=args, args2=args2, args3=args3, args4=args4
    )
    base = prefix.rsplit(".", 1)[0]
    return DedupResult(
        alignment=await pick(results, prefix, "Dedup"),
        index=await pick(results, index_name(prefix), "Dedup"),
        score=await pick(results, f"{base}.score", "Dedup"),
        metrics=await pick(results, f"{prefix}.metrics", "Dedup"),
        metrics_multiqc_tsv=await pick(results, f"{prefix}.metrics.multiqc.tsv", "Dedup"),
    )


@dataclass
class DataMetricsResult:
    mq_metrics: File
    qd_metrics: File
    gc_summary: File
    gc_metrics: File
    aln_metrics: File
    is_metrics: File
    mq_plot: File | None  # the PDFs, with plot_results
    qd_plot: File | None
    is_plot: File | None
    gc_plot: File | None


async def sentieon_datametrics(
    bams: list[File],
    bais: list[File],
    fasta: File,
    fai: File,
    prefix: str,
    plot_results: bool = False,
    args: str = "",
    secrets: list[flyte.Secret] | None = None,
) -> DataMetricsResult:
    """GC bias, quality by cycle, quality distribution, insert size and alignment metrics, ``<prefix>_*``."""
    results = await licensed(sentieon_datametrics_cmd, secrets)(
        bams=bams, bais=bais, fasta=fasta, fai=fai, plot_results=plot_results, prefix=prefix, args=args
    )

    async def metrics(name: str) -> File:
        return await pick(results, f"{prefix}_{name}.txt", "metrics")

    async def plot(name: str) -> File | None:
        return await pick(results, f"{prefix}_{name}.pdf", "plot") if plot_results else None

    return DataMetricsResult(
        mq_metrics=await metrics("mq_metrics"),
        qd_metrics=await metrics("qd_metrics"),
        gc_summary=await metrics("gc_summary"),
        gc_metrics=await metrics("gc_metrics"),
        aln_metrics=await metrics("aln_metrics"),
        is_metrics=await metrics("is_metrics"),
        mq_plot=await plot("mq_metrics"),
        qd_plot=await plot("qd_metrics"),
        is_plot=await plot("is_metrics"),
        gc_plot=await plot("gc_metrics"),
    )


async def sentieon_wgsmetrics(
    bams: list[File],
    bais: list[File],
    fasta: File,
    fai: File,
    prefix: str,
    intervals: File | None = None,
    args: str = "",
    secrets: list[flyte.Secret] | None = None,
) -> File:
    """Whole-genome coverage metrics, ``<prefix>.txt``."""
    results = await licensed(sentieon_wgsmetrics_cmd, secrets)(
        bams=bams, bais=bais, fasta=fasta, fai=fai, intervals=intervals, prefix=prefix, args=args
    )
    return await pick(results, f"{prefix}.txt", "WgsMetricsAlgo")


@dataclass
class QualCalResult:
    table: File | None  # <prefix>.table, without a recalibration table
    table_post: File | None  # <prefix>.table.post, with one
    recal_alignment: File | None  # <prefix>.recalibrated.{bam,cram}, with generate_recalibrated_bams
    csv: File | None  # <prefix>.csv, with a recalibration table
    pdf: File | None  # <prefix>.pdf, with a recalibration table


async def sentieon_qualcal(
    alignments: list[File],
    indexes: list[File],
    fasta: File,
    fai: File,
    prefix: str,
    known_sites: list[File] | None = None,
    known_sites_tbi: list[File] | None = None,
    recalibration_table: File | None = None,
    generate_recalibrated_bams: bool = False,
    args: str = "",
    secrets: list[flyte.Secret] | None = None,
) -> QualCalResult:
    """A base-quality recalibration table; given ``recalibration_table``, its before/after report."""
    results = await licensed(sentieon_qualcal_cmd, secrets)(
        alignments=alignments,
        indexes=indexes,
        fasta=fasta,
        fai=fai,
        known_sites=known_sites or [],
        known_sites_tbi=known_sites_tbi or [],
        recalibration_table=recalibration_table,
        generate_recalibrated_bams=generate_recalibrated_bams,
        prefix=prefix,
        args=args,
    )
    recal = await matching(results, f"{prefix}.recalibrated.*")
    return QualCalResult(
        table=await results.get_file(f"{prefix}.table"),
        table_post=await results.get_file(f"{prefix}.table.post"),
        recal_alignment=recal[0] if recal else None,
        csv=await results.get_file(f"{prefix}.csv"),
        pdf=await results.get_file(f"{prefix}.pdf"),
    )


@dataclass
class GermlineCallsResult:
    vcf: File | None  # <prefix>.unfiltered.vcf.gz, with emit_vcf
    vcf_tbi: File | None
    gvcf: File | None  # <prefix>.g.vcf.gz, with emit_gvcf
    gvcf_tbi: File | None


async def germline_calls(results: Dir, prefix: str) -> GermlineCallsResult:
    return GermlineCallsResult(
        vcf=await results.get_file(f"{prefix}.unfiltered.vcf.gz"),
        vcf_tbi=await results.get_file(f"{prefix}.unfiltered.vcf.gz.tbi"),
        gvcf=await results.get_file(f"{prefix}.g.vcf.gz"),
        gvcf_tbi=await results.get_file(f"{prefix}.g.vcf.gz.tbi"),
    )


async def sentieon_haplotyper(
    alignments: list[File],
    indexes: list[File],
    fasta: File,
    fai: File,
    prefix: str,
    emit_vcf: str = "",
    emit_gvcf: bool = False,
    intervals: File | None = None,
    recal_table: File | None = None,
    dbsnp: File | None = None,
    dbsnp_tbi: File | None = None,
    args: str = "",
    args2: str = "",
    args3: str = "",
    secrets: list[flyte.Secret] | None = None,
) -> GermlineCallsResult:
    """Germline calls: a VCF with ``emit_vcf`` (its ``--emit_mode``, e.g. ``variant``), a GVCF with ``emit_gvcf``.

    ``args`` go to the driver, ``args2`` to the VCF pass, ``args3`` to the GVCF pass.
    """
    results = await licensed(sentieon_haplotyper_cmd, secrets)(
        alignments=alignments,
        indexes=indexes,
        intervals=intervals,
        recal_table=recal_table,
        fasta=fasta,
        fai=fai,
        dbsnp=dbsnp,
        dbsnp_tbi=dbsnp_tbi,
        emit_vcf=emit_vcf,
        emit_gvcf=emit_gvcf,
        prefix=prefix,
        args=args,
        args2=args2,
        args3=args3,
    )
    return await germline_calls(results, prefix)


@dataclass
class VcfResult:
    vcf: File  # <prefix>.vcf.gz
    tbi: File  # <prefix>.vcf.gz.tbi


async def vcf_result(results: Dir, prefix: str, tool: str) -> VcfResult:
    return VcfResult(
        vcf=await pick(results, f"{prefix}.vcf.gz", tool), tbi=await pick(results, f"{prefix}.vcf.gz.tbi", tool)
    )


async def sentieon_gvcftyper(
    gvcfs: list[File],
    tbis: list[File],
    fasta: File,
    fai: File,
    prefix: str,
    intervals: File | None = None,
    dbsnp: File | None = None,
    dbsnp_tbi: File | None = None,
    args: str = "",
    secrets: list[flyte.Secret] | None = None,
) -> VcfResult:
    """Joint-genotype GVCFs into ``<prefix>.vcf.gz`` (upstream's prefix: ``<id>_genotyped``)."""
    results = await licensed(sentieon_gvcftyper_cmd, secrets)(
        gvcfs=gvcfs,
        tbis=tbis,
        intervals=intervals,
        fasta=fasta,
        fai=fai,
        dbsnp=dbsnp,
        dbsnp_tbi=dbsnp_tbi,
        prefix=prefix,
        args=args,
    )
    return await vcf_result(results, prefix, "GVCFtyper")


async def sentieon_dnascope(
    bam: File,
    bai: File,
    fasta: File,
    fai: File,
    prefix: str,
    emit_vcf: str = "",
    emit_gvcf: bool = False,
    pcr_indel_model: str = "",
    ml_model: File | None = None,
    intervals: File | None = None,
    dbsnp: File | None = None,
    dbsnp_tbi: File | None = None,
    args: str = "",
    args2: str = "",
    args3: str = "",
    secrets: list[flyte.Secret] | None = None,
) -> GermlineCallsResult:
    """DNAscope germline calls, as :func:`sentieon_haplotyper`, optionally with a model and PCR indel model."""
    results = await licensed(sentieon_dnascope_cmd, secrets)(
        bam=bam,
        bai=bai,
        intervals=intervals,
        fasta=fasta,
        fai=fai,
        dbsnp=dbsnp,
        dbsnp_tbi=dbsnp_tbi,
        ml_model=ml_model,
        pcr_indel_model=pcr_indel_model,
        emit_vcf=emit_vcf,
        emit_gvcf=emit_gvcf,
        prefix=prefix,
        args=args,
        args2=args2,
        args3=args3,
    )
    return await germline_calls(results, prefix)


async def sentieon_dnamodelapply(
    vcf: File,
    idx: File,
    fasta: File,
    fai: File,
    ml_model: File,
    prefix: str,
    args: str = "",
    secrets: list[flyte.Secret] | None = None,
) -> VcfResult:
    """Apply a DNAscope model to ``vcf`` (upstream's prefix: ``<id>_applied``)."""
    results = await licensed(sentieon_dnamodelapply_cmd, secrets)(
        vcf=vcf, idx=idx, fasta=fasta, fai=fai, ml_model=ml_model, prefix=prefix, args=args
    )
    return await vcf_result(results, prefix, "DNAModelApply")


async def sentieon_tnscope(
    alignments: list[File],
    indexes: list[File],
    fasta: File,
    fai: File,
    prefix: str,
    intervals: File | None = None,
    dbsnp: File | None = None,
    dbsnp_tbi: File | None = None,
    pon: File | None = None,
    pon_tbi: File | None = None,
    cosmic: File | None = None,
    cosmic_tbi: File | None = None,
    args: str = "",
    args2: str = "",
    secrets: list[flyte.Secret] | None = None,
) -> VcfResult:
    """Somatic calls (``--tumor_sample``/``--normal_sample`` go in ``args2``) into ``<prefix>.vcf.gz``."""
    results = await licensed(sentieon_tnscope_cmd, secrets)(
        alignments=alignments,
        indexes=indexes,
        intervals=intervals,
        fasta=fasta,
        fai=fai,
        dbsnp=dbsnp,
        dbsnp_tbi=dbsnp_tbi,
        pon=pon,
        pon_tbi=pon_tbi,
        cosmic=cosmic,
        cosmic_tbi=cosmic_tbi,
        prefix=prefix,
        args=args,
        args2=args2,
    )
    return await vcf_result(results, prefix, "TNscope")


def varcal_resource_options(labels: list[str]) -> str:
    """VarCal's resource options from upstream's labels, as upstream builds them.

    Each label holds one or more ``--resource:<name>,<params> <file>`` directives;
    each becomes ``--resource <file> --resource_param <name>,<params>``.
    """
    options = []
    for label in labels:
        for directive in (part for part in label.split("--resource:") if part):
            items = directive.split(" ", 1)
            if len(items) != 2:
                raise ValueError(
                    f"Expected the resource string {directive!r} to contain two elements separated by a space."
                )
            options.append(f"--resource {items[1]} --resource_param {items[0]}")
    return " ".join(options)


@dataclass
class VarCalResult:
    recal: File  # <prefix>.recal
    idx: File  # <prefix>.recal.idx
    tranches: File  # <prefix>.tranches
    plots: File | None  # *plots.R, when args ask for it


async def sentieon_varcal(
    vcf: File,
    tbi: File,
    resource_vcfs: list[File],
    resource_tbis: list[File],
    labels: list[str],
    fasta: File,
    fai: File,
    prefix: str,
    args: str = "",
    secrets: list[flyte.Secret] | None = None,
) -> VarCalResult:
    """A VQSR model of ``vcf`` trained on the labelled resources (labels name the resource files)."""
    results = await licensed(sentieon_varcal_cmd, secrets)(
        vcf=vcf,
        tbi=tbi,
        resource_vcfs=resource_vcfs,
        resource_tbis=resource_tbis,
        resources=varcal_resource_options(labels),
        fasta=fasta,
        fai=fai,
        prefix=prefix,
        args=args,
    )
    plots = await matching(results, "*plots.R")
    return VarCalResult(
        recal=await pick(results, f"{prefix}.recal", "VarCal"),
        idx=await pick(results, f"{prefix}.recal.idx", "VarCal"),
        tranches=await pick(results, f"{prefix}.tranches", "VarCal"),
        plots=plots[0] if plots else None,
    )


async def sentieon_applyvarcal(
    vcf: File,
    vcf_tbi: File,
    recal: File,
    recal_index: File,
    tranches: File,
    fasta: File,
    fai: File,
    prefix: str,
    args: str = "",
    args2: str = "",
    secrets: list[flyte.Secret] | None = None,
) -> VcfResult:
    """Apply a VQSR model to ``vcf`` (upstream's prefix: ``<id>_applyvarcal``)."""
    results = await licensed(sentieon_applyvarcal_cmd, secrets)(
        vcf=vcf,
        vcf_tbi=vcf_tbi,
        recal=recal,
        recal_index=recal_index,
        tranches=tranches,
        fasta=fasta,
        fai=fai,
        prefix=prefix,
        args=args,
        args2=args2,
    )
    return await vcf_result(results, prefix, "ApplyVarCal")


@dataclass
class StarAlignResult:
    log_final: File  # <prefix>.Log.final.out
    log_out: File  # <prefix>.Log.out
    log_progress: File  # <prefix>.Log.progress.out
    bam: list[File]  # *d.out.bam
    bam_sorted_aligned: File | None  # <prefix>.Aligned.sortedByCoord.out.bam
    bam_transcript: File | None  # *toTranscriptome.out.bam
    bam_unsorted: File | None  # *Aligned.unsort.out.bam
    fastq: list[File]  # unmapped mates, <prefix>.unmapped_{1,2}.fastq.gz
    tab: list[File]  # *.tab
    spl_junc_tab: File | None  # <prefix>.SJ.out.tab
    read_per_gene_tab: File | None  # <prefix>.ReadsPerGene.out.tab
    junction: File | None  # *.out.junction
    sam: File | None  # *.out.sam
    wig: list[File]
    bedgraph: list[File]
    results: Dir  # everything STAR wrote


async def sentieon_staralign(
    reads_1: list[File],
    index: Dir,
    prefix: str,
    reads_2: list[File] | None = None,
    gtf: File | None = None,
    ignore_sjdbgtf: bool = False,
    args: str = "",
    secrets: list[flyte.Secret] | None = None,
) -> StarAlignResult:
    """Align RNA-seq reads with Sentieon's STAR; ``reads_2`` is None for single-end."""
    results = await licensed(sentieon_staralign_cmd, secrets)(
        reads_1=reads_1,
        reads_2=reads_2 or [],
        index=index,
        gtf=gtf,
        ignore_sjdbgtf=ignore_sjdbgtf,
        prefix=prefix,
        args=args,
    )

    async def one(pattern: str) -> File | None:
        found = await matching(results, pattern)
        return found[0] if found else None

    return StarAlignResult(
        log_final=await pick(results, f"{prefix}.Log.final.out", "STAR"),
        log_out=await pick(results, f"{prefix}.Log.out", "STAR"),
        log_progress=await pick(results, f"{prefix}.Log.progress.out", "STAR"),
        bam=await matching(results, "*d.out.bam"),
        bam_sorted_aligned=await results.get_file(f"{prefix}.Aligned.sortedByCoord.out.bam"),
        bam_transcript=await one("*toTranscriptome.out.bam"),
        bam_unsorted=await one("*Aligned.unsort.out.bam"),
        fastq=await matching(results, "*fastq.gz"),
        tab=await matching(results, "*.tab"),
        spl_junc_tab=await one("*.SJ.out.tab"),
        read_per_gene_tab=await one("*.ReadsPerGene.out.tab"),
        junction=await one("*.out.junction"),
        sam=await one("*.out.sam"),
        wig=await matching(results, "*.wig"),
        bedgraph=await matching(results, "*.bg"),
        results=results,
    )


@dataclass
class RsemReferenceResult:
    index: Dir  # upstream's rsem/: the FASTA, genome.* and, with --star, STAR's index
    transcript_fasta: File  # genome.transcripts.fa


async def sentieon_rsempreparereference(
    fasta: File, gtf: File, args: str = "", args2: str = "", secrets: list[flyte.Secret] | None = None
) -> RsemReferenceResult:
    """An RSEM reference; ``--star`` in ``args`` also builds STAR's index (``args2`` goes to STAR)."""
    index = await licensed(sentieon_rsempreparereference_cmd, secrets)(fasta=fasta, gtf=gtf, args=args, args2=args2)
    return RsemReferenceResult(index=index, transcript_fasta=await pick(index, "genome.transcripts.fa", "RSEM"))


@dataclass
class RsemExpressionResult:
    counts_gene: File  # <prefix>.genes.results
    counts_transcript: File  # <prefix>.isoforms.results
    stat: list[File]  # <prefix>.stat/: <prefix>.{cnt,model,theta}
    logs: list[File]  # *.log
    bam_star: File | None  # <prefix>.STAR.genome.bam
    bam_genome: File | None  # <prefix>.genome.bam
    bam_transcript: File | None  # <prefix>.transcript.bam


async def sentieon_rsemcalculateexpression(
    reads_1: File,
    index: Dir,
    prefix: str,
    reads_2: File | None = None,
    strandedness: str = "",
    single_end: bool | None = None,
    args: str = "",
    secrets: list[flyte.Secret] | None = None,
) -> RsemExpressionResult:
    """Quantify expression from FASTQ reads (``reads_2`` for paired-end) or one BAM (``--alignments``).

    ``strandedness`` is ``forward``, ``reverse`` or empty; ``single_end`` None detects pairing.
    """
    results = await licensed(sentieon_rsemcalculateexpression_cmd, secrets)(
        reads_1=reads_1,
        reads_2=reads_2,
        index=index,
        strandedness=strandedness,
        single_end="" if single_end is None else str(single_end).lower(),
        prefix=prefix,
        args=args,
    )
    stat = sorted(
        [f async for f in results.walk() if f.path.rstrip("/").rsplit("/", 2)[-2] == f"{prefix}.stat"],
        key=lambda f: f.name,
    )
    return RsemExpressionResult(
        counts_gene=await pick(results, f"{prefix}.genes.results", "RSEM"),
        counts_transcript=await pick(results, f"{prefix}.isoforms.results", "RSEM"),
        stat=stat,
        logs=await matching(results, "*.log"),
        bam_star=await results.get_file(f"{prefix}.STAR.genome.bam"),
        bam_genome=await results.get_file(f"{prefix}.genome.bam"),
        bam_transcript=await results.get_file(f"{prefix}.transcript.bam"),
    )
