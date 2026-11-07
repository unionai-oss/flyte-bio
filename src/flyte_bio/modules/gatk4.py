"""gatk4 — GATK tools.

Exposes:

- :func:`gatk4_createsequencedictionary` — the ``.dict`` sequence dictionary
  of a FASTA.
- :func:`gatk4_intervallisttobed` — convert a Picard interval list to BED.
- :func:`gatk4_bedtointervallist` — convert a BED to a Picard interval list
  against a sequence dictionary.
- :func:`gatk4_markduplicates` — mark duplicates across one or more BAMs
  (e.g. a sample's lanes), writing BAM or, as upstream does, CRAM via
  samtools.
- :func:`gatk4_baserecalibrator` — a base-quality recalibration table from
  known sites (optionally within intervals).
- :func:`gatk4_gatherbqsrreports` — merge per-interval recalibration tables.
- :func:`gatk4_applybqsr` — apply a recalibration table (optionally within
  intervals), writing CRAM or BAM.
- :func:`gatk4_haplotypecaller` — call germline SNVs and indels (optionally
  within intervals, annotating dbSNP IDs).
- :func:`gatk4_mergevcfs` — merge per-interval VCFs.
- :func:`gatk4_cnnscorevariants` — score variants with GATK's 1D CNN.
- :func:`gatk4_filtervarianttranches` — filter CNN-scored variants by
  sensitivity tranches against known-sites resources.
- :func:`gatk4_mutect2` — call somatic SNVs and indels (tumour-only, tumour
  with a matched normal, mitochondrial mode, or force-calling given alleles):
  ``<prefix>.vcf.gz`` with its ``.tbi`` and ``.stats``, and
  ``<prefix>.f1r2.tar.gz`` when ``args`` asks for it.
- :func:`gatk4_mergemutectstats` — merge per-interval Mutect2 stats into
  ``<prefix>.vcf.gz.stats``.
- :func:`gatk4_filtermutectcalls` — filter Mutect2 calls (with optional
  orientation-bias priors, tumour segmentation and contamination), writing
  ``<prefix>.vcf.gz`` (+ ``.tbi``) and ``<prefix>.vcf.gz.filteringStats.tsv``.
- :func:`gatk4_variantfiltration` — apply hard filters (named JEXL expressions) to
  a VCF, writing ``<prefix>.vcf.gz`` (+ ``.tbi``).
- :func:`gatk4_intervallisttools` — split an interval list for scattering,
  as upstream into ``<prefix>_split/temp_*/<n>scattered.interval_list``.

The Mutect2 and VariantFiltration references may be
bgzipped FASTA (pass the ``.gzi``); the dictionary is linked under the name
GATK looks for (``genome.dict`` for ``genome.fasta.gz``), as upstream.

As upstream, each tool gets 80% of the task's memory as Java heap. Inputs
whose names GATK records (the dictionary's ``UR`` field) are linked into the
working directory under their own names first, so outputs match upstream's.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import JAVA_HEAP_MB, TASK_LIMITS

GATK4_IMAGE = "community.wave.seqera.io/library/gatk4_gcnvkernel:edb12e4f0bf02cd3"
# MarkDuplicates' image adds samtools, for the CRAM conversion.
GATK4_SAMTOOLS_IMAGE = "community.wave.seqera.io/library/gatk4_gcnvkernel_htslib_samtools:d3becb6465454c35"
# CNNScoreVariants needs GATK's python/CNN stack, which upstream gets from this image.
GATK4_CNN_IMAGE = "quay.io/nf-core/gatk:4.5.0.0"
# BedToIntervalList is pinned upstream to the newer, lighter GATK image.
GATK4_LITE_IMAGE = "community.wave.seqera.io/library/gatk4-lite:4.7.0.0--79918b8b2632f4b1"
# IntervalListTools' upstream image adds python, which upstream uses to rename the split lists.
GATK4_LITE_PYTHON_IMAGE = "community.wave.seqera.io/library/gatk4-lite_python:fb04467550379e2f"

DEFAULT_MEMORY_MB = 6144
DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory=f"{DEFAULT_MEMORY_MB}Mi")
MARKDUPLICATES_MEMORY_MB = 16384
MARKDUPLICATES_RESOURCES = flyte.Resources(cpu=2, memory=f"{MARKDUPLICATES_MEMORY_MB}Mi")
HAPLOTYPECALLER_CPUS = 2
HAPLOTYPECALLER_MEMORY_MB = 12288
HAPLOTYPECALLER_RESOURCES = flyte.Resources(cpu=HAPLOTYPECALLER_CPUS, memory=f"{HAPLOTYPECALLER_MEMORY_MB}Mi")
# Upstream's process_low; enough for the test data. Whole-genome samples need more memory.
MUTECT2_MEMORY_MB = 12288
MUTECT2_RESOURCES = flyte.Resources(cpu=2, memory=f"{MUTECT2_MEMORY_MB}Mi")

# GATK writes the dictionary next to the reference, named after it.
gatk4_createsequencedictionary_cmd = shell.create(
    name="gatk4_createsequencedictionary",
    image=GATK4_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"fasta": File, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        JAVA_OPTIONS="-Xmx{JAVA_HEAP_MB}M -XX:-UsePerfData"
        FA=({{inputs.fasta}})
        NAME=$(basename "${{FA[0]}}")
        ARGS={{inputs.args}}
        W=$(mktemp -d); cd "$W"
        ln -s "${{FA[0]}}" "$NAME"
        gatk --java-options "$JAVA_OPTIONS" CreateSequenceDictionary \
            --REFERENCE "$NAME" --URI "$NAME" --TMP_DIR . $ARGS
        mv ./*.dict {{outputs.results}}/
    """,
)

gatk4_intervallisttobed_cmd = shell.create(
    name="gatk4_intervallisttobed",
    image=GATK4_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"intervals": File, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        JAVA_OPTIONS="-Xmx{JAVA_HEAP_MB}M -XX:-UsePerfData"
        IN=({{inputs.intervals}})
        ARGS={{inputs.args}}
        W=$(mktemp -d); cd "$W"
        gatk --java-options "$JAVA_OPTIONS" IntervalListToBed \
            --INPUT "${{IN[0]}}" --OUTPUT {{outputs.results}}/{{inputs.prefix}}.bed --TMP_DIR . $ARGS
    """,
)

gatk4_bedtointervallist_cmd = shell.create(
    name="gatk4_bedtointervallist",
    image=GATK4_LITE_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"bed": File, "dict": File, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        JAVA_OPTIONS="-Xmx{JAVA_HEAP_MB}M -XX:-UsePerfData"
        BED=({{inputs.bed}}); D=({{inputs.dict}})
        ARGS={{inputs.args}}
        W=$(mktemp -d); cd "$W"
        gatk --java-options "$JAVA_OPTIONS" BedToIntervalList \
            --INPUT "${{BED[0]}}" --OUTPUT {{outputs.results}}/{{inputs.prefix}}.interval_list \
            --SEQUENCE_DICTIONARY "${{D[0]}}" --TMP_DIR . $ARGS
    """,
)


# As upstream: every BAM is an --INPUT, the FASTA (with its index beside it)
# is the reference, and a `.cram` prefix makes MarkDuplicates write BAM that
# samtools then converts to CRAM and indexes (.crai).
gatk4_markduplicates_cmd = shell.create(
    name="gatk4_markduplicates",
    image=GATK4_SAMTOOLS_IMAGE,
    resources=MARKDUPLICATES_RESOURCES,
    inputs={"bams": list[File], "fasta": File | None, "fai": File | None, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        JAVA_OPTIONS="-Xmx{JAVA_HEAP_MB}M -XX:-UsePerfData"
        BAMS=({{inputs.bams}})
        shopt -s nullglob; FA=({{inputs.fasta}}); FAI=({{inputs.fai}}); shopt -u nullglob
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        INPUTS=(); for b in "${{BAMS[@]}}"; do INPUTS+=(--INPUT "$b"); done
        REF=()
        if [ ${{#FA[@]}} -gt 0 ]; then
            ln -s "${{FA[0]}}" genome.fa
            if [ ${{#FAI[@]}} -gt 0 ]; then ln -s "${{FAI[0]}}" genome.fa.fai; fi
            REF=(--REFERENCE_SEQUENCE genome.fa)
        fi
        if [[ "$P" == *.cram ]]; then PBAM="${{P%.cram}}.bam"; else PBAM="$OUT/$P"; fi
        gatk --java-options "$JAVA_OPTIONS" MarkDuplicates \
            "${{INPUTS[@]}}" --OUTPUT "$PBAM" --METRICS_FILE "$OUT/$P.metrics" \
            --TMP_DIR . "${{REF[@]}}" $ARGS
        if [[ "$P" == *.cram ]]; then
            samtools view -Ch -T genome.fa -o "$OUT/$P" "$PBAM"
            samtools index "$OUT/$P"
        fi
    """,
)


# Links the reference (FASTA + .fai + .dict, named so GATK pairs them), the
# alignment with its index, and each known-sites VCF with its .tbi into the
# working directory. Expects REF_FA/REF_FAI/REF_DICT, ALN/ALN_IDX and the
# SITES/SITES_TBI arrays; sets KNOWN_SITES and ALN_LINK.
LINK_INPUTS = r"""
        ln -s "$REF_FA" genome.fa; ln -s "$REF_FAI" genome.fa.fai; ln -s "$REF_DICT" genome.dict
        ALN_LINK=$(basename "$ALN"); IDX_EXT=bai; [[ "$ALN_LINK" == *.cram ]] && IDX_EXT=crai
        ln -s "$ALN" "$ALN_LINK"; ln -s "$ALN_IDX" "$ALN_LINK.$IDX_EXT"
        KNOWN_SITES=()
        for v in "${SITES[@]}"; do ln -s "$v" "$(basename "$v")"; KNOWN_SITES+=(--known-sites "$(basename "$v")"); done
        for t in "${SITES_TBI[@]}"; do ln -s "$t" "$(basename "$t")"; done
"""

gatk4_baserecalibrator_cmd = shell.create(
    name="gatk4_baserecalibrator",
    image=GATK4_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "alignment": File,
        "index": File,
        "intervals": File | None,
        "fasta": File,
        "fai": File,
        "dict": File,
        "known_sites": list[File],
        "known_sites_tbi": list[File],
        "prefix": str,
        "args": str,
    },
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        JAVA_OPTIONS="-Xmx{JAVA_HEAP_MB}M -XX:-UsePerfData"
        A=({{inputs.alignment}}); I=({{inputs.index}}); FA=({{inputs.fasta}}); FAI=({{inputs.fai}}); D=({{inputs.dict}})
        shopt -s nullglob
        BED=({{inputs.intervals}}); SITES=({{inputs.known_sites}}); SITES_TBI=({{inputs.known_sites_tbi}})
        shopt -u nullglob
        REF_FA="${{FA[0]}}"; REF_FAI="${{FAI[0]}}"; REF_DICT="${{D[0]}}"; ALN="${{A[0]}}"; ALN_IDX="${{I[0]}}"
        ARGS={{inputs.args}}
        W=$(mktemp -d); cd "$W"
        {LINK_INPUTS}
        INTERVALS=(); [ ${{#BED[@]}} -gt 0 ] && INTERVALS=(--intervals "${{BED[0]}}")
        gatk --java-options "$JAVA_OPTIONS" BaseRecalibrator \
            --input "$ALN_LINK" --output {{outputs.results}}/{{inputs.prefix}}.table --reference genome.fa \
            "${{INTERVALS[@]}}" "${{KNOWN_SITES[@]}}" --tmp-dir . $ARGS
    """,
)

gatk4_gatherbqsrreports_cmd = shell.create(
    name="gatk4_gatherbqsrreports",
    image=GATK4_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"tables": list[File], "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        JAVA_OPTIONS="-Xmx{JAVA_HEAP_MB}M -XX:-UsePerfData"
        TABLES=({{inputs.tables}})
        ARGS={{inputs.args}}
        INPUTS=(); for t in "${{TABLES[@]}}"; do INPUTS+=(--input "$t"); done
        W=$(mktemp -d); cd "$W"
        gatk --java-options "$JAVA_OPTIONS" GatherBQSRReports \
            "${{INPUTS[@]}}" --output {{outputs.results}}/{{inputs.prefix}}.table --tmp-dir . $ARGS
    """,
)

gatk4_applybqsr_cmd = shell.create(
    name="gatk4_applybqsr",
    image=GATK4_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "alignment": File,
        "index": File,
        "table": File,
        "intervals": File | None,
        "fasta": File,
        "fai": File,
        "dict": File,
        "prefix": str,
        "suffix": str,
        "args": str,
    },
    defaults={"suffix": "cram", "args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        JAVA_OPTIONS="-Xmx{JAVA_HEAP_MB}M -XX:-UsePerfData"
        A=({{inputs.alignment}}); I=({{inputs.index}}); T=({{inputs.table}})
        FA=({{inputs.fasta}}); FAI=({{inputs.fai}}); D=({{inputs.dict}})
        shopt -s nullglob; BED=({{inputs.intervals}}); shopt -u nullglob
        REF_FA="${{FA[0]}}"; REF_FAI="${{FAI[0]}}"; REF_DICT="${{D[0]}}"; ALN="${{A[0]}}"; ALN_IDX="${{I[0]}}"
        SITES=(); SITES_TBI=()
        ARGS={{inputs.args}}
        W=$(mktemp -d); cd "$W"
        {LINK_INPUTS}
        INTERVALS=(); [ ${{#BED[@]}} -gt 0 ] && INTERVALS=(--intervals "${{BED[0]}}")
        gatk --java-options "$JAVA_OPTIONS" ApplyBQSR \
            --input "$ALN_LINK" --output {{outputs.results}}/{{inputs.prefix}}.{{inputs.suffix}} --reference genome.fa \
            --bqsr-recal-file "${{T[0]}}" "${{INTERVALS[@]}}" --tmp-dir . $ARGS
    """,
)


# As upstream: --dbsnp annotates IDs, HMM threads follow the task's CPUs, and a
# .vcf.gz output gets its .tbi from GATK.
gatk4_haplotypecaller_cmd = shell.create(
    name="gatk4_haplotypecaller",
    image=GATK4_IMAGE,
    resources=HAPLOTYPECALLER_RESOURCES,
    inputs={
        "alignment": File,
        "index": File,
        "intervals": File | None,
        "fasta": File,
        "fai": File,
        "dict": File,
        "dbsnp": File | None,
        "dbsnp_tbi": File | None,
        "prefix": str,
        "args": str,
    },
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        JAVA_OPTIONS="-Xmx{JAVA_HEAP_MB}M -XX:-UsePerfData"
        A=({{inputs.alignment}}); I=({{inputs.index}}); FA=({{inputs.fasta}}); FAI=({{inputs.fai}}); D=({{inputs.dict}})
        shopt -s nullglob
        BED=({{inputs.intervals}}); SITES=({{inputs.dbsnp}}); SITES_TBI=({{inputs.dbsnp_tbi}})
        shopt -u nullglob
        REF_FA="${{FA[0]}}"; REF_FAI="${{FAI[0]}}"; REF_DICT="${{D[0]}}"; ALN="${{A[0]}}"; ALN_IDX="${{I[0]}}"
        ARGS={{inputs.args}}
        W=$(mktemp -d); cd "$W"
        {LINK_INPUTS}
        DBSNP=(); [ ${{#SITES[@]}} -gt 0 ] && DBSNP=(--dbsnp "$(basename "${{SITES[0]}}")")
        INTERVALS=(); [ ${{#BED[@]}} -gt 0 ] && INTERVALS=(--intervals "${{BED[0]}}")
        gatk --java-options "$JAVA_OPTIONS" HaplotypeCaller \
            --input "$ALN_LINK" --output {{outputs.results}}/{{inputs.prefix}}.vcf.gz --reference genome.fa \
            --native-pair-hmm-threads $TASK_CPUS "${{DBSNP[@]}}" "${{INTERVALS[@]}}" --tmp-dir . $ARGS
    """,
)

gatk4_mergevcfs_cmd = shell.create(
    name="gatk4_mergevcfs",
    image=GATK4_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"vcfs": list[File], "dict": File | None, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        JAVA_OPTIONS="-Xmx{JAVA_HEAP_MB}M -XX:-UsePerfData"
        VCFS=({{inputs.vcfs}})
        shopt -s nullglob; D=({{inputs.dict}}); shopt -u nullglob
        ARGS={{inputs.args}}
        INPUTS=(); for v in "${{VCFS[@]}}"; do INPUTS+=(--INPUT "$v"); done
        DICT=(); [ ${{#D[@]}} -gt 0 ] && DICT=(--SEQUENCE_DICTIONARY "${{D[0]}}")
        W=$(mktemp -d); cd "$W"
        gatk --java-options "$JAVA_OPTIONS" MergeVcfs \
            "${{INPUTS[@]}}" --OUTPUT {{outputs.results}}/{{inputs.prefix}}.vcf.gz "${{DICT[@]}}" --TMP_DIR . $ARGS
    """,
)

# Links a VCF (VCF/VCF_TBI) and the reference (REF_*) into the working
# directory side by side; sets VCF_LINK.
LINK_VCF = r"""
        ln -s "$REF_FA" genome.fa; ln -s "$REF_FAI" genome.fa.fai; ln -s "$REF_DICT" genome.dict
        VCF_LINK=$(basename "$VCF"); ln -s "$VCF" "$VCF_LINK"; ln -s "$VCF_TBI" "$VCF_LINK.tbi"
"""

gatk4_cnnscorevariants_cmd = shell.create(
    name="gatk4_cnnscorevariants",
    image=GATK4_CNN_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "vcf": File,
        "tbi": File,
        "intervals": File | None,
        "fasta": File,
        "fai": File,
        "dict": File,
        "prefix": str,
        "args": str,
    },
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        JAVA_OPTIONS="-Xmx{JAVA_HEAP_MB}M -XX:-UsePerfData"
        V=({{inputs.vcf}}); T=({{inputs.tbi}}); FA=({{inputs.fasta}}); FAI=({{inputs.fai}}); D=({{inputs.dict}})
        shopt -s nullglob; BED=({{inputs.intervals}}); shopt -u nullglob
        REF_FA="${{FA[0]}}"; REF_FAI="${{FAI[0]}}"; REF_DICT="${{D[0]}}"; VCF="${{V[0]}}"; VCF_TBI="${{T[0]}}"
        ARGS={{inputs.args}}
        W=$(mktemp -d); cd "$W"
        {LINK_VCF}
        INTERVALS=(); [ ${{#BED[@]}} -gt 0 ] && INTERVALS=(--intervals "${{BED[0]}}")
        export THEANO_FLAGS="base_compiledir=$PWD"
        gatk --java-options "$JAVA_OPTIONS" CNNScoreVariants \
            --variant "$VCF_LINK" --output {{outputs.results}}/{{inputs.prefix}}.cnn.vcf.gz --reference genome.fa \
            "${{INTERVALS[@]}}" --tmp-dir . $ARGS
    """,
)

# The resources (known sites) are linked beside their .tbi files.
gatk4_filtervarianttranches_cmd = shell.create(
    name="gatk4_filtervarianttranches",
    image=GATK4_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "vcf": File,
        "tbi": File,
        "resources": list[File],
        "resources_tbi": list[File],
        "fasta": File,
        "fai": File,
        "dict": File,
        "prefix": str,
        "args": str,
    },
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        JAVA_OPTIONS="-Xmx{JAVA_HEAP_MB}M -XX:-UsePerfData"
        V=({{inputs.vcf}}); T=({{inputs.tbi}}); FA=({{inputs.fasta}}); FAI=({{inputs.fai}}); D=({{inputs.dict}})
        RES=({{inputs.resources}}); RES_TBI=({{inputs.resources_tbi}})
        REF_FA="${{FA[0]}}"; REF_FAI="${{FAI[0]}}"; REF_DICT="${{D[0]}}"; VCF="${{V[0]}}"; VCF_TBI="${{T[0]}}"
        ARGS={{inputs.args}}
        W=$(mktemp -d); cd "$W"
        {LINK_VCF}
        RESOURCES=()
        for r in "${{RES[@]}}"; do ln -s "$r" "$(basename "$r")"; RESOURCES+=(--resource "$(basename "$r")"); done
        for t in "${{RES_TBI[@]}}"; do ln -s "$t" "$(basename "$t")"; done
        gatk --java-options "$JAVA_OPTIONS" FilterVariantTranches \
            --variant "$VCF_LINK" "${{RESOURCES[@]}}" --output {{outputs.results}}/{{inputs.prefix}}.filtered.vcf.gz \
            --tmp-dir . $ARGS
    """,
)


# Links the reference (REF_FA, REF_FAI, REF_DICT, and REF_GZI, empty for a
# plain FASTA) under the FASTA's own name, with the dictionary named as GATK
# looks for it (`genome.dict` for `genome.fa` or `genome.fasta.gz`); sets REF.
# Defines link_vcf VCF [INDEX]: links a VCF under its own name, with its
# index (.tbi or .idx) beside it.
LINK_REFERENCE = r"""
        REF=$(basename "$REF_FA"); ln -s "$REF_FA" "$REF"; ln -s "$REF_FAI" "$REF.fai"
        if [ -n "$REF_GZI" ]; then ln -s "$REF_GZI" "$REF.gzi"; fi
        REF_STEM="${REF%.gz}"; ln -s "$REF_DICT" "${REF_STEM%.*}.dict"
        link_vcf() {
            ln -s "$1" "$(basename "$1")"
            if [ -n "${2:-}" ]; then ln -s "$2" "$(basename "$1").${2##*.}"; fi
        }
"""

# As upstream: every alignment is an --input (a matched normal is named with
# --normal-sample in args), and the panel of normals, germline resource and
# force-called alleles are optional. Alignments and their indices are linked
# under their own names, so GATK pairs them. An f1r2 archive that args ask for
# (--f1r2-tar-gz <prefix>.f1r2.tar.gz) is kept.
gatk4_mutect2_cmd = shell.create(
    name="gatk4_mutect2",
    image=GATK4_LITE_IMAGE,
    resources=MUTECT2_RESOURCES,
    inputs={
        "alignments": list[File],
        "indices": list[File],
        "intervals": File | None,
        "fasta": File,
        "fai": File,
        "gzi": File | None,
        "dict": File,
        "alleles": File | None,
        "alleles_tbi": File | None,
        "germline_resource": File | None,
        "germline_resource_tbi": File | None,
        "panel_of_normals": File | None,
        "panel_of_normals_tbi": File | None,
        "prefix": str,
        "args": str,
    },
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        JAVA_OPTIONS="-Xmx{JAVA_HEAP_MB}M -XX:-UsePerfData"
        ALNS=({{inputs.alignments}}); IDXS=({{inputs.indices}})
        FA=({{inputs.fasta}}); FAI=({{inputs.fai}}); D=({{inputs.dict}})
        shopt -s nullglob
        GZI=({{inputs.gzi}}); BED=({{inputs.intervals}})
        AL=({{inputs.alleles}}); AL_TBI=({{inputs.alleles_tbi}})
        GR=({{inputs.germline_resource}}); GR_TBI=({{inputs.germline_resource_tbi}})
        PON=({{inputs.panel_of_normals}}); PON_TBI=({{inputs.panel_of_normals_tbi}})
        shopt -u nullglob
        REF_FA="${{FA[0]}}"; REF_FAI="${{FAI[0]}}"; REF_DICT="${{D[0]}}"; REF_GZI="${{GZI[0]:-}}"
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        {LINK_REFERENCE}
        INPUTS=()
        for a in "${{ALNS[@]}}"; do ln -s "$a" "$(basename "$a")"; INPUTS+=(--input "$(basename "$a")"); done
        for i in "${{IDXS[@]}}"; do ln -s "$i" "$(basename "$i")"; done
        OPTIONAL=()
        if [ ${{#PON[@]}} -gt 0 ]; then
            link_vcf "${{PON[0]}}" "${{PON_TBI[0]:-}}"; OPTIONAL+=(--panel-of-normals "$(basename "${{PON[0]}}")")
        fi
        if [ ${{#GR[@]}} -gt 0 ]; then
            link_vcf "${{GR[0]}}" "${{GR_TBI[0]:-}}"; OPTIONAL+=(--germline-resource "$(basename "${{GR[0]}}")")
        fi
        if [ ${{#AL[@]}} -gt 0 ]; then
            link_vcf "${{AL[0]}}" "${{AL_TBI[0]:-}}"; OPTIONAL+=(--alleles "$(basename "${{AL[0]}}")")
        fi
        if [ ${{#BED[@]}} -gt 0 ]; then
            ln -s "${{BED[0]}}" "$(basename "${{BED[0]}}")"; OPTIONAL+=(--intervals "$(basename "${{BED[0]}}")")
        fi
        gatk --java-options "$JAVA_OPTIONS" Mutect2 \
            "${{INPUTS[@]}}" --output "$OUT/$P.vcf.gz" --reference "$REF" \
            "${{OPTIONAL[@]}}" --tmp-dir . $ARGS
        if [ -e "$P.f1r2.tar.gz" ]; then mv "$P.f1r2.tar.gz" "$OUT/"; fi
    """,
)

gatk4_mergemutectstats_cmd = shell.create(
    name="gatk4_mergemutectstats",
    image=GATK4_LITE_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"stats": list[File], "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        JAVA_OPTIONS="-Xmx{JAVA_HEAP_MB}M -XX:-UsePerfData"
        STATS=({{inputs.stats}})
        ARGS={{inputs.args}}
        INPUTS=(); for s in "${{STATS[@]}}"; do INPUTS+=(--stats "$s"); done
        W=$(mktemp -d); cd "$W"
        gatk --java-options "$JAVA_OPTIONS" MergeMutectStats \
            "${{INPUTS[@]}}" --output {{outputs.results}}/{{inputs.prefix}}.vcf.gz.stats --tmp-dir . $ARGS
    """,
)

# As upstream: the stats are passed with --stats (without them GATK looks for
# <vcf>.stats beside the VCF, which isn't there, and fails), each orientation
# prior, segmentation and contamination table is repeated, and a contamination
# estimate (a fraction, as text) is passed when set.
gatk4_filtermutectcalls_cmd = shell.create(
    name="gatk4_filtermutectcalls",
    image=GATK4_LITE_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "vcf": File,
        "tbi": File,
        "stats": File | None,
        "orientation_priors": list[File],
        "segmentations": list[File],
        "contamination_tables": list[File],
        "contamination_estimate": str,
        "fasta": File,
        "fai": File,
        "dict": File,
        "prefix": str,
        "args": str,
    },
    defaults={"contamination_estimate": "", "args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        JAVA_OPTIONS="-Xmx{JAVA_HEAP_MB}M -XX:-UsePerfData"
        V=({{inputs.vcf}}); T=({{inputs.tbi}}); FA=({{inputs.fasta}}); FAI=({{inputs.fai}}); D=({{inputs.dict}})
        shopt -s nullglob
        STATS=({{inputs.stats}}); PRIORS=({{inputs.orientation_priors}})
        SEGMENTS=({{inputs.segmentations}}); TABLES=({{inputs.contamination_tables}})
        shopt -u nullglob
        REF_FA="${{FA[0]}}"; REF_FAI="${{FAI[0]}}"; REF_DICT="${{D[0]}}"; REF_GZI=""
        ESTIMATE={{inputs.contamination_estimate}}
        ARGS={{inputs.args}}
        W=$(mktemp -d); cd "$W"
        {LINK_REFERENCE}
        link_vcf "${{V[0]}}" "${{T[0]}}"
        OPTIONAL=()
        if [ ${{#STATS[@]}} -gt 0 ]; then OPTIONAL+=(--stats "${{STATS[0]}}"); fi
        for f in "${{PRIORS[@]}}"; do OPTIONAL+=(--orientation-bias-artifact-priors "$f"); done
        for f in "${{SEGMENTS[@]}}"; do OPTIONAL+=(--tumor-segmentation "$f"); done
        if [ -n "$ESTIMATE" ]; then OPTIONAL+=(--contamination-estimate "$ESTIMATE"); fi
        for f in "${{TABLES[@]}}"; do OPTIONAL+=(--contamination-table "$f"); done
        gatk --java-options "$JAVA_OPTIONS" FilterMutectCalls \
            --variant "$(basename "${{V[0]}}")" --output {{outputs.results}}/{{inputs.prefix}}.vcf.gz \
            --reference "$REF" "${{OPTIONAL[@]}}" --tmp-dir . $ARGS
    """,
)

# The VCF's index may be a .tbi or, for a plain VCF, a .idx. Filters (name to
# JEXL expression) are their own input because expressions contain spaces.
gatk4_variantfiltration_cmd = shell.create(
    name="gatk4_variantfiltration",
    image=GATK4_LITE_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "vcf": File,
        "index": File,
        "fasta": File,
        "fai": File,
        "gzi": File | None,
        "dict": File,
        "filters": dict[str, str],
        "prefix": str,
        "args": str,
    },
    defaults={"filters": {}, "args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        JAVA_OPTIONS="-Xmx{JAVA_HEAP_MB}M -XX:-UsePerfData"
        V=({{inputs.vcf}}); I=({{inputs.index}}); FA=({{inputs.fasta}}); FAI=({{inputs.fai}}); D=({{inputs.dict}})
        shopt -s nullglob; GZI=({{inputs.gzi}}); shopt -u nullglob
        REF_FA="${{FA[0]}}"; REF_FAI="${{FAI[0]}}"; REF_DICT="${{D[0]}}"; REF_GZI="${{GZI[0]:-}}"
        ARGS={{inputs.args}}
        W=$(mktemp -d); cd "$W"
        {LINK_REFERENCE}
        link_vcf "${{V[0]}}" "${{I[0]}}"
        PAIRS=({{inputs.filters}}); FILTERS=()
        for ((k = 0; k < ${{#PAIRS[@]}}; k += 2)); do
            FILTERS+=(--filter-name "${{PAIRS[k]}}" --filter-expression "${{PAIRS[k + 1]}}")
        done
        gatk --java-options "$JAVA_OPTIONS" VariantFiltration \
            --variant "$(basename "${{V[0]}}")" --output {{outputs.results}}/{{inputs.prefix}}.vcf.gz \
            --reference "$REF" "${{FILTERS[@]}}" --tmp-dir . $ARGS
    """,
)

# As upstream: the lists go to <prefix>_split/temp_*/, and each is renamed
# with its 1-based position in path order (`1scattered.interval_list`, ...) so
# their names stay distinct once collected. The lists' @PG header records the
# command line, so the input and output are named relative to the working
# directory, as upstream's are.
gatk4_intervallisttools_cmd = shell.create(
    name="gatk4_intervallisttools",
    image=GATK4_LITE_PYTHON_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"intervals": File, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        JAVA_OPTIONS="-Xmx{JAVA_HEAP_MB}M -XX:-UsePerfData"
        IN=({{inputs.intervals}})
        ARGS={{inputs.args}}
        SPLIT={{inputs.prefix}}_split
        W=$(mktemp -d); cd "$W"
        ln -s "${{IN[0]}}" "$(basename "${{IN[0]}}")"
        mkdir "$SPLIT"
        gatk --java-options "$JAVA_OPTIONS" IntervalListTools \
            --INPUT "$(basename "${{IN[0]}}")" --OUTPUT "$SPLIT" --TMP_DIR . $ARGS
        N=0
        for f in $(printf '%s\n' "$SPLIT"/*/*.interval_list | LC_ALL=C sort); do
            N=$((N + 1)); mv "$f" "$(dirname "$f")/$N$(basename "$f")"
        done
        mv "$SPLIT" {{outputs.results}}/
    """,
)

# One image per environment: MarkDuplicates', CNNScoreVariants', the lite tools' and IntervalListTools' differ,
# so `env` has five.
gatk4_env = flyte.TaskEnvironment.from_task(
    "gatk4_tools",
    gatk4_createsequencedictionary_cmd.as_task(),
    gatk4_intervallisttobed_cmd.as_task(),
    gatk4_baserecalibrator_cmd.as_task(),
    gatk4_gatherbqsrreports_cmd.as_task(),
    gatk4_applybqsr_cmd.as_task(),
    gatk4_haplotypecaller_cmd.as_task(),
    gatk4_mergevcfs_cmd.as_task(),
    gatk4_filtervarianttranches_cmd.as_task(),
)
gatk4_samtools_env = flyte.TaskEnvironment.from_task("gatk4_samtools", gatk4_markduplicates_cmd.as_task())
gatk4_cnn_env = flyte.TaskEnvironment.from_task("gatk4_cnn", gatk4_cnnscorevariants_cmd.as_task())
gatk4_lite_env = flyte.TaskEnvironment.from_task(
    "gatk4_lite",
    gatk4_bedtointervallist_cmd.as_task(),
    gatk4_mutect2_cmd.as_task(),
    gatk4_mergemutectstats_cmd.as_task(),
    gatk4_filtermutectcalls_cmd.as_task(),
    gatk4_variantfiltration_cmd.as_task(),
)
gatk4_lite_python_env = flyte.TaskEnvironment.from_task("gatk4_lite_python", gatk4_intervallisttools_cmd.as_task())
env = flyte.TaskEnvironment(
    name="gatk4",
    depends_on=[gatk4_env, gatk4_samtools_env, gatk4_cnn_env, gatk4_lite_env, gatk4_lite_python_env],
)


async def only_file(results: Dir, suffix: str, tool: str) -> File:
    found = [f async for f in results.walk() if f.path.endswith(suffix)]
    if len(found) != 1:
        raise FileNotFoundError(f"GATK {tool} wrote {len(found)} *{suffix} files, expected one")
    return found[0]


async def gatk4_createsequencedictionary(fasta: File, args: str = "") -> File:
    """The sequence dictionary of ``fasta``: ``<fasta name without extension>.dict``."""
    results = await gatk4_createsequencedictionary_cmd(fasta=fasta, args=args)
    return await only_file(results, ".dict", "CreateSequenceDictionary")


async def gatk4_intervallisttobed(intervals: File, prefix: str, args: str = "") -> File:
    """``intervals`` (a Picard ``.interval_list``) as ``<prefix>.bed``."""
    results = await gatk4_intervallisttobed_cmd(intervals=intervals, prefix=prefix, args=args)
    return await only_file(results, f"{prefix}.bed", "IntervalListToBed")


async def gatk4_bedtointervallist(bed: File, dict: File, prefix: str, args: str = "") -> File:
    """``bed`` as a Picard interval list, ``<prefix>.interval_list``, headed by ``dict``."""
    results = await gatk4_bedtointervallist_cmd(bed=bed, dict=dict, prefix=prefix, args=args)
    return await only_file(results, f"{prefix}.interval_list", "BedToIntervalList")


@dataclass
class MarkDuplicatesResult:
    alignment: File  # <prefix> (.cram or .bam)
    index: File | None  # .crai for CRAM; BAM output is indexed only if args ask (--CREATE_INDEX true)
    metrics: File  # <prefix>.metrics


async def gatk4_markduplicates(
    bams: list[File],
    prefix: str,
    fasta: File | None = None,
    fai: File | None = None,
    args: str = "",
) -> MarkDuplicatesResult:
    """Mark duplicates across ``bams`` into ``<prefix>``: CRAM when it ends in ``.cram`` (needs ``fasta``).

    BAMs are staged under their own names, so they must be distinct.
    """
    if prefix.endswith(".cram") and fasta is None:
        raise ValueError("CRAM output needs the reference fasta")
    results = await gatk4_markduplicates_cmd(bams=bams, fasta=fasta, fai=fai, prefix=prefix, args=args)
    alignment = await results.get_file(prefix)
    metrics = await results.get_file(f"{prefix}.metrics")
    if alignment is None or metrics is None:
        raise FileNotFoundError(f"MarkDuplicates wrote no {prefix} or {prefix}.metrics")
    index = await results.get_file(f"{prefix}.crai") if prefix.endswith(".cram") else None
    return MarkDuplicatesResult(alignment=alignment, index=index, metrics=metrics)


async def gatk4_baserecalibrator(
    alignment: File,
    index: File,
    fasta: File,
    fai: File,
    dict: File,
    known_sites: list[File],
    known_sites_tbi: list[File],
    prefix: str,
    intervals: File | None = None,
    args: str = "",
) -> File:
    """``<prefix>.table``: recalibration covariates of ``alignment`` (BAM/CRAM) against the known sites."""
    results = await gatk4_baserecalibrator_cmd(
        alignment=alignment,
        index=index,
        intervals=intervals,
        fasta=fasta,
        fai=fai,
        dict=dict,
        known_sites=known_sites,
        known_sites_tbi=known_sites_tbi,
        prefix=prefix,
        args=args,
    )
    return await only_file(results, f"{prefix}.table", "BaseRecalibrator")


async def gatk4_gatherbqsrreports(tables: list[File], prefix: str, args: str = "") -> File:
    """``<prefix>.table``: per-interval recalibration tables merged into one.

    Tables are staged under their own names, so they must be distinct.
    """
    results = await gatk4_gatherbqsrreports_cmd(tables=tables, prefix=prefix, args=args)
    return await only_file(results, f"{prefix}.table", "GatherBQSRReports")


async def gatk4_applybqsr(
    alignment: File,
    index: File,
    table: File,
    fasta: File,
    fai: File,
    dict: File,
    prefix: str,
    intervals: File | None = None,
    suffix: str = "cram",
    args: str = "",
) -> File:
    """``<prefix>.<suffix>``: ``alignment`` with recalibrated base qualities (CRAM, or BAM with ``suffix="bam"``)."""
    if suffix not in ("cram", "bam"):
        raise ValueError(f"suffix must be cram or bam, got {suffix!r}")
    results = await gatk4_applybqsr_cmd(
        alignment=alignment,
        index=index,
        table=table,
        intervals=intervals,
        fasta=fasta,
        fai=fai,
        dict=dict,
        prefix=prefix,
        suffix=suffix,
        args=args,
    )
    return await only_file(results, f"{prefix}.{suffix}", "ApplyBQSR")


@dataclass
class IndexedVcf:
    vcf: File  # .vcf.gz
    tbi: File


async def vcf_and_index(results: Dir, name: str, tool: str) -> IndexedVcf:
    vcf, tbi = await results.get_file(name), await results.get_file(f"{name}.tbi")
    if vcf is None or tbi is None:
        raise FileNotFoundError(f"GATK {tool} wrote no {name} or {name}.tbi")
    return IndexedVcf(vcf=vcf, tbi=tbi)


async def gatk4_haplotypecaller(
    alignment: File,
    index: File,
    fasta: File,
    fai: File,
    dict: File,
    prefix: str,
    intervals: File | None = None,
    dbsnp: File | None = None,
    dbsnp_tbi: File | None = None,
    args: str = "",
) -> IndexedVcf:
    """Germline calls from ``alignment`` (BAM/CRAM): ``<prefix>.vcf.gz`` and its ``.tbi``."""
    results = await gatk4_haplotypecaller_cmd(
        alignment=alignment,
        index=index,
        intervals=intervals,
        fasta=fasta,
        fai=fai,
        dict=dict,
        dbsnp=dbsnp,
        dbsnp_tbi=dbsnp_tbi,
        prefix=prefix,
        args=args,
    )
    return await vcf_and_index(results, f"{prefix}.vcf.gz", "HaplotypeCaller")


async def gatk4_mergevcfs(vcfs: list[File], prefix: str, dict: File | None = None, args: str = "") -> IndexedVcf:
    """``<prefix>.vcf.gz`` (+ ``.tbi``): ``vcfs`` merged, ordered by ``dict`` when given.

    VCFs are staged under their own names, so they must be distinct.
    """
    results = await gatk4_mergevcfs_cmd(vcfs=vcfs, dict=dict, prefix=prefix, args=args)
    return await vcf_and_index(results, f"{prefix}.vcf.gz", "MergeVcfs")


async def gatk4_cnnscorevariants(
    vcf: File,
    tbi: File,
    fasta: File,
    fai: File,
    dict: File,
    prefix: str,
    intervals: File | None = None,
    args: str = "",
) -> IndexedVcf:
    """``<prefix>.cnn.vcf.gz`` (+ ``.tbi``): ``vcf`` with CNN scores (INFO ``CNN_1D`` by default)."""
    results = await gatk4_cnnscorevariants_cmd(
        vcf=vcf, tbi=tbi, intervals=intervals, fasta=fasta, fai=fai, dict=dict, prefix=prefix, args=args
    )
    return await vcf_and_index(results, f"{prefix}.cnn.vcf.gz", "CNNScoreVariants")


async def gatk4_filtervarianttranches(
    vcf: File,
    tbi: File,
    resources: list[File],
    resources_tbi: list[File],
    fasta: File,
    fai: File,
    dict: File,
    prefix: str,
    args: str = "",
) -> IndexedVcf:
    """``<prefix>.filtered.vcf.gz`` (+ ``.tbi``): CNN-scored ``vcf`` filtered by tranches of the resources' sites.

    ``args`` must name the score, e.g. ``--info-key CNN_1D``.
    """
    results = await gatk4_filtervarianttranches_cmd(
        vcf=vcf,
        tbi=tbi,
        resources=resources,
        resources_tbi=resources_tbi,
        fasta=fasta,
        fai=fai,
        dict=dict,
        prefix=prefix,
        args=args,
    )
    return await vcf_and_index(results, f"{prefix}.filtered.vcf.gz", "FilterVariantTranches")


@dataclass
class Mutect2Result:
    vcf: File  # <prefix>.vcf.gz
    tbi: File  # <prefix>.vcf.gz.tbi
    stats: File  # <prefix>.vcf.gz.stats, for FilterMutectCalls
    f1r2: File | None  # <prefix>.f1r2.tar.gz, when args pass --f1r2-tar-gz <prefix>.f1r2.tar.gz


async def gatk4_mutect2(
    alignments: list[File],
    indices: list[File],
    fasta: File,
    fai: File,
    dict: File,
    prefix: str,
    intervals: File | None = None,
    gzi: File | None = None,
    alleles: File | None = None,
    alleles_tbi: File | None = None,
    germline_resource: File | None = None,
    germline_resource_tbi: File | None = None,
    panel_of_normals: File | None = None,
    panel_of_normals_tbi: File | None = None,
    args: str = "",
) -> Mutect2Result:
    """Somatic calls from ``alignments`` (BAM/CRAM): ``<prefix>.vcf.gz`` with its ``.tbi`` and ``.stats``.

    A tumour-normal pair passes both alignments and names the normal's sample
    in ``args`` (``--normal-sample <SM>``). Indices are paired with their
    alignments by name (``x.bam.bai`` or ``x.bai``, ``x.cram.crai``), so all
    the names must be distinct. ``gzi`` is the ``.gzi`` of a bgzipped
    ``fasta``.
    """
    if len(indices) != len(alignments):
        raise ValueError(f"{len(alignments)} alignments but {len(indices)} indices")
    results = await gatk4_mutect2_cmd(
        alignments=alignments,
        indices=indices,
        intervals=intervals,
        fasta=fasta,
        fai=fai,
        gzi=gzi,
        dict=dict,
        alleles=alleles,
        alleles_tbi=alleles_tbi,
        germline_resource=germline_resource,
        germline_resource_tbi=germline_resource_tbi,
        panel_of_normals=panel_of_normals,
        panel_of_normals_tbi=panel_of_normals_tbi,
        prefix=prefix,
        args=args,
    )
    calls = await vcf_and_index(results, f"{prefix}.vcf.gz", "Mutect2")
    stats = await results.get_file(f"{prefix}.vcf.gz.stats")
    if stats is None:
        raise FileNotFoundError(f"GATK Mutect2 wrote no {prefix}.vcf.gz.stats")
    f1r2 = await results.get_file(f"{prefix}.f1r2.tar.gz")
    return Mutect2Result(vcf=calls.vcf, tbi=calls.tbi, stats=stats, f1r2=f1r2)


async def gatk4_mergemutectstats(stats: list[File], prefix: str, args: str = "") -> File:
    """``<prefix>.vcf.gz.stats``: per-interval Mutect2 stats merged into one.

    Stats are staged under their own names, so they must be distinct.
    """
    results = await gatk4_mergemutectstats_cmd(stats=stats, prefix=prefix, args=args)
    return await only_file(results, f"{prefix}.vcf.gz.stats", "MergeMutectStats")


@dataclass
class FilterMutectCallsResult:
    vcf: File  # <prefix>.vcf.gz
    tbi: File  # <prefix>.vcf.gz.tbi
    stats: File  # <prefix>.vcf.gz.filteringStats.tsv


async def gatk4_filtermutectcalls(
    vcf: File,
    tbi: File,
    stats: File | None,
    fasta: File,
    fai: File,
    dict: File,
    prefix: str,
    orientation_priors: list[File] | None = None,
    segmentations: list[File] | None = None,
    contamination_tables: list[File] | None = None,
    contamination_estimate: float | None = None,
    args: str = "",
) -> FilterMutectCallsResult:
    """Mutect2 calls in ``vcf`` with filters applied: ``<prefix>.vcf.gz`` (+ ``.tbi``) and its filtering stats.

    ``stats`` is Mutect2's ``.stats`` for ``vcf``; GATK fails without it.
    ``orientation_priors`` come from LearnReadOrientationModel,
    ``segmentations`` and ``contamination_tables`` from
    CalculateContamination; ``contamination_estimate`` is a fixed fraction
    instead of a table.
    """
    results = await gatk4_filtermutectcalls_cmd(
        vcf=vcf,
        tbi=tbi,
        stats=stats,
        orientation_priors=orientation_priors or [],
        segmentations=segmentations or [],
        contamination_tables=contamination_tables or [],
        contamination_estimate="" if contamination_estimate is None else str(contamination_estimate),
        fasta=fasta,
        fai=fai,
        dict=dict,
        prefix=prefix,
        args=args,
    )
    calls = await vcf_and_index(results, f"{prefix}.vcf.gz", "FilterMutectCalls")
    filtering_stats = await results.get_file(f"{prefix}.vcf.gz.filteringStats.tsv")
    if filtering_stats is None:
        raise FileNotFoundError(f"GATK FilterMutectCalls wrote no {prefix}.vcf.gz.filteringStats.tsv")
    return FilterMutectCallsResult(vcf=calls.vcf, tbi=calls.tbi, stats=filtering_stats)


async def gatk4_variantfiltration(
    vcf: File,
    index: File,
    fasta: File,
    fai: File,
    dict: File,
    prefix: str,
    gzi: File | None = None,
    filters: dict[str, str] | None = None,
    args: str = "",
) -> IndexedVcf:
    """``<prefix>.vcf.gz`` (+ ``.tbi``): ``vcf`` with ``filters`` applied.

    ``filters`` maps each FILTER name to its JEXL expression, e.g.
    ``{"FS": "FS > 30.0"}`` (``--filter-name``/``--filter-expression``).
    ``args`` carries other options (``--window 35 --cluster 3``); it is split
    on whitespace, so it can't hold expressions. ``index`` is the VCF's
    ``.tbi`` (or ``.idx`` for a plain VCF); ``gzi`` is the ``.gzi`` of a
    bgzipped ``fasta``.
    """
    results = await gatk4_variantfiltration_cmd(
        vcf=vcf,
        index=index,
        fasta=fasta,
        fai=fai,
        gzi=gzi,
        dict=dict,
        filters=filters or {},
        prefix=prefix,
        args=args,
    )
    return await vcf_and_index(results, f"{prefix}.vcf.gz", "VariantFiltration")


async def gatk4_intervallisttools(intervals: File, prefix: str, args: str = "") -> list[File]:
    """``intervals`` split as ``args`` say (e.g. ``--SCATTER_COUNT 6``): the lists, in order.

    Each is named ``<n>scattered.interval_list``, numbered from 1 in path order
    (the order returned).
    """
    results = await gatk4_intervallisttools_cmd(intervals=intervals, prefix=prefix, args=args)
    lists = [f async for f in results.walk() if f.path.endswith(".interval_list")]
    return sorted(lists, key=lambda f: (len(f.name), f.name))
