"""parabricks — NVIDIA Parabricks, GPU-accelerated alignment and variant calling.

Exposes:

- :func:`parabricks_fq2bam` — align DNA reads with GPU BWA-MEM, sort and mark
  duplicates (``pbrun fq2bam``), optionally with a BQSR report from known
  sites: ``<prefix>.bam`` + ``.bai`` (or ``.cram`` + ``.crai``), and
  ``<prefix>.table`` with known sites.
- :func:`parabricks_fq2bammeth` — the same for bisulfite reads against a
  bwa-meth index (``pbrun fq2bam_meth``): ``<prefix>.bam`` + ``.bai``.
- :func:`parabricks_rnafq2bam` — align RNA reads with GPU STAR
  (``pbrun rna_fq2bam``): ``<prefix>.bam`` + ``.bai`` and STAR's
  ``<prefix>.Log.final.out``, ``.Log.out``, ``.Log.progress.out``,
  ``.SJ.out.tab`` (and chimeric junctions, counts, ... as ``args`` ask).
- :func:`parabricks_star_index` — a STAR genome index built with STAR 2.7.2a,
  the STAR version whose indexes ``rna_fq2bam`` reads (as upstream's test).
- :func:`parabricks_applybqsr` — apply a BQSR report (``pbrun applybqsr``):
  ``<prefix>.bam`` + ``.bai`` (or ``.cram`` + ``.crai``).
- :func:`parabricks_haplotypecaller` — GPU GATK HaplotypeCaller:
  ``<prefix>.vcf``, or ``<prefix>.g.vcf.gz`` with ``--gvcf`` in ``args``.
- :func:`parabricks_deepvariant` — GPU DeepVariant: ``<prefix>.vcf.gz``, or
  ``<prefix>.g.vcf.gz`` with ``--gvcf`` in ``args``.
- :func:`parabricks_deepsomatic` — GPU DeepSomatic on a tumour/normal pair:
  ``<prefix>.vcf.gz`` + ``.tbi``.

Every tool also returns its whole output directory as ``results`` (QC metrics
directories, duplicate metrics and STAR's extra outputs are picked from it).

GPU: every Parabricks task requests one NVIDIA GPU (``gpu=1``) and passes
Parabricks ``--num-gpus`` as the number of GPUs the task actually sees, as
upstream does from the task's accelerator request. A task without a GPU fails
with a message. On the test cluster ``gpu=1`` schedules on an AWS g4dn.4xlarge
(one 16 GB T4, 16 vCPUs, 64 GiB); other GPU types (``"A10G:1"``,
``"L4:1"``, ...) work wherever the cluster offers them, via
``<cmd>.override(resources=flyte.Resources(..., gpu="A10G:1"))``.

Resources: upstream labels every tool ``process_high`` (12 CPUs, 72 GB) plus
``process_gpu``. The defaults here (12 CPUs, 48 GiB, 1 GPU) fit a single-T4
node and are enough for test-sized data. Real human genomes need more host
memory: NVIDIA recommends at least 64 GB of host memory per GPU for fq2bam,
DeepVariant and HaplotypeCaller, and ``rna_fq2bam`` needs the human STAR index
in memory (about 40 GB). On 16 GB GPUs such as the T4, pass Parabricks'
low-memory options in ``args`` (``--low-memory`` for fq2bam,
``--htvc-low-memory`` for HaplotypeCaller). Raise them with
``<cmd>.override(resources=...)``.

As upstream, the BWA index and the FASTA are put side by side (fq2bam copies
the FASTA to the index's prefix, fq2bammeth into the bwa-meth index as
``<index>/<fasta name>``), and fq2bam gets the task's CPUs as
``--bwa-cpu-thread-pool``, applybqsr as ``--num-threads``. Parabricks follows
symlinks and looks for an input's index (``.bai``, ``.crai``, ``.tbi``) next
to the real file, so alignments, references and known sites are hard-linked
(or, across filesystems, copied) into the task's work directory with their
indexes, as upstream stages its inputs by copying.

Reads are given as a list in their upstream order (R1, R2 for a pair;
R1, R2, R1, R2, ... for several pairs in rna_fq2bam): the wrapper passes their
names so the task keeps that order, whatever their names sort as.

Licence: Parabricks is free of charge but proprietary, not open source. The
container ships NVIDIA's own terms: the NVIDIA AI Product Agreement
(``/usr/local/parabricks/AI_EULA.txt``, which you accept by using it), the
Parabricks licence (``NVIDIA-Parabricks-License-29-09-22.pdf``) and NGC's
container licence (``/NGC-DL-CONTAINER-LICENSE``). They allow use and
deployment on infrastructure you own or lease, but not redistributing the
software.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

PARABRICKS_IMAGE = "nvcr.io/nvidia/clara/clara-parabricks:4.8.0-1"
# STAR 2.7.2a: the STAR whose index format rna_fq2bam reads (upstream's test builds its index with this image).
STAR_272A_IMAGE = "community.wave.seqera.io/library/htslib_samtools_star_gawk:4de2f983041d42e6"

DEFAULT_RESOURCES = flyte.Resources(cpu=12, memory="48Gi", gpu=1)
STAR_INDEX_RESOURCES = flyte.Resources(cpu=4, memory="16Gi")

# The GPUs this task sees, as Parabricks' --num-gpus; Parabricks can't run without one.
GPU_CHECK = r"""
        NUM_GPUS=$(nvidia-smi -L 2>/dev/null | grep -c '^GPU' || true)
        if [ "${NUM_GPUS:-0}" -eq 0 ]; then echo "Parabricks needs an NVIDIA GPU; this task has none" >&2; exit 1; fi
        # pbrun, with its log's last lines on stderr and in the pod's termination message if it fails
        # (pbrun reports errors on stdout).
        run_pbrun() {
            pbrun "$@" > pbrun.out 2>&1
            local STATUS=$?
            cat pbrun.out
            if [ $STATUS -ne 0 ]; then
                tail -n 40 pbrun.out >&2
                # Kubernetes reports this file as the pod's termination message.
                {
                    grep -iE "error|exception|fail|not found|could not|cannot" pbrun.out | tail -n 12
                    tail -c 600 pbrun.out
                } | tail -c 1800 > /dev/termination-log 2> /dev/null
                exit 1
            fi
            rm -f pbrun.out
        }
        # Parabricks resolves symlinks and looks for companion files (indexes) next to the real file, so inputs
        # that need companions are hard-linked, or copied, side by side (upstream stages them by copying).
        stage() { ln "$1" "$2"/ 2> /dev/null || cp "$1" "$2"/; }
""".strip()

# The reads, in the order of their names in {inputs.read_names} (a list input is staged in name order).
ORDERED_READS = r"""
        shopt -s nullglob; STAGED=({inputs.reads}); shopt -u nullglob
        NAMES={inputs.read_names}
        READS=()
        for N in $NAMES; do
            for F in "${STAGED[@]}"; do if [ "$(basename "$F")" = "$N" ]; then READS+=("$F"); fi; done
        done
        if [ ${#READS[@]} -ne ${#STAGED[@]} ]; then
            echo "reads: their names don't match the staged files" >&2
            exit 1
        fi
""".strip()

READS_INPUTS = {"reads": list[File], "read_names": str, "single_end": bool}

fq2bam_cmd = shell.create(
    name="parabricks_fq2bam",
    image=PARABRICKS_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        **READS_INPUTS,
        "fasta": File,
        "index": Dir,
        "intervals": list[File],
        "known_sites": list[File],
        "known_sites_index": list[File],
        "output_fmt": str,
        "prefix": str,
        "args": str,
    },
    defaults={"output_fmt": "bam", "args": ""},
    outputs={"results": Dir},
    script=f"""
        {TASK_LIMITS}
        {GPU_CHECK}
        {ORDERED_READS}
"""
    r"""
        FA=({inputs.fasta})
        shopt -s nullglob
        INTERVALS=({inputs.intervals}); SITES=({inputs.known_sites}); SITE_INDEXES=({inputs.known_sites_index})
        shopt -u nullglob
        P={inputs.prefix}
        FMT={inputs.output_fmt}
        ARGS={inputs.args}
        OUT={outputs.results}
        W=$(mktemp -d); cd "$W"

        # The BWA index's files and the FASTA, copied to the index's prefix, side by side.
        mkdir ref sites
        for F in {inputs.index}/*; do ln -s "$F" ref/; done
        AMB=(ref/*.amb)
        INDEX="$W/${AMB[0]%.amb}"
        cp "${FA[0]}" "$INDEX"

        if [ {inputs.single_end} = true ]; then IN_FQ=(--in-se-fq "${READS[@]}"); else IN_FQ=(--in-fq "${READS[@]}"); fi
        EXTRA=()
        for F in "${SITE_INDEXES[@]}"; do stage "$F" sites; done
        for F in "${SITES[@]}"; do stage "$F" sites; EXTRA+=(--knownSites "sites/$(basename "$F")"); done
        if [ ${#SITES[@]} -gt 0 ]; then EXTRA+=(--out-recal-file "$P.table"); fi
        for F in "${INTERVALS[@]}"; do EXTRA+=(--interval-file "$F"); done

        run_pbrun fq2bam \
            --ref "$INDEX" \
            "${IN_FQ[@]}" \
            --out-bam "$P.$FMT" \
            "${EXTRA[@]}" \
            --num-gpus $NUM_GPUS \
            --bwa-cpu-thread-pool $TASK_CPUS \
            --monitor-usage \
            $ARGS
        shopt -s nullglob
        mv ./*.bam ./*.bai ./*.cram ./*.crai ./*.table ./*_qc_metrics ./*.duplicate-metrics.txt "$OUT"/
    """,
)

fq2bammeth_cmd = shell.create(
    name="parabricks_fq2bammeth",
    image=PARABRICKS_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        **READS_INPUTS,
        "fasta": File,
        "index": Dir,
        "known_sites": list[File],
        "prefix": str,
        "args": str,
    },
    defaults={"args": ""},
    outputs={"results": Dir},
    script=f"""
        {TASK_LIMITS}
        {GPU_CHECK}
        {ORDERED_READS}
"""
    r"""
        FA=({inputs.fasta})
        shopt -s nullglob; SITES=({inputs.known_sites}); shopt -u nullglob
        P={inputs.prefix}
        ARGS={inputs.args}
        OUT={outputs.results}
        W=$(mktemp -d); cd "$W"

        # As upstream: the FASTA linked into the bwa-meth index, as <index>/<fasta name>.
        NAME=$(basename "${FA[0]}")
        mkdir index
        for F in {inputs.index}/*; do stage "$F" index; done
        cp "${FA[0]}" "index/$NAME"

        if [ {inputs.single_end} = true ]; then IN_FQ=(--in-se-fq "${READS[@]}"); else IN_FQ=(--in-fq "${READS[@]}"); fi
        EXTRA=()
        for F in "${SITES[@]}"; do EXTRA+=(--knownSites "$F"); done
        if [ ${#SITES[@]} -gt 0 ]; then EXTRA+=(--out-recal-file "$P.table"); fi

        run_pbrun fq2bam_meth \
            --ref "index/$NAME" \
            "${IN_FQ[@]}" \
            --out-bam "$P.bam" \
            "${EXTRA[@]}" \
            --num-gpus $NUM_GPUS \
            $ARGS
        shopt -s nullglob
        for F in ./*.bam ./*.bai ./*.table ./qc_metrics ./duplicate-metrics.txt; do
            if [ -e "$F" ]; then mv "$F" "$OUT"/; fi
        done
    """,
)

rnafq2bam_cmd = shell.create(
    name="parabricks_rnafq2bam",
    image=PARABRICKS_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        **READS_INPUTS,
        "fasta": File,
        "index": Dir,
        "qc_metrics": bool,
        "mark_duplicates": bool,
        "prefix": str,
        "args": str,
    },
    defaults={"qc_metrics": False, "mark_duplicates": False, "args": ""},
    outputs={"results": Dir},
    script=f"""
        {TASK_LIMITS}
        {GPU_CHECK}
        {ORDERED_READS}
"""
    r"""
        FA=({inputs.fasta})
        P={inputs.prefix}
        ARGS={inputs.args}
        OUT={outputs.results}
        W=$(mktemp -d); cd "$W"

        IN_FQ=()
        if [ {inputs.single_end} = true ]; then
            for F in "${READS[@]}"; do IN_FQ+=(--in-se-fq "$F"); done
        else
            if [ $(( ${#READS[@]} % 2 )) -ne 0 ]; then
                echo "Paired-end input needs an even number of FASTQ files ordered R1, R2 per pair;" \
                    "got ${#READS[@]}." >&2
                exit 1
            fi
            for (( I = 0; I < ${#READS[@]}; I += 2 )); do IN_FQ+=(--in-fq "${READS[$I]}" "${READS[$I + 1]}"); done
        fi
        EXTRA=()
        if [ {inputs.qc_metrics} = true ]; then EXTRA+=(--out-qc-metrics-dir "${P}_qc_metrics"); fi
        if [ {inputs.mark_duplicates} = true ]; then
            EXTRA+=(--out-duplicate-metrics "$P.duplicate-metrics.txt")
        else
            EXTRA+=(--no-markdups)
        fi

        run_pbrun rna_fq2bam \
            --ref "${FA[0]}" \
            "${IN_FQ[@]}" \
            --output-dir . \
            --genome-lib-dir {inputs.index} \
            --out-bam "$P.bam" \
            --logfile "$P.pbrun.log" \
            --out-prefix "$P." \
            --num-gpus $NUM_GPUS \
            "${EXTRA[@]}" \
            $ARGS
        rm -f "$P.pbrun.log"
        shopt -s nullglob
        for F in ./*; do if [ -f "$F" ] || [ "$F" = "./${P}_qc_metrics" ]; then mv "$F" "$OUT"/; fi; done
    """,
)

# Upstream's star/genomegenerate script, run with STAR 2.7.2a.
star_index_cmd = shell.create(
    name="parabricks_star_index",
    image=STAR_272A_IMAGE,
    resources=STAR_INDEX_RESOURCES,
    inputs={"fasta": File, "gtf": File | None, "args": str},
    defaults={"args": ""},
    outputs={"index": Dir},
    script=f"""
        {TASK_LIMITS}
"""
    r"""
        FA=({inputs.fasta})
        shopt -s nullglob; GTF=({inputs.gtf}); shopt -u nullglob
        ARGS={inputs.args}
        W=$(mktemp -d); cd "$W"
        NAME=$(basename "${FA[0]}")
        ln -s "${FA[0]}" "$NAME"
        EXTRA=()
        if [ ${#GTF[@]} -gt 0 ]; then EXTRA+=(--sjdbGTFfile "${GTF[0]}"); fi
        if [[ " $ARGS " != *" --genomeSAindexNbases "* ]]; then
            samtools faidx "$NAME"
            NUM_BASES=$(gawk '{sum = sum + $2}END{if ((log(sum)/log(2))/2 - 1 > 14) {printf "%.0f", 14}
                else {printf "%d", int((log(sum)/log(2))/2 - 1)}}' "$NAME.fai")
            EXTRA+=(--genomeSAindexNbases $NUM_BASES)
        fi
        STAR \
            --runMode genomeGenerate \
            --genomeDir {outputs.index}/ \
            --genomeFastaFiles "$NAME" \
            "${EXTRA[@]}" \
            --runThreadN $TASK_CPUS \
            --limitGenomeGenerateRAM $(( TASK_MEMORY_MB * 1048576 * 9 / 10 )) \
            $ARGS
    """,
)

# The reference, staged into a writable directory: Parabricks writes the FASTA's .fai next to it when missing.
LINK_REFERENCE = r"""
        FA=({inputs.fasta})
        W=$(mktemp -d); cd "$W"
        mkdir ref
        stage "${FA[0]}" ref
        REF="ref/$(basename "${FA[0]}")"
""".strip()

applybqsr_cmd = shell.create(
    name="parabricks_applybqsr",
    image=PARABRICKS_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "input": File,
        "input_index": File | None,
        "bqsr_table": File,
        "intervals": list[File],
        "fasta": File,
        "output_fmt": str,
        "prefix": str,
        "args": str,
    },
    defaults={"output_fmt": "bam", "args": ""},
    outputs={"results": Dir},
    script=f"""
        {TASK_LIMITS}
        {GPU_CHECK}
"""
    r"""
        IN=({inputs.input}); TABLE=({inputs.bqsr_table})
        shopt -s nullglob; IDX=({inputs.input_index}); INTERVALS=({inputs.intervals}); shopt -u nullglob
        P={inputs.prefix}
        FMT={inputs.output_fmt}
        ARGS={inputs.args}
        OUT={outputs.results}
"""
    f"        {LINK_REFERENCE}\n"
    r"""
        mkdir aln
        stage "${IN[0]}" aln
        if [ ${#IDX[@]} -gt 0 ]; then stage "${IDX[0]}" aln; fi
        EXTRA=()
        for F in "${INTERVALS[@]}"; do EXTRA+=(--interval-file "$F"); done

        run_pbrun applybqsr \
            --ref "$REF" \
            --in-bam "aln/$(basename "${IN[0]}")" \
            --in-recal-file "${TABLE[0]}" \
            "${EXTRA[@]}" \
            --out-bam "$P.$FMT" \
            --num-threads $TASK_CPUS \
            --num-gpus $NUM_GPUS \
            $ARGS
        shopt -s nullglob
        mv ./*.bam ./*.bai ./*.cram ./*.crai "$OUT"/
    """,
)

# The alignment and its index (when given), side by side.
LINK_ALIGNMENT = r"""
        mkdir aln
        stage "${IN[0]}" aln
        if [ ${#IDX[@]} -gt 0 ]; then stage "${IDX[0]}" aln; fi
        EXTRA=()
        for F in "${INTERVALS[@]}"; do EXTRA+=(--interval-file "$F"); done
""".strip()

CALLER_INPUTS = {
    "input": File,
    "input_index": File | None,
    "intervals": list[File],
    "fasta": File,
    "prefix": str,
    "args": str,
}

CALLER_PREAMBLE = (
    f"""
        {TASK_LIMITS}
        {GPU_CHECK}
"""
    r"""
        IN=({inputs.input})
        shopt -s nullglob; IDX=({inputs.input_index}); INTERVALS=({inputs.intervals}); shopt -u nullglob
        P={inputs.prefix}
        ARGS={inputs.args}
        OUT={outputs.results}
"""
    f"""
        {LINK_REFERENCE}
        {LINK_ALIGNMENT}
"""
)

haplotypecaller_cmd = shell.create(
    name="parabricks_haplotypecaller",
    image=PARABRICKS_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs=CALLER_INPUTS,
    defaults={"args": ""},
    outputs={"results": Dir},
    script=CALLER_PREAMBLE
    + r"""
        OUTPUT="$P.vcf"
        if [[ " $ARGS " == *" --gvcf "* ]]; then OUTPUT="$P.g.vcf.gz"; fi
        run_pbrun haplotypecaller \
            --ref "$REF" \
            --in-bam "aln/$(basename "${IN[0]}")" \
            --out-variants "$OUTPUT" \
            "${EXTRA[@]}" \
            --num-gpus $NUM_GPUS \
            $ARGS
        shopt -s nullglob
        mv ./*.vcf ./*.g.vcf.gz ./*.tbi "$OUT"/
    """,
)

deepvariant_cmd = shell.create(
    name="parabricks_deepvariant",
    image=PARABRICKS_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs=CALLER_INPUTS,
    defaults={"args": ""},
    outputs={"results": Dir},
    script=CALLER_PREAMBLE
    + r"""
        OUTPUT="$P.vcf.gz"
        if [[ " $ARGS " == *" --gvcf "* ]]; then OUTPUT="$P.g.vcf.gz"; fi
        run_pbrun deepvariant \
            --ref "$REF" \
            --in-bam "aln/$(basename "${IN[0]}")" \
            --out-variants "$OUTPUT" \
            "${EXTRA[@]}" \
            --num-gpus $NUM_GPUS \
            $ARGS
        shopt -s nullglob
        mv ./*.vcf.gz ./*.tbi "$OUT"/
    """,
)

deepsomatic_cmd = shell.create(
    name="parabricks_deepsomatic",
    image=PARABRICKS_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "tumor": File,
        "tumor_index": File | None,
        "normal": File,
        "normal_index": File | None,
        "intervals": list[File],
        "fasta": File,
        "prefix": str,
        "args": str,
    },
    defaults={"args": ""},
    outputs={"results": Dir},
    script=f"""
        {TASK_LIMITS}
        {GPU_CHECK}
"""
    r"""
        TUMOR=({inputs.tumor}); NORMAL=({inputs.normal})
        shopt -s nullglob
        TUMOR_IDX=({inputs.tumor_index}); NORMAL_IDX=({inputs.normal_index}); INTERVALS=({inputs.intervals})
        shopt -u nullglob
        P={inputs.prefix}
        ARGS={inputs.args}
        OUT={outputs.results}
"""
    f"        {LINK_REFERENCE}\n"
    r"""
        mkdir tumor normal
        stage "${TUMOR[0]}" tumor
        stage "${NORMAL[0]}" normal
        if [ ${#TUMOR_IDX[@]} -gt 0 ]; then stage "${TUMOR_IDX[0]}" tumor; fi
        if [ ${#NORMAL_IDX[@]} -gt 0 ]; then stage "${NORMAL_IDX[0]}" normal; fi
        EXTRA=()
        for F in "${INTERVALS[@]}"; do EXTRA+=(--interval-file "$F"); done
        run_pbrun deepsomatic \
            --ref "$REF" \
            --in-tumor-bam "tumor/$(basename "${TUMOR[0]}")" \
            --in-normal-bam "normal/$(basename "${NORMAL[0]}")" \
            --out-variants "$P.vcf.gz" \
            "${EXTRA[@]}" \
            --num-gpus $NUM_GPUS \
            $ARGS
        shopt -s nullglob
        mv ./*.vcf.gz ./*.tbi "$OUT"/
    """,
)

parabricks_env = flyte.TaskEnvironment.from_task(
    "parabricks_gpu",
    fq2bam_cmd.as_task(),
    fq2bammeth_cmd.as_task(),
    rnafq2bam_cmd.as_task(),
    applybqsr_cmd.as_task(),
    haplotypecaller_cmd.as_task(),
    deepvariant_cmd.as_task(),
    deepsomatic_cmd.as_task(),
)
star_index_env = flyte.TaskEnvironment.from_task("parabricks_star_index", star_index_cmd.as_task())
env = flyte.TaskEnvironment(name="parabricks", depends_on=[parabricks_env, star_index_env])


async def pick(results: Dir, name: str, tool: str) -> File:
    found = await results.get_file(name)
    if found is None:
        raise FileNotFoundError(f"parabricks {tool} wrote no {name}")
    return found


async def pick_optional(results: Dir, name: str) -> File | None:
    return await results.get_file(name)


def read_names(reads: list[File]) -> str:
    names = [read.name for read in reads]
    if len(set(names)) != len(names):
        raise ValueError(f"parabricks: read files need distinct names, got {names}")
    return " ".join(names)


@dataclass
class AlignmentResult:
    alignment: File  # <prefix>.bam or <prefix>.cram
    index: File  # .bai or .crai
    bqsr_table: File | None  # <prefix>.table, with known sites
    results: Dir  # everything the tool wrote (QC metrics, duplicate metrics, ...)


async def alignment_result(results: Dir, prefix: str, output_fmt: str, tool: str, table: str) -> AlignmentResult:
    index_ext = "crai" if output_fmt == "cram" else "bai"
    return AlignmentResult(
        alignment=await pick(results, f"{prefix}.{output_fmt}", tool),
        index=await pick(results, f"{prefix}.{output_fmt}.{index_ext}", tool),
        bqsr_table=await pick_optional(results, table),
        results=results,
    )


async def parabricks_fq2bam(
    reads: list[File],
    fasta: File,
    index: Dir,
    single_end: bool = False,
    intervals: list[File] | None = None,
    known_sites: list[File] | None = None,
    known_sites_index: list[File] | None = None,
    output_fmt: str = "bam",
    prefix: str = "sample",
    args: str = "",
) -> AlignmentResult:
    """Align ``reads`` (one file, or R1 and R2) to the BWA ``index`` of ``fasta``; ``output_fmt`` is bam or cram."""
    results = await fq2bam_cmd(
        reads=reads,
        read_names=read_names(reads),
        single_end=single_end,
        fasta=fasta,
        index=index,
        intervals=intervals or [],
        known_sites=known_sites or [],
        known_sites_index=known_sites_index or [],
        output_fmt=output_fmt,
        prefix=prefix,
        args=args,
    )
    return await alignment_result(results, prefix, output_fmt, "fq2bam", f"{prefix}.table")


async def parabricks_fq2bammeth(
    reads: list[File],
    fasta: File,
    index: Dir,
    single_end: bool = False,
    known_sites: list[File] | None = None,
    prefix: str = "sample",
    args: str = "",
) -> AlignmentResult:
    """Align bisulfite ``reads`` (one file, or R1 and R2) to the bwa-meth ``index`` of ``fasta``."""
    results = await fq2bammeth_cmd(
        reads=reads,
        read_names=read_names(reads),
        single_end=single_end,
        fasta=fasta,
        index=index,
        known_sites=known_sites or [],
        prefix=prefix,
        args=args,
    )
    return await alignment_result(results, prefix, "bam", "fq2bammeth", f"{prefix}.table")


@dataclass
class RnaFq2BamResult:
    bam: File | None  # <prefix>.bam
    bai: File | None  # <prefix>.bam.bai
    log_final: File  # <prefix>.Log.final.out
    log_out: File  # <prefix>.Log.out
    log_progress: File  # <prefix>.Log.progress.out
    splice_junctions: File | None  # <prefix>.SJ.out.tab
    junction: File | None  # <prefix>.Chimeric.out.junction, with --out-chim-type Junctions
    duplicate_metrics: File | None  # <prefix>.duplicate-metrics.txt, with mark_duplicates
    results: Dir  # everything rna_fq2bam wrote (QC metrics, counts, transcriptome BAM, ...)


async def parabricks_rnafq2bam(
    reads: list[File],
    fasta: File,
    index: Dir,
    single_end: bool = False,
    qc_metrics: bool = False,
    mark_duplicates: bool = False,
    prefix: str = "sample",
    args: str = "",
) -> RnaFq2BamResult:
    """Align RNA ``reads`` (R1, R2 per pair, or single-end files) to a STAR 2.7.2a ``index`` of ``fasta``."""
    results = await rnafq2bam_cmd(
        reads=reads,
        read_names=read_names(reads),
        single_end=single_end,
        fasta=fasta,
        index=index,
        qc_metrics=qc_metrics,
        mark_duplicates=mark_duplicates,
        prefix=prefix,
        args=args,
    )
    return RnaFq2BamResult(
        bam=await pick_optional(results, f"{prefix}.bam"),
        bai=await pick_optional(results, f"{prefix}.bam.bai"),
        log_final=await pick(results, f"{prefix}.Log.final.out", "rna_fq2bam"),
        log_out=await pick(results, f"{prefix}.Log.out", "rna_fq2bam"),
        log_progress=await pick(results, f"{prefix}.Log.progress.out", "rna_fq2bam"),
        splice_junctions=await pick_optional(results, f"{prefix}.SJ.out.tab"),
        junction=await pick_optional(results, f"{prefix}.Chimeric.out.junction"),
        duplicate_metrics=await pick_optional(results, f"{prefix}.duplicate-metrics.txt"),
        results=results,
    )


async def parabricks_star_index(fasta: File, gtf: File | None = None, args: str = "") -> Dir:
    """A STAR 2.7.2a genome index of ``fasta`` (with ``gtf``'s splice junctions), as ``rna_fq2bam`` needs."""
    return await star_index_cmd(fasta=fasta, gtf=gtf, args=args)


async def parabricks_applybqsr(
    input: File,
    bqsr_table: File,
    fasta: File,
    input_index: File | None = None,
    intervals: list[File] | None = None,
    output_fmt: str = "bam",
    prefix: str = "sample",
    args: str = "",
) -> AlignmentResult:
    """Apply ``bqsr_table`` to a BAM or CRAM; ``output_fmt`` is bam or cram."""
    results = await applybqsr_cmd(
        input=input,
        input_index=input_index,
        bqsr_table=bqsr_table,
        intervals=intervals or [],
        fasta=fasta,
        output_fmt=output_fmt,
        prefix=prefix,
        args=args,
    )
    return await alignment_result(results, prefix, output_fmt, "applybqsr", f"{prefix}.table")


@dataclass
class VariantsResult:
    vcf: File  # <prefix>.vcf (haplotypecaller) or <prefix>.vcf.gz; <prefix>.g.vcf.gz with --gvcf
    gvcf: bool  # whether vcf is a GVCF
    results: Dir


async def call_variants(cmd, tool: str, plain_name: str, prefix: str, args: str, **inputs) -> VariantsResult:
    results = await cmd(prefix=prefix, args=args, **inputs)
    gvcf = "--gvcf" in args.split()
    name = f"{prefix}.g.vcf.gz" if gvcf else plain_name
    return VariantsResult(vcf=await pick(results, name, tool), gvcf=gvcf, results=results)


async def parabricks_haplotypecaller(
    input: File,
    fasta: File,
    input_index: File | None = None,
    intervals: list[File] | None = None,
    prefix: str = "sample",
    args: str = "",
) -> VariantsResult:
    """Call germline variants in a BAM/CRAM with GPU HaplotypeCaller: ``<prefix>.vcf`` (or ``.g.vcf.gz``)."""
    return await call_variants(
        haplotypecaller_cmd,
        "haplotypecaller",
        f"{prefix}.vcf",
        prefix,
        args,
        input=input,
        input_index=input_index,
        intervals=intervals or [],
        fasta=fasta,
    )


async def parabricks_deepvariant(
    input: File,
    fasta: File,
    input_index: File | None = None,
    intervals: list[File] | None = None,
    prefix: str = "sample",
    args: str = "",
) -> VariantsResult:
    """Call germline variants in a BAM/CRAM with GPU DeepVariant: ``<prefix>.vcf.gz`` (or ``.g.vcf.gz``)."""
    return await call_variants(
        deepvariant_cmd,
        "deepvariant",
        f"{prefix}.vcf.gz",
        prefix,
        args,
        input=input,
        input_index=input_index,
        intervals=intervals or [],
        fasta=fasta,
    )


@dataclass
class DeepSomaticResult:
    vcf: File  # <prefix>.vcf.gz
    tbi: File | None  # <prefix>.vcf.gz.tbi
    results: Dir


async def parabricks_deepsomatic(
    tumor: File,
    normal: File,
    fasta: File,
    tumor_index: File | None = None,
    normal_index: File | None = None,
    intervals: list[File] | None = None,
    prefix: str = "sample",
    args: str = "",
) -> DeepSomaticResult:
    """Call somatic variants in a tumour/normal pair with GPU DeepSomatic: ``<prefix>.vcf.gz`` + ``.tbi``."""
    results = await deepsomatic_cmd(
        tumor=tumor,
        tumor_index=tumor_index,
        normal=normal,
        normal_index=normal_index,
        intervals=intervals or [],
        fasta=fasta,
        prefix=prefix,
        args=args,
    )
    return DeepSomaticResult(
        vcf=await pick(results, f"{prefix}.vcf.gz", "deepsomatic"),
        tbi=await pick_optional(results, f"{prefix}.vcf.gz.tbi"),
        results=results,
    )
