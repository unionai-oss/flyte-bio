"""stitch — STITCH imputation of low-coverage sequencing without a reference panel.

Exposes :func:`stitch`, which imputes the samples of BAM or CRAM files (named
in a reads list, ``--bamlist``/``--cramlist``) at the positions of a
``posfile``, over a chromosome or a region of it, with ``K`` ancestral
haplotypes and ``ngen`` generations. It can also run in two stages: with
``--generateInputOnly TRUE`` in ``args`` it only converts the reads into
STITCH's input files; a second run given those ``input/`` and ``RData/``
directories (and no reads list) imputes from them (``--regenerateInput FALSE``).

Outputs keep upstream's names, under one results directory:

- ``input/`` (the per-sample input RData files) and ``RData/``;
- ``plots/``, unless only generating input;
- ``<prefix>.vcf.gz``, or ``<prefix>.bgen`` with ``--output_format bgen`` in
  ``args``.

As upstream, a ``seed`` runs STITCH with ``set.seed(<seed>)`` first, which
makes its output reproducible; ``--nCores`` follows the task's CPUs. Unlike
upstream, ``args`` holds STITCH's options (upstream's ``ext.args2``):
upstream's ``ext.args`` only extends the command that copies ``rdata``.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

STITCH_IMAGE = "quay.io/biocontainers/r-stitch:1.7.3--r44h64f727c_0"
# Upstream labels STITCH process_medium; its tests ran on 4 CPUs.
DEFAULT_RESOURCES = flyte.Resources(cpu=4, memory="8Gi")

# As upstream, a first stage's input files are staged as input/ and its RData
# copied into RData/ (copies: STITCH writes into both).
stitch_cmd = shell.create(
    name="stitch",
    image=STITCH_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "reads": list[File],
        "indexes": list[File],
        "reads_list": File | None,
        "reads_type": str,
        "sample_names": File | None,
        "posfile": File,
        "input_dir": Dir | None,
        "genetic_map": File | None,
        "rdata_dir": Dir | None,
        "chromosome": str,
        "start": str,
        "end": str,
        "buffer": str,
        "k": str,
        "ngen": str,
        "fasta": File | None,
        "fasta_index": File | None,
        "seed": str,
        "prefix": str,
        "suffix": str,
        "args": str,
    },
    defaults={"start": "", "end": "", "buffer": "", "seed": "", "args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        POS=({{inputs.posfile}})
        shopt -s nullglob
        READS=({{inputs.reads}}); IDX=({{inputs.indexes}}); LIST=({{inputs.reads_list}})
        NAMES_FILE=({{inputs.sample_names}}); M=({{inputs.genetic_map}})
        FA=({{inputs.fasta}}); FAI=({{inputs.fasta_index}})
        shopt -u nullglob
        INPUT_DIR={{inputs.input_dir}}; RDATA_DIR={{inputs.rdata_dir}}
        READS_TYPE={{inputs.reads_type}}
        CHR={{inputs.chromosome}}; START={{inputs.start}}; END={{inputs.end}}; BUFFER={{inputs.buffer}}
        K={{inputs.k}}; NGEN={{inputs.ngen}}; SEED={{inputs.seed}}
        P={{inputs.prefix}}; SUFFIX={{inputs.suffix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        link() {{ ln -s "$1" "$(basename "$1")"; basename "$1"; }}
        for F in "${{READS[@]}}" "${{IDX[@]}}"; do link "$F" > /dev/null; done
        if [ -n "$INPUT_DIR" ] && [ -d "$INPUT_DIR" ]; then mkdir input; cp -rL "$INPUT_DIR"/. input/; fi
        if [ -n "$RDATA_DIR" ] && [ -d "$RDATA_DIR" ]; then mkdir RData; cp -rL "$RDATA_DIR"/. RData/; fi
        EXTRA=()
        if [ ${{#LIST[@]}} -gt 0 ]; then
            if [ "$READS_TYPE" = cram ]; then EXTRA+=(--cramlist "$(link "${{LIST[0]}}")"); fi
            if [ "$READS_TYPE" = bam ]; then EXTRA+=(--bamlist "$(link "${{LIST[0]}}")"); fi
        fi
        if [ ${{#FA[@]}} -gt 0 ]; then
            FA_NAME=$(link "${{FA[0]}}"); EXTRA+=(--reference "$FA_NAME")
            if [ ${{#FAI[@]}} -gt 0 ]; then ln -s "${{FAI[0]}}" "$FA_NAME.fai"; fi
        fi
        REGION_NAME=$CHR
        if [ -n "$START" ] && [ -n "$END" ]; then REGION_NAME=$CHR.$START.$END; fi
        if [ -e input ] && [ -e RData ] && [ ${{#LIST[@]}} -eq 0 ]; then
            EXTRA+=(--regenerateInput FALSE --originalRegionName "$REGION_NAME")
        fi
        if [ ${{#NAMES_FILE[@]}} -gt 0 ]; then EXTRA+=(--sampleNames_file "$(link "${{NAMES_FILE[0]}}")"); fi
        if [ ${{#M[@]}} -gt 0 ]; then EXTRA+=(--genetic_map_file="$(link "${{M[0]}}")"); fi
        if [ -n "$START" ] && [ -n "$END" ]; then
            EXTRA+=(--regionStart "$START" --regionEnd "$END")
            if [ -n "$BUFFER" ]; then EXTRA+=(--buffer "$BUFFER"); fi
        fi
        STITCH=(STITCH.R)
        if [ -n "$SEED" ]; then
            tail -n +2 "$(which STITCH.R)" | cat <(echo "set.seed($SEED)") - > stitch_seeded.R
            STITCH=(Rscript stitch_seeded.R)
        fi
        "${{STITCH[@]}}" \
            --chr "$CHR" \
            --posfile "$(link "${{POS[0]}}")" \
            --outputdir . \
            --nCores $TASK_CPUS \
            --K "$K" \
            --nGen "$NGEN" \
            "${{EXTRA[@]}}" \
            --output_filename "$P.$SUFFIX" \
            $ARGS || exit 1
        for D in input RData plots; do if [ -d "$D" ]; then mv "$D" "$OUT"/; fi; done
        if [ -f "$P.$SUFFIX" ]; then mv "$P.$SUFFIX" "$OUT"/; fi
        true
    """,
)

env = flyte.TaskEnvironment.from_task("stitch", stitch_cmd.as_task())


@dataclass
class StitchResult:
    vcf: File | None  # <prefix>.vcf.gz, unless only generating input or writing BGEN
    bgen: File | None  # <prefix>.bgen, with --output_format bgen
    input: list[File]  # input/*: the per-sample input files
    rdata: list[File]  # RData/*
    plots: list[File]  # plots/*
    input_dir: Dir | None  # input/, for a second stage
    rdata_dir: Dir | None  # RData/, for a second stage
    results: Dir  # everything above, under input/, RData/ and plots/


async def stitch(
    posfile: File,
    chromosome: str,
    k: int,
    ngen: int,
    prefix: str,
    reads: list[File] | None = None,
    indexes: list[File] | None = None,
    reads_list: File | None = None,
    sample_names: File | None = None,
    input_dir: Dir | None = None,
    rdata_dir: Dir | None = None,
    genetic_map: File | None = None,
    start: int | None = None,
    end: int | None = None,
    buffer: int | None = None,
    fasta: File | None = None,
    fasta_index: File | None = None,
    seed: int | None = None,
    args: str = "",
) -> StitchResult:
    """Impute the samples in ``reads`` at ``posfile``'s positions of ``chromosome`` (from ``start`` to ``end``).

    ``reads_list`` names the reads files (``--bamlist`` or ``--cramlist``);
    ``sample_names`` renames the samples. For a second stage, pass a first
    run's ``input/`` and ``RData/`` directories as ``input_dir`` and
    ``rdata_dir``, and no reads.
    """
    reads, indexes = list(reads or []), list(indexes or [])
    extensions = sorted({f.name.rsplit(".", 1)[-1] for f in reads})
    if len(extensions) > 1:
        raise ValueError(
            f"STITCH process: Mixed input read file types detected: {extensions}. "
            "Please provide either all BAM or all CRAM files."
        )
    bgen = "--output_format bgen" in args
    input_only = "--generateInputOnly TRUE" in args
    suffix = "bgen" if bgen else "vcf.gz"
    results = await stitch_cmd(
        reads=reads,
        indexes=indexes,
        reads_list=reads_list,
        reads_type=extensions[0] if extensions else "",
        sample_names=sample_names,
        posfile=posfile,
        input_dir=input_dir,
        genetic_map=genetic_map,
        rdata_dir=rdata_dir,
        chromosome=chromosome,
        start="" if start is None else str(start),
        end="" if end is None else str(end),
        buffer="" if buffer is None else str(buffer),
        k=str(k),
        ngen=str(ngen),
        fasta=fasta,
        fasta_index=fasta_index,
        seed="" if seed is None else str(seed),
        prefix=prefix,
        suffix=suffix,
        args=args,
    )
    files = [f async for f in results.walk()]

    def under(directory: str) -> list[File]:
        return sorted((f for f in files if f"/{directory}/" in f.path), key=lambda f: f.name)

    def subdir(directory: str) -> Dir | None:
        return Dir.from_existing_remote(results.path.rstrip("/") + f"/{directory}") if under(directory) else None

    output = await results.get_file(f"{prefix}.{suffix}")
    if output is None and not input_only:
        raise FileNotFoundError(f"STITCH wrote no {prefix}.{suffix}")
    return StitchResult(
        vcf=None if bgen or input_only else output,
        bgen=output if bgen else None,
        input=under("input"),
        rdata=under("RData"),
        plots=under("plots"),
        input_dir=subdir("input"),
        rdata_dir=subdir("RData"),
        results=results,
    )
