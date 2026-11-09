"""coverm — CoverM, read coverage of contigs and genomes.

Exposes:

- :func:`coverm_contig` — coverage per contig (``coverm contig``) from reads
  mapped to a reference FASTA, or from BAMs mapped already.
- :func:`coverm_genome` — coverage per genome (``coverm genome``) from reads
  or BAMs, with the genomes given as one or more FASTA files or as a
  directory of FASTAs.

Both return ``<prefix>.depth.tsv`` and, with ``enable_bam_output``, the BAMs
CoverM mapped (``<reference>.<reads>.bam``, from its BAM cache).

As upstream:

- reads are ``--single`` (``single_end``), ``--interleaved``
  (``interleaved``) or ``--coupled`` pairs (R1, R2, R1, R2, …); BAMs are
  ``--bam-files`` (``bam_input``);
- ``coverm genome`` reads the genomes from a directory
  (``--genome-fasta-directory``) when ``ref_mode`` is ``dir``, or ``auto``
  with a directory; otherwise from files (``--genome-fasta-files``). Any
  other ``ref_mode`` is an error;
- ``--threads`` follows the task's CPUs, and ``TMPDIR`` is the working
  directory; ``args`` come last but for the output file.

Every input is linked into the working directory under its own name, so the
BAMs are named as upstream's. A list of inputs is staged, and so passed to
CoverM, in name order, so name paired files so that R1 sorts before R2.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

COVERM_IMAGE = "community.wave.seqera.io/library/coverm:0.7.0--c57cbb9db4d660fb"
# Upstream: process_medium. The 4 CPUs match upstream's test CPUs.
DEFAULT_RESOURCES = flyte.Resources(cpu=4, memory="16Gi")


# The subcommands' reference options, inserted into the script after its braces are unescaped.
CONTIG_REFERENCE = r"""
        REFERENCE=()
        if [ {inputs.bam_input} = false ]; then REFERENCE=(--reference "$(name "${REF[0]}")"); fi
"""
GENOME_REFERENCE = r"""
        REF_DIR={inputs.reference_dir}
        # A missing optional Dir is staged as a file, not a directory.
        if [ -d "$REF_DIR" ]; then
            REFERENCE=(--genome-fasta-directory "$REF_DIR")
        else
            REFERENCE=(--genome-fasta-files)
            for F in "${REF[@]}"; do REFERENCE+=("$(name "$F")"); done
        fi
"""


def coverm_script(subcommand: str, reference: str) -> str:
    """The script shared by the subcommands; ``reference`` sets ``REFERENCE``, their reference options."""
    return rf"""
        {TASK_LIMITS}
        INPUT=({{inputs.input}})
        shopt -s nullglob; REF=({{inputs.reference}}); shopt -u nullglob
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        TMPDIR=.
        export TMPDIR

        name() {{ basename "$1"; }}
        INPUTS=()
        for F in "${{INPUT[@]}}"; do ln -s "$F" "$(name "$F")"; INPUTS+=("$(name "$F")"); done
        for F in "${{REF[@]}}"; do ln -s "$F" "$(name "$F")"; done

        if [ {{inputs.bam_input}} = true ]; then
            INPUT_TYPE=--bam-files
        elif [ {{inputs.single_end}} = true ]; then
            INPUT_TYPE=--single
        elif [ {{inputs.interleaved}} = true ]; then
            INPUT_TYPE=--interleaved
        else
            INPUT_TYPE=--coupled
        fi
        {reference}
        BAM_OUTPUT=()
        if [ {{inputs.enable_bam_output}} = true ]; then BAM_OUTPUT=(--bam-file-cache-directory _bam_cache/); fi

        coverm {subcommand} \
            --threads $TASK_CPUS \
            $INPUT_TYPE "${{INPUTS[@]}}" \
            "${{REFERENCE[@]}}" \
            "${{BAM_OUTPUT[@]}}" \
            $ARGS \
            --output-file "$P".depth.tsv || exit 1

        mv "$P".depth.tsv "$OUT"/
        shopt -s nullglob
        for F in _bam_cache/*.bam; do mv "$F" "$OUT"/; done
    """


COMMON_INPUTS = {
    "input": list[File],
    "reference": list[File],
    "bam_input": bool,
    "single_end": bool,
    "interleaved": bool,
    "enable_bam_output": bool,
    "prefix": str,
    "args": str,
}

coverm_contig_cmd = shell.create(
    name="coverm_contig",
    image=COVERM_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs=COMMON_INPUTS,
    defaults={"args": ""},
    outputs={"results": Dir},
    script=coverm_script("contig", CONTIG_REFERENCE),
)

coverm_genome_cmd = shell.create(
    name="coverm_genome",
    image=COVERM_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={**COMMON_INPUTS, "reference_dir": Dir | None},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=coverm_script("genome", GENOME_REFERENCE),
)

env = flyte.TaskEnvironment.from_task("coverm", coverm_contig_cmd.as_task(), coverm_genome_cmd.as_task())


@dataclass
class CovermResult:
    coverage: File  # <prefix>.depth.tsv
    bams: list[File]  # with enable_bam_output: CoverM's BAMs, sorted by name
    results: Dir  # everything above


async def read_coverm_results(results: Dir, prefix: str) -> CovermResult:
    files = await results.list_files()
    coverage = [f for f in files if f.name == f"{prefix}.depth.tsv"]
    if not coverage:
        raise FileNotFoundError(f"CoverM wrote no {prefix}.depth.tsv")
    bams = sorted((f for f in files if f.name.endswith(".bam")), key=lambda f: f.name)
    return CovermResult(coverage=coverage[0], bams=bams, results=results)


def as_list(files: File | list[File]) -> list[File]:
    return files if isinstance(files, list) else [files]


async def coverm_contig(
    input: File | list[File],
    reference: File | None = None,
    bam_input: bool = False,
    single_end: bool = False,
    interleaved: bool = False,
    enable_bam_output: bool = False,
    prefix: str = "sample",
    args: str = "",
) -> CovermResult:
    """Coverage per contig of ``reference`` from reads (mapped here) or, with ``bam_input``, BAMs."""
    if not bam_input and reference is None:
        raise ValueError("coverm contig needs a reference for reads")
    results = await coverm_contig_cmd(
        input=as_list(input),
        reference=[reference] if reference is not None and not bam_input else [],
        bam_input=bam_input,
        single_end=single_end,
        interleaved=interleaved,
        enable_bam_output=enable_bam_output,
        prefix=prefix,
        args=args,
    )
    return await read_coverm_results(results, prefix)


async def coverm_genome(
    input: File | list[File],
    reference: File | list[File] | Dir,
    bam_input: bool = False,
    single_end: bool = False,
    interleaved: bool = False,
    ref_mode: str = "auto",
    enable_bam_output: bool = False,
    prefix: str = "sample",
    args: str = "",
) -> CovermResult:
    """Coverage per genome from reads (mapped here) or BAMs; ``reference`` is genome FASTAs or a directory of them.

    ``ref_mode`` is ``dir``, ``file`` or ``auto`` (empty means ``auto``): a directory with ``auto``.
    """
    mode = ref_mode or "auto"
    allowed = ["dir", "file", "auto"]
    if mode not in allowed:
        raise ValueError(f"Invalid ref_mode='{mode}'. Allowed: {', '.join(allowed)}")
    use_dir = mode == "dir" or (mode == "auto" and isinstance(reference, Dir))
    if use_dir != isinstance(reference, Dir):
        raise ValueError(f"ref_mode {mode!r} needs the reference as {'a Dir' if use_dir else 'FASTA files'}")
    results = await coverm_genome_cmd(
        input=as_list(input),
        reference=[] if isinstance(reference, Dir) else as_list(reference),
        reference_dir=reference if isinstance(reference, Dir) else None,
        bam_input=bam_input,
        single_end=single_end,
        interleaved=interleaved,
        enable_bam_output=enable_bam_output,
        prefix=prefix,
        args=args,
    )
    return await read_coverm_results(results, prefix)
