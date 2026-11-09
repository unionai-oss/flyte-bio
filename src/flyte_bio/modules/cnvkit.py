"""cnvkit — CNVkit, copy-number calling from targeted and whole-genome sequencing.

Exposes:

- :func:`cnvkit_target` — tidy a baits BED into a target BED (optionally
  labelling targets with gene names from an annotation): ``<prefix>.bed``.
- :func:`cnvkit_batch` — CNVkit's whole pipeline (``cnvkit.py batch``): bin
  the targets and antitargets, compute coverages, build or apply a reference,
  and call copy number. Its outputs keep CNVkit's names: target/antitarget
  BEDs (``<targets>.target.bed``, ...), coverages and references (``.cnn``),
  copy ratios (``<sample>.cnr``), segments (``<sample>*.cns``) and any plots.

:func:`cnvkit_batch` follows upstream's modes:

- **tumour-normal**: a tumour with a matched normal (the reference is built
  from the normal against ``fasta`` and ``targets``);
- **tumour-only**: a tumour with a prebuilt ``reference``, or, without one,
  against a flat reference (``--normal`` with no files);
- **germline**: a normal alone, called against a flat reference;
- **panel of normals**: several normals and ``panel_of_normals=True``, which
  builds a reference from them all (name it with ``--output-reference``).

As upstream, CRAMs are converted to BAM with samtools first (CNVkit reads
CRAM natively, but slowly), using ``fasta`` and its ``fai``; every CRAM input
is converted, including a panel of normals' (upstream converts only a single
normal). The task's CPUs go to ``--processes``.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

CNVKIT_IMAGE = "community.wave.seqera.io/library/cnvkit:0.9.14--288e98d6210b7304"
# batch's image adds samtools, for the CRAM conversion.
CNVKIT_SAMTOOLS_IMAGE = "community.wave.seqera.io/library/cnvkit_htslib_samtools:10abfe6987f6efe6"

# Upstream labels both process_low.
TARGET_RESOURCES = flyte.Resources(cpu=1, memory="4Gi")
BATCH_RESOURCES = flyte.Resources(cpu=2, memory="12Gi")

cnvkit_target_cmd = shell.create(
    name="cnvkit_target",
    image=CNVKIT_IMAGE,
    resources=TARGET_RESOURCES,
    inputs={"baits": File, "annotation": File | None, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        BAITS=({inputs.baits})
        shopt -s nullglob; ANNOTATION=({inputs.annotation}); shopt -u nullglob
        ARGS={inputs.args}
        W=$(mktemp -d); cd "$W"
        ANNOTATE=(); if [ ${#ANNOTATION[@]} -gt 0 ]; then ANNOTATE=(--annotate "${ANNOTATION[0]}"); fi
        cnvkit.py target "${BAITS[0]}" "${ANNOTATE[@]}" $ARGS --output {outputs.results}/{inputs.prefix}.bed || exit 1
    """,
)

# Inputs are linked into the working directory (CNVkit's output directory)
# under their own names, as upstream stages them, and are left out of the
# collected outputs.
cnvkit_batch_cmd = shell.create(
    name="cnvkit_batch",
    image=CNVKIT_SAMTOOLS_IMAGE,
    resources=BATCH_RESOURCES,
    inputs={
        "tumor": File | None,
        "tumor_index": File | None,
        "normals": list[File],
        "normal_indices": list[File],
        "fasta": File | None,
        "fai": File | None,
        "targets": File | None,
        "reference": File | None,
        "panel_of_normals": bool,
        "args": str,
    },
    defaults={"panel_of_normals": False, "args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        shopt -s nullglob
        T=({{inputs.tumor}}); TI=({{inputs.tumor_index}}); N=({{inputs.normals}}); NI=({{inputs.normal_indices}})
        FA=({{inputs.fasta}}); FAI=({{inputs.fai}}); TG=({{inputs.targets}}); REF=({{inputs.reference}})
        shopt -u nullglob
        PON={{inputs.panel_of_normals}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        INPUT_NAMES=()
        ALL=("${{T[@]}}" "${{TI[@]}}" "${{N[@]}}" "${{NI[@]}}" "${{FA[@]}}" "${{FAI[@]}}" "${{TG[@]}}" "${{REF[@]}}")
        for f in "${{ALL[@]}}"; do
            ln -s "$f" "$(basename "$f")"; INPUT_NAMES+=("$(basename "$f")")
        done
        FASTA_NAME=""; if [ ${{#FA[@]}} -gt 0 ]; then FASTA_NAME=$(basename "${{FA[0]}}"); fi
        FAI_REFERENCE=()
        if [ ${{#FAI[@]}} -gt 0 ]; then FAI_REFERENCE=(--fai-reference "$(basename "${{FAI[0]}}")"); fi
        # as_bam NAME: sets BAM to NAME, or to its BAM conversion (indexed) if NAME is a CRAM
        as_bam() {{
            BAM="$1"
            if [[ "$1" == *.cram ]]; then
                BAM="${{1%.cram}}.bam"
                samtools view -T "$FASTA_NAME" "${{FAI_REFERENCE[@]}}" "$1" -@ "$TASK_CPUS" -o "$BAM" || exit 1
                samtools index "$BAM" || exit 1
            fi
        }}
        TUMOR=()
        if [ ${{#T[@]}} -gt 0 ]; then as_bam "$(basename "${{T[0]}}")"; TUMOR=("$BAM"); fi
        NORMALS=(); for n in "${{N[@]}}"; do as_bam "$(basename "$n")"; NORMALS+=("$BAM"); done
        NORMAL_ARGS=(); FASTA_ARGS=()
        if [ ${{#N[@]}} -gt 0 ]; then
            if [ -n "$FASTA_NAME" ]; then FASTA_ARGS=(--fasta "$FASTA_NAME"); fi
            if [ ${{#T[@]}} -eq 0 ]; then
                # germline: the normal is the sample, against a flat reference
                TUMOR=("${{NORMALS[0]}}"); NORMAL_ARGS=(--normal)
            else
                NORMAL_ARGS=(--normal "${{NORMALS[0]}}")
            fi
            if [ ${{#REF[@]}} -gt 0 ]; then FASTA_ARGS=(); NORMAL_ARGS=(); fi
        fi
        if [ "$PON" = true ] && [ ${{#T[@]}} -eq 0 ]; then NORMAL_ARGS=(--normal "${{NORMALS[@]}}"); TUMOR=(); fi
        if [ ${{#REF[@]}} -eq 0 ] && [ ${{#N[@]}} -eq 0 ] && [ ${{#T[@]}} -gt 0 ]; then NORMAL_ARGS=(--normal); fi
        TARGET_ARGS=()
        if [ ${{#TG[@]}} -gt 0 ] && [ ${{#REF[@]}} -eq 0 ]; then TARGET_ARGS=(--targets "$(basename "${{TG[0]}}")"); fi
        REFERENCE_ARGS=(); if [ ${{#REF[@]}} -gt 0 ]; then REFERENCE_ARGS=(--reference "$(basename "${{REF[0]}}")"); fi
        cnvkit.py batch \
            "${{TUMOR[@]}}" "${{NORMAL_ARGS[@]}}" "${{FASTA_ARGS[@]}}" "${{REFERENCE_ARGS[@]}}" "${{TARGET_ARGS[@]}}" \
            --processes "$TASK_CPUS" $ARGS || exit 1
        shopt -s nullglob
        for f in *.bed *.cnn *.cnr *.cns *.pdf *.png; do
            keep=1; for i in "${{INPUT_NAMES[@]}}"; do if [ "$f" = "$i" ]; then keep=0; fi; done
            if [ $keep = 1 ]; then mv "$f" "$OUT/"; fi
        done
    """,
)

cnvkit_env = flyte.TaskEnvironment.from_task("cnvkit_tools", cnvkit_target_cmd.as_task())
cnvkit_samtools_env = flyte.TaskEnvironment.from_task("cnvkit_samtools", cnvkit_batch_cmd.as_task())
env = flyte.TaskEnvironment(name="cnvkit", depends_on=[cnvkit_env, cnvkit_samtools_env])


async def cnvkit_target(baits: File, prefix: str, annotation: File | None = None, args: str = "") -> File:
    """``<prefix>.bed``: the target regions of ``baits``, labelled from ``annotation`` (GFF/refFlat/BED) if given."""
    results = await cnvkit_target_cmd(baits=baits, annotation=annotation, prefix=prefix, args=args)
    bed = await results.get_file(f"{prefix}.bed")
    if bed is None:
        raise FileNotFoundError(f"cnvkit target wrote no {prefix}.bed")
    return bed


@dataclass
class CnvkitBatchResult:
    bed: list[File]  # target/antitarget BEDs (e.g. <targets>.target.bed, <targets>.antitarget.bed)
    cnn: list[File]  # per-sample (anti)target coverages, and the reference built (--output-reference)
    cnr: list[File]  # <sample>.cnr bin-level copy ratios
    cns: list[File]  # <sample>.cns segments, <sample>.call.cns and <sample>.bintest.cns
    pdf: list[File]  # plots, when args ask (--scatter, --diagram)
    png: list[File]
    results: Dir  # all of the above


async def cnvkit_batch(
    tumor: File | None = None,
    tumor_index: File | None = None,
    normals: list[File] | None = None,
    normal_indices: list[File] | None = None,
    fasta: File | None = None,
    fai: File | None = None,
    targets: File | None = None,
    reference: File | None = None,
    panel_of_normals: bool = False,
    args: str = "",
) -> CnvkitBatchResult:
    """Run ``cnvkit.py batch`` in the mode the inputs select (see the module docstring).

    ``args`` carries CNVkit options such as ``--method wgs`` and
    ``--output-reference reference.cnn``. Indices are optional (CNVkit and
    samtools index what they need); CRAMs need ``fasta``. Inputs are staged
    under their own names, so they must be distinct.
    """
    normals = normals or []
    if tumor is None and not normals:
        raise ValueError("cnvkit batch needs a tumor or normals")
    if panel_of_normals and tumor is not None:
        raise ValueError("a panel of normals is built from normals alone")
    alignments = [*([tumor] if tumor else []), *normals]
    if fasta is None and any(a.path.endswith(".cram") for a in alignments):
        raise ValueError("CRAM input needs the reference fasta")
    results = await cnvkit_batch_cmd(
        tumor=tumor,
        tumor_index=tumor_index,
        normals=normals,
        normal_indices=normal_indices or [],
        fasta=fasta,
        fai=fai,
        targets=targets,
        reference=reference,
        panel_of_normals=panel_of_normals,
        args=args,
    )
    files = sorted([f async for f in results.walk()], key=lambda f: f.name)

    def with_extension(extension: str) -> list[File]:
        return [f for f in files if f.name.endswith(extension)]

    return CnvkitBatchResult(
        bed=with_extension(".bed"),
        cnn=with_extension(".cnn"),
        cnr=with_extension(".cnr"),
        cns=with_extension(".cns"),
        pdf=with_extension(".pdf"),
        png=with_extension(".png"),
        results=results,
    )
