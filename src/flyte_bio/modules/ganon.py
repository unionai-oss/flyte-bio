"""ganon — taxonomic classification databases (interleaved Bloom filters).

Exposes :func:`ganon_buildcustom`, which builds a ganon database from
reference FASTA files (``ganon build-custom``): either the files themselves
(``--input``), or a TSV that lists them with their targets and taxa
(``--input-file``; it names the files by their basenames). Outputs keep
upstream's names:

- the filter, ``<prefix>.hibf`` (or ``<prefix>.ibf`` with
  ``--filter-type ibf``);
- the taxonomy, ``<prefix>.tax`` (unless ``--taxonomy skip``);
- with ``--write-info-file``, ``<prefix>.info.tsv``.

Without ``taxonomy_files`` ganon downloads the NCBI taxonomy itself, so the
task needs internet access. NCBI taxonomy files are passed in the order ganon
reads them, ``nodes``, ``names`` then ``merged``, whatever order they're given
in. ``--threads`` is the task's CPUs; options such as
``--input-target`` go in ``args``.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

GANON_IMAGE = "quay.io/biocontainers/ganon:2.1.0--py310hab1bfa5_1"
# Upstream labels it process_high; real databases need far more memory.
DEFAULT_RESOURCES = flyte.Resources(cpu=4, memory="16Gi")

ganon_buildcustom_cmd = shell.create(
    name="ganon_buildcustom",
    image=GANON_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "fasta": list[File],
        "input_tsv": File | None,
        "taxonomy_files": list[File],
        "genome_size_files": list[File],
        "prefix": str,
        "args": str,
    },
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        FASTA=({{inputs.fasta}})
        shopt -s nullglob
        TSV=({{inputs.input_tsv}}); TAXONOMY=({{inputs.taxonomy_files}}); SIZES=({{inputs.genome_size_files}})
        shopt -u nullglob
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        NAMES=()
        for F in "${{FASTA[@]}}"; do NAMES+=("$(basename "$F")"); ln -s "$F" "$(basename "$F")"; done
        if [ ${{#TSV[@]}} -gt 0 ]; then
            ln -s "${{TSV[0]}}" "$(basename "${{TSV[0]}}")"
            INPUT=(--input-file "$(basename "${{TSV[0]}}")")
        else
            INPUT=(--input "${{NAMES[@]}}")
        fi
        # NCBI taxonomy files go in ganon's order: nodes, names, merged, then any others.
        ORDERED=()
        for PATTERN in '*nodes*' '*names*' '*merged*'; do
            for F in "${{TAXONOMY[@]}}"; do
                case "$(basename "$F")" in $PATTERN) ORDERED+=("$F") ;; esac
            done
        done
        for F in "${{TAXONOMY[@]}}"; do
            case "$(basename "$F")" in *nodes*|*names*|*merged*) ;; *) ORDERED+=("$F") ;; esac
        done
        TAXONOMY_ARGS=()
        if [ ${{#ORDERED[@]}} -gt 0 ]; then TAXONOMY_ARGS=(--taxonomy-files "${{ORDERED[@]}}"); fi
        SIZE_ARGS=()
        if [ ${{#SIZES[@]}} -gt 0 ]; then SIZE_ARGS=(--genome-size-files "${{SIZES[@]}}"); fi
        ganon \
            build-custom \
            --threads $TASK_CPUS \
            "${{INPUT[@]}}" \
            --db-prefix "$P" \
            "${{TAXONOMY_ARGS[@]}}" \
            "${{SIZE_ARGS[@]}}" \
            $ARGS \
            > ganon.log 2>&1 \
            || {{ tail -n 20 ganon.log >&2; exit 1; }}
        cat ganon.log >&2
        for F in "$P".hibf "$P".ibf "$P".tax "$P".info.tsv; do
            if [ -e "$F" ]; then mv "$F" "$OUT"/; fi
        done
    """,
)

env = flyte.TaskEnvironment.from_task("ganon", ganon_buildcustom_cmd.as_task())


@dataclass
class GanonBuildResult:
    filter: File  # <prefix>.hibf, or <prefix>.ibf
    tax: File | None  # <prefix>.tax
    info: File | None  # <prefix>.info.tsv, with --write-info-file
    db: Dir  # everything above


async def ganon_buildcustom(
    fasta: File | list[File],
    input_tsv: File | None = None,
    taxonomy_files: list[File] | None = None,
    genome_size_files: list[File] | None = None,
    prefix: str = "sample",
    args: str = "",
) -> GanonBuildResult:
    """Build a ganon database of ``fasta``, listed by ``input_tsv`` if given."""
    fastas = fasta if isinstance(fasta, list) else [fasta]
    db = await ganon_buildcustom_cmd(
        fasta=fastas,
        input_tsv=input_tsv,
        taxonomy_files=taxonomy_files or [],
        genome_size_files=genome_size_files or [],
        prefix=prefix,
        args=args,
    )
    bloom = await db.get_file(f"{prefix}.hibf") or await db.get_file(f"{prefix}.ibf")
    if bloom is None:
        raise FileNotFoundError(f"ganon build-custom wrote no {prefix}.hibf or {prefix}.ibf")
    return GanonBuildResult(
        filter=bloom,
        tax=await db.get_file(f"{prefix}.tax"),
        info=await db.get_file(f"{prefix}.info.tsv"),
        db=db,
    )
