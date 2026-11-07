"""chelae — adapter trimming of FASTQ reads.

Exposes :func:`chelae_trim`, which trims one sample's reads (single-end, or
paired-end mates) with ``chelae trim`` and writes, under upstream's names:

- trimmed reads: ``<prefix>.chelae.fastq.gz`` (single-end) or
  ``<prefix>_R1.chelae.fastq.gz`` / ``<prefix>_R2.chelae.fastq.gz``;
- the JSON report ``<prefix>.chelae.json`` and metrics ``<prefix>.chelae.tsv``.

Inputs may be plain, gzip- or BGZF-compressed FASTQ. An optional adapter
FASTA adds known adapter sequences.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

CHELAE_IMAGE = "quay.io/biocontainers/chelae:0.1.0--hfa8f182_0"

# Upstream's process_medium, capped at the 4 CPUs its tests run with.
DEFAULT_RESOURCES = flyte.Resources(cpu=4, memory="8Gi")

# `single_end` sets the outputs, as upstream's meta.single_end does; passing
# two reads with single_end=True fails, as upstream's mismatched-inputs case.
chelae_trim_cmd = shell.create(
    name="chelae_trim",
    image=CHELAE_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "reads_1": File,
        "reads_2": File | None,
        "adapter_fasta": File | None,
        "single_end": bool,
        "prefix": str,
        "args": str,
    },
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        R1=({{inputs.reads_1}})
        shopt -s nullglob; R2=({{inputs.reads_2}}); ADAPTER=({{inputs.adapter_fasta}}); shopt -u nullglob
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); mkdir "$W/in"; cd "$W"
        INPUTS=()
        for r in "${{R1[@]}}" "${{R2[@]}}"; do ln -s "$r" in/; INPUTS+=("in/$(basename "$r")"); done
        if [ {{inputs.single_end}} = true ]; then
            OUTPUTS=("$P.chelae.fastq.gz")
        else
            OUTPUTS=("${{P}}_R1.chelae.fastq.gz" "${{P}}_R2.chelae.fastq.gz")
        fi
        EXTRA=()
        if [ ${{#ADAPTER[@]}} -gt 0 ]; then EXTRA+=(--adapter-fasta "${{ADAPTER[0]}}"); fi
        chelae trim \
            --inputs "${{INPUTS[@]}}" \
            --outputs "${{OUTPUTS[@]}}" \
            --threads $TASK_CPUS \
            --json "$P.chelae.json" \
            --metrics "$P.chelae.tsv" \
            "${{EXTRA[@]}}" \
            $ARGS || exit 1
        mv "${{OUTPUTS[@]}}" "$P.chelae.json" "$P.chelae.tsv" "$OUT"/
    """,
)

env = flyte.TaskEnvironment.from_task("chelae", chelae_trim_cmd.as_task())


@dataclass
class ChelaeTrimResult:
    reads: list[File]  # the trimmed reads: one file, or R1 and R2
    json: File  # <prefix>.chelae.json
    metrics: File  # <prefix>.chelae.tsv
    results: Dir  # everything the task wrote


async def chelae_trim(
    reads_1: File,
    reads_2: File | None,
    prefix: str,
    adapter_fasta: File | None = None,
    args: str = "",
) -> ChelaeTrimResult:
    """Trim adapters from ``reads_1`` (and mate ``reads_2``) with ``chelae trim``."""
    single_end = reads_2 is None
    results = await chelae_trim_cmd(
        reads_1=reads_1,
        reads_2=reads_2,
        adapter_fasta=adapter_fasta,
        single_end=single_end,
        prefix=prefix,
        args=args,
    )
    names = [f"{prefix}.chelae.fastq.gz"] if single_end else [f"{prefix}_R{i}.chelae.fastq.gz" for i in (1, 2)]

    async def pick(name: str) -> File:
        file = await results.get_file(name)
        if file is None:
            raise FileNotFoundError(f"chelae wrote no {name}")
        return file

    reads = [await pick(name) for name in names]
    return ChelaeTrimResult(
        reads=reads,
        json=await pick(f"{prefix}.chelae.json"),
        metrics=await pick(f"{prefix}.chelae.tsv"),
        results=results,
    )
