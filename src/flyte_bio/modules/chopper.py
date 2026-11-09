"""chopper — quality and length filtering of long reads.

Exposes :func:`chopper`, which filters the reads in ``fastq`` (gzipped) with
chopper, optionally dropping reads that map to a contaminant reference,
``contaminants`` (``--contam``, a FASTA, may be gzipped). As upstream, the
reads are decompressed with pigz, filtered and recompressed with pigz, all on
the task's CPUs, to ``<prefix>.fastq.gz`` (the prefix defaults to
``<sample>.chopper``); an input of that same name is an error.

chopper's options (``-q 10 -l 500``) go in ``args``; ``decompress_args`` and
``compress_args`` pass options to the two pigz steps (upstream's ``args`` and
``args3``). chopper writes reads in the order its threads finish them, so the
output's read order isn't reproducible.
"""

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

# Upstream's Docker image: chopper 0.14.1 and pigz 2.8.
CHOPPER_IMAGE = "community.wave.seqera.io/library/chopper_pigz:4ae8a43e228724c3"
# Upstream's process_medium, sized down.
DEFAULT_RESOURCES = flyte.Resources(cpu=4, memory="8Gi")

chopper_cmd = shell.create(
    name="chopper",
    image=CHOPPER_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "fastq": File,
        "contaminants": File | None,
        "prefix": str,
        "args": str,
        "decompress_args": str,
        "compress_args": str,
    },
    defaults={"args": "", "decompress_args": "", "compress_args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        set -o pipefail
        FQ=({{inputs.fastq}})
        shopt -s nullglob; CONTAM=({{inputs.contaminants}}); shopt -u nullglob
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        UNZIP_ARGS={{inputs.decompress_args}}
        ZIP_ARGS={{inputs.compress_args}}
        OUT={{outputs.results}}
        CONTAM_OPT=()
        if [ ${{#CONTAM[@]}} -gt 0 ]; then CONTAM_OPT=(--contam "${{CONTAM[0]}}"); fi
        pigz \
            --decompress \
            --processes $TASK_CPUS \
            --stdout \
            $UNZIP_ARGS \
            "${{FQ[0]}}" | \
        chopper \
            --threads $TASK_CPUS \
            "${{CONTAM_OPT[@]}}" \
            $ARGS | \
        pigz \
            --processes $TASK_CPUS \
            --stdout \
            $ZIP_ARGS > "$OUT/$P.fastq.gz" || exit 1
    """,
)

env = flyte.TaskEnvironment.from_task("chopper", chopper_cmd.as_task())


async def chopper(
    fastq: File,
    contaminants: File | None = None,
    prefix: str = "sample.chopper",
    args: str = "",
    decompress_args: str = "",
    compress_args: str = "",
) -> File:
    """``fastq``'s reads filtered by chopper (and against ``contaminants``): ``<prefix>.fastq.gz``."""
    filename = f"{prefix}.fastq.gz"
    if fastq.path.rstrip("/").rsplit("/", 1)[-1] == filename:
        raise ValueError("Input and output names are the same, set prefix in module configuration to disambiguate!")
    results = await chopper_cmd(
        fastq=fastq,
        contaminants=contaminants,
        prefix=prefix,
        args=args,
        decompress_args=decompress_args,
        compress_args=compress_args,
    )
    found = await results.get_file(filename)
    if found is None:
        raise FileNotFoundError(f"chopper wrote no {filename}")
    return found
