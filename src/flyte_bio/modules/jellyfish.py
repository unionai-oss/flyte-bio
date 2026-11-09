"""jellyfish — fast k-mer counting.

Exposes :func:`jellyfish_count`, which counts the k-mers of a FASTA/FASTQ
(plain or gzipped; a gzipped input is decompressed first, as upstream) into
a jellyfish hash, ``<prefix>.jf``. ``kmer_length`` is ``-m`` and ``size``
the initial hash size ``-s`` (e.g. ``100M``); ``-t`` follows the task's CPUs.
"""

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

JELLYFISH_IMAGE = "quay.io/biocontainers/kmer-jellyfish:2.3.1--py310h184ae93_5"
DEFAULT_RESOURCES = flyte.Resources(cpu=4, memory="16Gi")

jellyfish_count_cmd = shell.create(
    name="jellyfish_count",
    image=JELLYFISH_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"fasta": File, "kmer_length": int, "size": str, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        FA=({{inputs.fasta}})
        ARGS={{inputs.args}}
        W=$(mktemp -d); cd "$W"
        NAME=$(basename "${{FA[0]}}")
        if [[ "$NAME" == *.gz ]]; then
            NAME=${{NAME%.gz}}
            gzip -c -d "${{FA[0]}}" > "$NAME"
        else
            ln -s "${{FA[0]}}" "$NAME"
        fi
        jellyfish count $ARGS -m {{inputs.kmer_length}} -s {{inputs.size}} -t $TASK_CPUS \
            -o {{outputs.results}}/{{inputs.prefix}}.jf "$NAME" || exit 1
    """,
)

env = flyte.TaskEnvironment.from_task("jellyfish", jellyfish_count_cmd.as_task())


async def jellyfish_count(fasta: File, kmer_length: int, size: str, prefix: str = "sample", args: str = "") -> File:
    """The k-mer counts of ``fasta``: ``<prefix>.jf``."""
    results = await jellyfish_count_cmd(fasta=fasta, kmer_length=kmer_length, size=size, prefix=prefix, args=args)
    jf = await results.get_file(f"{prefix}.jf")
    if jf is None:
        raise FileNotFoundError(f"jellyfish count wrote no {prefix}.jf")
    return jf
