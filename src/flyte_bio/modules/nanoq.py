"""nanoq — filtering and summary statistics of nanopore reads.

Exposes :func:`nanoq`, which filters one sample's reads (by length or
quality, with ``-l``/``-q`` etc. in ``args``) and writes them in
``output_format``: ``fastq``, ``fasta``, or either compressed (``.gz``,
``.bz2``, ``.lzma``), as ``<prefix>.<output_format>``, with the read summary
``<prefix>.stats`` (``-r``). As upstream, the default prefix is
``<sample>_filtered``.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

NANOQ_IMAGE = "quay.io/biocontainers/nanoq:0.10.0--h031d066_2"
# Upstream labels it process_low; nanoq is single-threaded.
DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory="2Gi")

OUTPUT_FORMATS = tuple(f"{kind}{suffix}" for kind in ("fastq", "fasta") for suffix in ("", ".gz", ".bz2", ".lzma"))

nanoq_cmd = shell.create(
    name="nanoq",
    image=NANOQ_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"reads": File, "output_format": str, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        READS=({inputs.reads})
        FORMAT={inputs.output_format}
        P={inputs.prefix}
        ARGS={inputs.args}
        OUT={outputs.results}
        W=$(mktemp -d); cd "$W"
        nanoq -i "${READS[0]}" \
            $ARGS \
            -r "$P.stats" \
            -o "$P.$FORMAT" \
            || exit 1
        mv "$P.stats" "$P.$FORMAT" "$OUT"/
    """,
)

env = flyte.TaskEnvironment.from_task("nanoq", nanoq_cmd.as_task())


@dataclass
class NanoqResult:
    reads: File  # <prefix>.<output_format>
    stats: File  # <prefix>.stats
    results: Dir  # both of the above


async def nanoq(
    reads: File, output_format: str = "fastq.gz", prefix: str = "sample_filtered", args: str = ""
) -> NanoqResult:
    """Filter one sample's nanopore ``reads``, written as ``<prefix>.<output_format>``."""
    if output_format not in OUTPUT_FORMATS:
        raise ValueError(f"output_format must be one of {', '.join(OUTPUT_FORMATS)}, got {output_format!r}")
    results = await nanoq_cmd(reads=reads, output_format=output_format, prefix=prefix, args=args)
    filtered = await results.get_file(f"{prefix}.{output_format}")
    stats = await results.get_file(f"{prefix}.stats")
    if filtered is None or stats is None:
        raise FileNotFoundError(f"nanoq wrote no {prefix}.{output_format} or {prefix}.stats")
    return NanoqResult(reads=filtered, stats=stats, results=results)
