"""art — ART read simulation (``art_illumina``).

Exposes :func:`art_illumina`, which simulates Illumina reads from a reference
FASTA for a sequencing system's error profile (``HS25``, ``MSv3``, …), fold
coverage and read length. Pass ART's other options in ``args``: ``-p -m
<mean fragment> -s <sd>`` for paired-end reads, ``-rs <seed>`` for a
reproducible run, ``-sam`` for SAM alignments, ``-na`` to skip the ``.aln``
files.

Outputs keep upstream's names, in one Dir:

- reads: ``<prefix>.fq.gz`` (single-end) or ``<prefix>1.fq.gz`` /
  ``<prefix>2.fq.gz`` (paired-end), gzipped with ``--no-name`` (plus
  ``args2``);
- alignments, when written: ``<prefix>.aln`` or ``<prefix>1.aln`` /
  ``<prefix>2.aln``, and ``<prefix>.sam``.

As upstream, the FASTA is staged under its own name, since ART records its
command line in the ``.aln`` header.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

ART_IMAGE = "community.wave.seqera.io/library/art:2016.06.05--5eb6337227f62cf5"
DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory="4Gi")  # upstream: process_single

art_illumina_cmd = shell.create(
    name="art_illumina",
    image=ART_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "fasta": File,
        "sequencing_system": str,
        "fold_coverage": str,
        "read_length": str,
        "prefix": str,
        "args": str,
        "args2": str,
    },
    defaults={"args": "", "args2": ""},
    outputs={"results": Dir},
    script=r"""
        FA=({inputs.fasta})
        P={inputs.prefix}
        ARGS={inputs.args}
        ARGS2={inputs.args2}
        OUT={outputs.results}
        W=$(mktemp -d); cd "$W"
        NAME=$(basename "${FA[0]}")
        ln -s "${FA[0]}" "$NAME"

        art_illumina \
            -ss {inputs.sequencing_system} \
            -i "$NAME" \
            -l {inputs.read_length} \
            -f {inputs.fold_coverage} \
            -o "$P" \
            $ARGS || exit 1

        gzip \
            --no-name \
            $ARGS2 \
            "$P"*.fq || exit 1
        shopt -s nullglob
        mv ./*.fq.gz ./*.aln ./*.sam "$OUT"/
    """,
)

env = flyte.TaskEnvironment.from_task("art", art_illumina_cmd.as_task())


@dataclass
class ArtIlluminaResult:
    fastq: list[File]  # <prefix>.fq.gz, or <prefix>1.fq.gz and <prefix>2.fq.gz
    aln: list[File]  # <prefix>.aln, or <prefix>1.aln and <prefix>2.aln (unless -na)
    sam: list[File]  # <prefix>.sam (with -sam)
    results: Dir  # everything above


async def art_illumina(
    fasta: File,
    sequencing_system: str,
    fold_coverage: float | int | str,
    read_length: int | str,
    prefix: str = "sample",
    args: str = "",
    args2: str = "",
) -> ArtIlluminaResult:
    """Simulate reads from ``fasta``; ``args`` takes ART's options, ``args2`` gzip's."""
    results = await art_illumina_cmd(
        fasta=fasta,
        sequencing_system=sequencing_system,
        fold_coverage=str(fold_coverage),
        read_length=str(read_length),
        prefix=prefix,
        args=args,
        args2=args2,
    )
    files = sorted(await results.list_files(), key=lambda f: f.name)

    def ending(suffix: str) -> list[File]:
        return [f for f in files if f.name.startswith(prefix) and f.name.endswith(suffix)]

    fastq = ending(".fq.gz")
    if not fastq:
        raise FileNotFoundError(f"art_illumina wrote no {prefix}*.fq.gz")
    return ArtIlluminaResult(fastq=fastq, aln=ending(".aln"), sam=ending(".sam"), results=results)
