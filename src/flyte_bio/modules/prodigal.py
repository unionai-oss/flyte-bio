"""prodigal — prokaryotic gene prediction with Prodigal.

Exposes :func:`prodigal`, which predicts genes in one genome (FASTA,
optionally gzipped) and writes, as upstream, gzipped with pigz:

- ``<prefix>.<output_format>.gz``: the gene annotations (``-f``: ``gff``,
  ``gbk``, ``sco`` or ``gca``);
- ``<prefix>.fna.gz``: the genes' nucleotide sequences;
- ``<prefix>.faa.gz``: their protein translations;
- ``<prefix>_all.txt.gz``: every potential start, with its scores (``-s``).

Prodigal is single-threaded; pigz uses the task's CPUs.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

PRODIGAL_IMAGE = "community.wave.seqera.io/library/prodigal_pigz:27b5cdb2d516f347"
DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory="2Gi")
OUTPUT_FORMATS = ("gff", "gbk", "sco", "gca")

prodigal_cmd = shell.create(
    name="prodigal",
    image=PRODIGAL_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"genome": File, "output_format": str, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        GENOME=({{inputs.genome}})
        FORMAT={{inputs.output_format}}
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        set -o pipefail
        pigz -cdf "${{GENOME[0]}}" | prodigal \
            $ARGS \
            -f "$FORMAT" \
            -d "$P.fna" \
            -o "$P.$FORMAT" \
            -a "$P.faa" \
            -s "${{P}}_all.txt"
        for F in "$P.fna" "$P.$FORMAT" "$P.faa" "${{P}}_all.txt"; do
            pigz -nm -p $TASK_CPUS "$F"
            mv "$F.gz" "$OUT"/
        done
    """,
)

env = flyte.TaskEnvironment.from_task("prodigal", prodigal_cmd.as_task())


@dataclass
class ProdigalResult:
    gene_annotations: File  # <prefix>.<output_format>.gz
    nucleotide_fasta: File  # <prefix>.fna.gz
    amino_acid_fasta: File  # <prefix>.faa.gz
    all_gene_annotations: File  # <prefix>_all.txt.gz
    results: Dir  # all of the above


async def prodigal(genome: File, output_format: str = "gff", prefix: str = "sample", args: str = "") -> ProdigalResult:
    """Predict genes in ``genome``, writing the annotations in ``output_format``."""
    if output_format not in OUTPUT_FORMATS:
        raise ValueError(f"output_format must be one of {OUTPUT_FORMATS}, got {output_format!r}")
    results = await prodigal_cmd(genome=genome, output_format=output_format, prefix=prefix, args=args)

    async def pick(name: str) -> File:
        found = await results.get_file(name)
        if found is None:
            raise FileNotFoundError(f"prodigal wrote no {name}")
        return found

    return ProdigalResult(
        gene_annotations=await pick(f"{prefix}.{output_format}.gz"),
        nucleotide_fasta=await pick(f"{prefix}.fna.gz"),
        amino_acid_fasta=await pick(f"{prefix}.faa.gz"),
        all_gene_annotations=await pick(f"{prefix}_all.txt.gz"),
        results=results,
    )
