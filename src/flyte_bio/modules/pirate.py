"""pirate — pangenome analysis of prokaryote annotations with PIRATE.

Exposes :func:`pirate`, which clusters the genes of several genomes' GFF3
annotations (with their sequences, as Prokka writes them) into gene families
across a range of amino-acid identity thresholds. As upstream, the GFFs are
gathered in one directory, a ``.gff3`` extension renamed ``.gff``, and PIRATE
writes into ``<prefix>_results/``: the gene-family tables
(``PIRATE.gene_families.tsv``, ``PIRATE.gene_families.ordered.tsv``,
``PIRATE.unique_alleles.tsv``, …), the pangenome graph
(``pangenome.gfa``), the presence/absence matrix
(``binary_presence_absence.{fasta,nwk}``), the representative sequences and
PIRATE's logs; with ``-a`` in ``args``, also the core-gene alignment
``core_alignment.fasta``. PIRATE uses the task's CPUs.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

PIRATE_IMAGE = "quay.io/biocontainers/pirate:1.0.5--hdfd78af_0"
# Upstream labels PIRATE process_medium.
DEFAULT_RESOURCES = flyte.Resources(cpu=4, memory="8Gi")

pirate_cmd = shell.create(
    name="pirate",
    image=PIRATE_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"gff": list[File], "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        GFF=({{inputs.gff}})
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        for F in "${{GFF[@]}}"; do
            NAME=$(basename "$F")
            # Rename .gff3 to .gff, as upstream.
            case "$NAME" in *.gff3) NAME="${{NAME%3}}" ;; esac
            ln -s "$F" "$NAME"
        done
        PIRATE \
            $ARGS \
            --threads $TASK_CPUS \
            --input ./ \
            --output "${{P}}_results/" || exit 1
        mv "${{P}}_results" "$OUT"/
    """,
)

env = flyte.TaskEnvironment.from_task("pirate", pirate_cmd.as_task())


@dataclass
class PirateResult:
    results: Dir  # <prefix>_results: everything PIRATE wrote
    aln: File | None  # <prefix>_results/core_alignment.fasta, with -a in args


async def pirate(gff: list[File], prefix: str = "sample", args: str = "") -> PirateResult:
    """PIRATE's pangenome of the genomes annotated in ``gff`` (GFF3 with sequences), in ``<prefix>_results``."""
    names = [f.name for f in gff]
    if len(set(names)) != len(names):
        raise ValueError(f"pirate needs distinct GFF names, got {names}")
    out = await pirate_cmd(gff=gff, prefix=prefix, args=args)
    results = Dir.from_existing_remote(f"{out.path.rstrip('/')}/{prefix}_results")
    return PirateResult(results=results, aln=await results.get_file("core_alignment.fasta"))
