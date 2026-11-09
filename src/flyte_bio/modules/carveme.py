"""carveme — CarveMe, genome-scale metabolic model reconstruction from a proteome.

Exposes :func:`carveme_carve`, which carves an SBML model (``<prefix>.xml``)
out of CarveMe's universal model for one protein FASTA (``carve``). With
``--debug`` in ``args``, CarveMe also writes its intermediate tables:
``<prefix>_gene_scores.tsv``, ``<prefix>_protein_scores.tsv``,
``<prefix>_reaction_scores.tsv``, ``<prefix>_milp_problem.lp`` and
``<prefix>_milp_solution.tsv``.

Optional inputs, as upstream: ``universe_file`` (``--universe-file``), a
media database (``mediadb``, ``--mediadb``, used with ``--gapfill`` or
``--init`` in ``args``), soft and hard constraints (``--soft``, ``--hard``)
and a reference model (``--reference``).

The container ships the open-source SCIP solver (Apache-2.0), which CarveMe
picks when neither CPLEX nor Gurobi is installed, so no solver licence is
needed. CarveMe aligns the proteome with DIAMOND; as upstream, it runs on one
CPU.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

CARVEME_IMAGE = "community.wave.seqera.io/library/carveme_scip:8b66c18ea3a4fb1b"
# Upstream's process_single: one CPU.
DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory="6Gi")

carveme_carve_cmd = shell.create(
    name="carveme_carve",
    image=CARVEME_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "fasta": File,
        "universe_file": File | None,
        "mediadb": File | None,
        "soft": File | None,
        "hard": File | None,
        "reference": File | None,
        "prefix": str,
        "args": str,
    },
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        FA=({{inputs.fasta}})
        shopt -s nullglob
        UNIVERSE=({{inputs.universe_file}}); MEDIADB=({{inputs.mediadb}}); SOFT=({{inputs.soft}})
        HARD=({{inputs.hard}}); REFERENCE=({{inputs.reference}})
        shopt -u nullglob
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        export HOME="$W"
        EXTRA=()
        if [ ${{#UNIVERSE[@]}} -gt 0 ]; then EXTRA+=(--universe-file "${{UNIVERSE[0]}}"); fi
        if [ ${{#MEDIADB[@]}} -gt 0 ]; then EXTRA+=(--mediadb "${{MEDIADB[0]}}"); fi
        if [ ${{#SOFT[@]}} -gt 0 ]; then EXTRA+=(--soft "${{SOFT[0]}}"); fi
        if [ ${{#HARD[@]}} -gt 0 ]; then EXTRA+=(--hard "${{HARD[0]}}"); fi
        if [ ${{#REFERENCE[@]}} -gt 0 ]; then EXTRA+=(--reference "${{REFERENCE[0]}}"); fi
        NAME=$(basename "${{FA[0]}}")
        ln -s "${{FA[0]}}" "$NAME"
        carve "$NAME" --output "$P.xml" "${{EXTRA[@]}}" $ARGS || exit 1
        shopt -s nullglob
        mv "$P.xml" "${{P}}"_*.tsv "${{P}}"_*.lp "$OUT"/
    """,
)

env = flyte.TaskEnvironment.from_task("carveme", carveme_carve_cmd.as_task())


@dataclass
class CarvemeResult:
    model: File  # <prefix>.xml, the SBML model
    gene_scores: File | None  # the --debug tables
    milp_problem: File | None
    milp_solution: File | None
    protein_scores: File | None
    reaction_scores: File | None
    results: Dir  # everything CarveMe wrote


async def carveme_carve(
    fasta: File,
    universe_file: File | None = None,
    mediadb: File | None = None,
    soft: File | None = None,
    hard: File | None = None,
    reference: File | None = None,
    prefix: str = "sample",
    args: str = "",
) -> CarvemeResult:
    """Carve a metabolic model for the proteome ``fasta``."""
    results = await carveme_carve_cmd(
        fasta=fasta,
        universe_file=universe_file,
        mediadb=mediadb,
        soft=soft,
        hard=hard,
        reference=reference,
        prefix=prefix,
        args=args,
    )
    model = await results.get_file(f"{prefix}.xml")
    if model is None:
        raise FileNotFoundError(f"carve wrote no {prefix}.xml")
    return CarvemeResult(
        model=model,
        gene_scores=await results.get_file(f"{prefix}_gene_scores.tsv"),
        milp_problem=await results.get_file(f"{prefix}_milp_problem.lp"),
        milp_solution=await results.get_file(f"{prefix}_milp_solution.tsv"),
        protein_scores=await results.get_file(f"{prefix}_protein_scores.tsv"),
        reaction_scores=await results.get_file(f"{prefix}_reaction_scores.tsv"),
        results=results,
    )
