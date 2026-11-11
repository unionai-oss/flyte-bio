"""plasmidid — maps contigs onto reference sequences and draws them with Circos.

Exposes :func:`plasmidid`, which runs plasmidID on a sample's contigs (or
scaffolds) against a FASTA of references (upstream's viral pipeline uses
the virus genome) and returns its report directory ``<prefix>/``:

- ``<prefix>_final_results.html`` and ``<prefix>_final_results.tab``;
- ``images/`` (Circos images and configs), ``logs/``, ``data/`` (alignments,
  annotations, coverage), ``database/`` (the filtered references, with
  Prokka's annotation), ``fasta_files/`` (the matched references) and
  ``kmer/`` (mash screening of the references).

``args`` carries plasmidID's options (upstream's test uses ``-k 0.8``).
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

PLASMIDID_IMAGE = "quay.io/biocontainers/plasmidid:1.6.5--hdfd78af_0"
# Upstream: process_medium.
DEFAULT_RESOURCES = flyte.Resources(cpu=2, memory="8Gi")

plasmidid_cmd = shell.create(
    name="plasmidid",
    image=PLASMIDID_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"scaffold": File, "fasta": File, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"report": Dir},
    script=r"""
        SCAFFOLD=({inputs.scaffold}); FA=({inputs.fasta})
        P={inputs.prefix}
        ARGS={inputs.args}
        OUT={outputs.report}
        W=$(mktemp -d); cd "$W"
        REF=$(basename "${FA[0]}"); CONTIGS=$(basename "${SCAFFOLD[0]}")
        ln -s "${FA[0]}" "$REF"
        ln -s "${SCAFFOLD[0]}" "$CONTIGS"
        plasmidID \
            -d "$REF" \
            -s "$P" \
            -c "$CONTIGS" \
            $ARGS \
            -o . || exit 1
        mv NO_GROUP/"$P"/* "$OUT"/
    """,
)

env = flyte.TaskEnvironment.from_task("plasmidid", plasmidid_cmd.as_task())


@dataclass
class PlasmididResult:
    html: File  # <prefix>_final_results.html
    tab: File  # <prefix>_final_results.tab
    report: Dir  # the whole <prefix>/ directory: images/, logs/, data/, database/, fasta_files/, kmer/


async def plasmidid(scaffold: File, fasta: File, prefix: str = "sample", args: str = "") -> PlasmididResult:
    """plasmidID's report on ``scaffold``'s contigs against the references in ``fasta``."""
    report = await plasmidid_cmd(scaffold=scaffold, fasta=fasta, prefix=prefix, args=args)
    picked = []
    for name in (f"{prefix}_final_results.html", f"{prefix}_final_results.tab"):
        found = await report.get_file(name)
        if found is None:
            raise FileNotFoundError(f"plasmidID wrote no {name}")
        picked.append(found)
    return PlasmididResult(html=picked[0], tab=picked[1], report=report)
