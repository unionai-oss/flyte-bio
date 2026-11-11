"""abacas — orders and orients contigs against a reference genome.

Exposes :func:`abacas`, which runs ABACAS (``abacas.pl``) on a sample's
scaffolds against a reference FASTA and returns everything it writes under
``<prefix>.*`` (upstream's default prefix is ``<id>.abacas``):

- ``<prefix>.fasta`` (the ordered pseudo-molecule), ``<prefix>.MULTIFASTA.fa``,
  ``<prefix>.tab``, ``<prefix>.crunch``, ``<prefix>.gaps``,
  ``<prefix>.gaps.tab`` and ``<prefix>.bin`` (sorted, as upstream, so it is
  stable);
- nucmer's ``<prefix>.nucmer.delta``, ``<prefix>.nucmer.filtered.delta`` and
  ``<prefix>.nucmer.tiling``, and ``<prefix>.unused.contigs.out``.

``args`` carries ABACAS's options (upstream's viral pipeline uses
``-m -p nucmer``).
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

ABACAS_IMAGE = "community.wave.seqera.io/library/abacas:1.3.1--0d1658c2eafe7036"
# Upstream: process_medium; ABACAS and nucmer run single-threaded here.
DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory="4Gi")

abacas_cmd = shell.create(
    name="abacas",
    image=ABACAS_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"scaffold": File, "fasta": File, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        SCAFFOLD=({inputs.scaffold}); FA=({inputs.fasta})
        P={inputs.prefix}
        ARGS={inputs.args}
        OUT={outputs.results}
        W=$(mktemp -d); cd "$W"
        REF=$(basename "${FA[0]}"); QUERY=$(basename "${SCAFFOLD[0]}")
        ln -s "${FA[0]}" "$REF"
        ln -s "${SCAFFOLD[0]}" "$QUERY"
        abacas.pl \
            -r "$REF" \
            -q "$QUERY" \
            $ARGS \
            -o "$P" || exit 1
        sort "$P.bin" > "$P.bin.tmp" && mv "$P.bin.tmp" "$P.bin"
        mv nucmer.delta "$P.nucmer.delta"
        mv nucmer.filtered.delta "$P.nucmer.filtered.delta"
        mv nucmer.tiling "$P.nucmer.tiling"
        mv unused_contigs.out "$P.unused.contigs.out"
        mv "$P".* "$OUT"/
    """,
)

env = flyte.TaskEnvironment.from_task("abacas", abacas_cmd.as_task())


@dataclass
class AbacasResult:
    fasta: File  # <prefix>.fasta, the ordered pseudo-molecule
    multifasta: File  # <prefix>.MULTIFASTA.fa
    results: Dir  # every <prefix>.* file


async def abacas(scaffold: File, fasta: File, prefix: str = "sample.abacas", args: str = "") -> AbacasResult:
    """``scaffold``'s contigs ordered and oriented against the reference ``fasta``."""
    results = await abacas_cmd(scaffold=scaffold, fasta=fasta, prefix=prefix, args=args)
    picked = []
    for name in (f"{prefix}.fasta", f"{prefix}.MULTIFASTA.fa"):
        found = await results.get_file(name)
        if found is None:
            raise FileNotFoundError(f"abacas wrote no {name}")
        picked.append(found)
    return AbacasResult(fasta=picked[0], multifasta=picked[1], results=results)
