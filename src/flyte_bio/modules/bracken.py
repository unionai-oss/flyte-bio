"""bracken — species (or other rank) abundance re-estimation from a Kraken2 report.

Exposes :func:`bracken`, which runs Bracken on one sample's Kraken2 report
(``--report`` of :func:`flyte_bio.modules.kraken2.kraken2`) against ``db``, a
Dir holding a Kraken2 database with Bracken's ``database<read length>mers.kmer_distrib``
files. Outputs keep upstream's names:

- ``<prefix>.tsv``: Bracken's abundance table (``-o``);
- ``<prefix>.kraken2.report_bracken.txt``: the Kraken2-style report with the
  re-estimated counts (``-w``).

The rank (``-l``), read length (``-r``) and threshold (``-t``) go in
``args``, as upstream's ``ext.args``; Bracken's defaults are ``-l S -t 10
-r 100``.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

BRACKEN_IMAGE = "community.wave.seqera.io/library/bracken:3.1--22a4e66ce04c5e01"
DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory="2Gi")

bracken_cmd = shell.create(
    name="bracken",
    image=BRACKEN_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"kraken_report": File, "db": Dir, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        REPORT=({inputs.kraken_report})
        DB={inputs.db}
        P={inputs.prefix}
        ARGS={inputs.args}
        OUT={outputs.results}
        W=$(mktemp -d); cd "$W"
        # Bracken's estimation step is Python and may exit 255 on errors, which
        # the cluster takes for an interruption; any failure exits 1 instead.
        bracken \
            $ARGS \
            -d "$DB" \
            -i "${REPORT[0]}" \
            -o "$P.tsv" \
            -w "$P.kraken2.report_bracken.txt" \
            || exit 1
        mv "$P.tsv" "$P.kraken2.report_bracken.txt" "$OUT"/
    """,
)

env = flyte.TaskEnvironment.from_task("bracken", bracken_cmd.as_task())


@dataclass
class BrackenResult:
    report: File  # <prefix>.tsv, the abundance table
    kraken_style_report: File  # <prefix>.kraken2.report_bracken.txt
    results: Dir  # both of the above


async def bracken(kraken_report: File, db: Dir, prefix: str = "sample", args: str = "") -> BrackenResult:
    """Re-estimate the abundances in ``kraken_report`` with Bracken against ``db``."""
    results = await bracken_cmd(kraken_report=kraken_report, db=db, prefix=prefix, args=args)

    async def pick(name: str) -> File:
        found = await results.get_file(name)
        if found is None:
            raise FileNotFoundError(f"bracken wrote no {name}")
        return found

    return BrackenResult(
        report=await pick(f"{prefix}.tsv"),
        kraken_style_report=await pick(f"{prefix}.kraken2.report_bracken.txt"),
        results=results,
    )
