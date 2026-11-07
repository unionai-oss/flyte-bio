"""krona — interactive taxonomy charts with KronaTools.

Exposes :func:`krona_ktimporttaxonomy`, which runs ``ktImportTaxonomy`` on
one classification (a Kraken2 report, a BLAST table, ...) and writes the
interactive chart ``<prefix>.html``. ``taxonomy`` is KronaTools' taxonomy
database, a ``taxonomy.tab`` (``updateTaxonomy.sh`` builds one); it's staged
as ``taxonomy.tab`` in a directory of its own, as upstream.

``args`` carries ktImportTaxonomy's options, as upstream's ``ext.args``: for a
Kraken2 report, ``-t 5 -m 3`` reads the taxon ID from column 5 and the count
from column 3.
"""

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

KRONA_IMAGE = "quay.io/biocontainers/krona:2.8.1--pl5321hdfd78af_1"
DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory="2Gi")

# KronaTools is Perl, whose die exits 255, which the cluster takes for an
# interruption and retries without end; any failure exits 1 instead.
krona_ktimporttaxonomy_cmd = shell.create(
    name="krona_ktimporttaxonomy",
    image=KRONA_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"report": File, "taxonomy": File, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        REPORT=({inputs.report})
        TAX=({inputs.taxonomy})
        P={inputs.prefix}
        ARGS={inputs.args}
        OUT={outputs.results}
        W=$(mktemp -d); cd "$W"
        mkdir taxonomy
        ln -s "${TAX[0]}" taxonomy/taxonomy.tab
        ktImportTaxonomy \
            $ARGS \
            -o "$P.html" \
            -tax taxonomy/ \
            "${REPORT[0]}" \
            || exit 1
        mv "$P.html" "$OUT"/
    """,
)

env = flyte.TaskEnvironment.from_task("krona", krona_ktimporttaxonomy_cmd.as_task())


async def krona_ktimporttaxonomy(report: File, taxonomy: File, prefix: str = "sample", args: str = "") -> File:
    """The Krona chart ``<prefix>.html`` of ``report``'s taxa, using the ``taxonomy.tab`` ``taxonomy``."""
    results = await krona_ktimporttaxonomy_cmd(report=report, taxonomy=taxonomy, prefix=prefix, args=args)
    html = await results.get_file(f"{prefix}.html")
    if html is None:
        raise FileNotFoundError(f"ktImportTaxonomy wrote no {prefix}.html")
    return html
