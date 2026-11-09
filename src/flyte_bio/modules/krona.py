"""krona — interactive taxonomy charts with KronaTools.

Exposes:

- :func:`krona_ktupdatetaxonomy`, which builds KronaTools' taxonomy database
  ``taxonomy.tab`` from NCBI's current taxonomy dump (``ktUpdateTaxonomy.sh``;
  it downloads a few hundred MB at run time). It's a cached task of its own,
  so the taxonomy is downloaded once; a fresh one needs a cache bust (a new
  ``args`` or a task override). It returns the single ``taxonomy.tab`` File,
  ready for :func:`krona_ktimporttaxonomy`.
- :func:`krona_ktimporttaxonomy`, which runs ``ktImportTaxonomy`` on
  one classification (a Kraken2 report, a BLAST table, ...) and writes the
  interactive chart ``<prefix>.html``. ``taxonomy`` is KronaTools' taxonomy
  database, a ``taxonomy.tab`` (``updateTaxonomy.sh`` builds one); it's staged
  as ``taxonomy.tab`` in a directory of its own, as upstream.
- :func:`krona_ktimporttext`, which runs ``ktImportText`` on one text file of
  counts and lineages (a count, then the lineage's names, tab-separated: as
  ``kaiju2krona`` and ``kreport2krona.py`` write) and writes the interactive
  chart ``<prefix>.html``. It needs no taxonomy database.

``args`` carries ktImportTaxonomy's (or ktImportText's) options, as
upstream's ``ext.args``: for a Kraken2 report, ``-t 5 -m 3`` reads the taxon
ID from column 5 and the count from column 3.
"""

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

KRONA_IMAGE = "quay.io/biocontainers/krona:2.8.1--pl5321hdfd78af_1"
# ktUpdateTaxonomy.sh builds the database with make, which upstream's image adds.
KRONA_MAKE_IMAGE = "community.wave.seqera.io/library/krona_make:7bb1fe2561793909"
DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory="2Gi")
# The NCBI dump is ~70 MB compressed and ~500 MB extracted.
UPDATE_TAXONOMY_RESOURCES = flyte.Resources(cpu=1, memory="4Gi", disk="10Gi")

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

krona_ktimporttext_cmd = shell.create(
    name="krona_ktimporttext",
    image=KRONA_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"report": File, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        REPORT=({inputs.report})
        P={inputs.prefix}
        ARGS={inputs.args}
        OUT={outputs.results}
        W=$(mktemp -d); cd "$W"
        ktImportText \
            $ARGS \
            -o "$P.html" \
            "${REPORT[0]}" \
            || exit 1
        mv "$P.html" "$OUT"/
    """,
)

krona_ktupdatetaxonomy_cmd = shell.create(
    name="krona_ktupdatetaxonomy",
    image=KRONA_MAKE_IMAGE,
    resources=UPDATE_TAXONOMY_RESOURCES,
    inputs={"args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        ARGS={inputs.args}
        OUT={outputs.results}
        W=$(mktemp -d); cd "$W"
        mkdir taxonomy
        ktUpdateTaxonomy.sh $ARGS taxonomy/ || exit 1
        [ -s taxonomy/taxonomy.tab ] || { echo "ktUpdateTaxonomy.sh wrote no taxonomy.tab" >&2; exit 1; }
        mv taxonomy/taxonomy.tab "$OUT"/
    """,
)

# One image per environment: ktUpdateTaxonomy.sh's adds make, so `env` has two.
krona_tools_env = flyte.TaskEnvironment.from_task(
    "krona_tools", krona_ktimporttaxonomy_cmd.as_task(), krona_ktimporttext_cmd.as_task()
)
krona_make_env = flyte.TaskEnvironment.from_task("krona_make", krona_ktupdatetaxonomy_cmd.as_task())
env = flyte.TaskEnvironment(name="krona", depends_on=[krona_tools_env, krona_make_env])


async def krona_ktupdatetaxonomy(args: str = "") -> File:
    """KronaTools' ``taxonomy.tab``, built from NCBI's current taxonomy dump."""
    results = await krona_ktupdatetaxonomy_cmd(args=args)
    taxonomy = await results.get_file("taxonomy.tab")
    if taxonomy is None:
        raise FileNotFoundError("ktUpdateTaxonomy.sh wrote no taxonomy.tab")
    return taxonomy


async def krona_ktimporttaxonomy(report: File, taxonomy: File, prefix: str = "sample", args: str = "") -> File:
    """The Krona chart ``<prefix>.html`` of ``report``'s taxa, using the ``taxonomy.tab`` ``taxonomy``."""
    results = await krona_ktimporttaxonomy_cmd(report=report, taxonomy=taxonomy, prefix=prefix, args=args)
    html = await results.get_file(f"{prefix}.html")
    if html is None:
        raise FileNotFoundError(f"ktImportTaxonomy wrote no {prefix}.html")
    return html


async def krona_ktimporttext(report: File, prefix: str = "sample", args: str = "") -> File:
    """The Krona chart ``<prefix>.html`` of the text file of counts and lineages ``report``."""
    results = await krona_ktimporttext_cmd(report=report, prefix=prefix, args=args)
    html = await results.get_file(f"{prefix}.html")
    if html is None:
        raise FileNotFoundError(f"ktImportText wrote no {prefix}.html")
    return html
