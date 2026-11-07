"""sourcepredict — predict the sources of a metagenomic sample (SourcePredict).

Exposes :func:`sourcepredict`, which embeds a sink sample's taxon counts
(``kraken_parse``: a TAXID count table, one column per sample) among labelled
source samples (``sources``, ``labels``) and predicts each sink sample's
proportion of every source and of an unknown source. It returns
``<prefix>.sourcepredict.csv``: one row per source (and ``unknown``), one
column per sink sample. ``-t`` is the task's CPUs.

SourcePredict reads the NCBI taxonomy through ete3, from
``~/.etetoolkit/taxa.sqlite`` and its ``taxa.sqlite.traverse.pkl``.
``taxa_sqlite`` and ``taxa_sqlite_traverse_pkl`` are linked there under those
names (upstream stages them into ``.etetoolkit/`` under their own names), so
ete3 uses them instead of downloading the taxonomy.
"""

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

SOURCEPREDICT_IMAGE = "quay.io/biocontainers/sourcepredict:0.5.1--pyhdfd78af_0"
# Upstream's process_medium is 6 CPUs and 36 GB; the test data needs far less.
DEFAULT_RESOURCES = flyte.Resources(cpu=4, memory="8Gi", disk="20Gi")

sourcepredict_cmd = shell.create(
    name="sourcepredict",
    image=SOURCEPREDICT_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "kraken_parse": File,
        "sources": File,
        "labels": File,
        "taxa_sqlite": File,
        "taxa_sqlite_traverse_pkl": File,
        "prefix": str,
        "args": str,
    },
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        SINK=({{inputs.kraken_parse}})
        SOURCES=({{inputs.sources}})
        LABELS=({{inputs.labels}})
        SQLITE=({{inputs.taxa_sqlite}})
        TRAVERSE=({{inputs.taxa_sqlite_traverse_pkl}})
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        W=$(mktemp -d); cd "$W"
        export HOME="$W" NUMBA_CACHE_DIR="$W/tmp" MPLCONFIGDIR="$W/tmp"
        mkdir .etetoolkit
        ln -s "${{SQLITE[0]}}" .etetoolkit/taxa.sqlite
        ln -s "${{TRAVERSE[0]}}" .etetoolkit/taxa.sqlite.traverse.pkl
        sourcepredict \
            -s "${{SOURCES[0]}}" \
            -l "${{LABELS[0]}}" \
            $ARGS \
            -t $TASK_CPUS \
            -o "$P.sourcepredict.csv" \
            "${{SINK[0]}}" || exit 1
        mv "$P.sourcepredict.csv" {{outputs.results}}/
    """,
)

env = flyte.TaskEnvironment.from_task("sourcepredict", sourcepredict_cmd.as_task())


async def sourcepredict(
    kraken_parse: File,
    sources: File,
    labels: File,
    taxa_sqlite: File,
    taxa_sqlite_traverse_pkl: File,
    prefix: str = "sample",
    args: str = "",
) -> File:
    """SourcePredict's source proportions for the sink samples in ``kraken_parse``: ``<prefix>.sourcepredict.csv``."""
    results = await sourcepredict_cmd(
        kraken_parse=kraken_parse,
        sources=sources,
        labels=labels,
        taxa_sqlite=taxa_sqlite,
        taxa_sqlite_traverse_pkl=taxa_sqlite_traverse_pkl,
        prefix=prefix,
        args=args,
    )
    name = f"{prefix}.sourcepredict.csv"
    out = await results.get_file(name)
    if out is None:
        raise FileNotFoundError(f"sourcepredict wrote no {name}")
    return out
