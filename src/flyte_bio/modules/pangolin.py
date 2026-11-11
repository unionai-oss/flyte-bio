"""pangolin — Pango lineage assignment of SARS-CoV-2 genomes.

Exposes:

- :func:`pangolin_updatedata` — installs the latest pangolin-data (the
  lineage models and designations) into a data directory. The data is live,
  so it changes over time; when the image's own pangolin-data is already
  the latest, the directory stays empty, as upstream.
- :func:`pangolin_run` — assigns lineages to the genomes in a FASTA:
  ``<prefix>.pangolin.csv``. With a data directory from
  :func:`pangolin_updatedata` it uses that data, else the image's. Its
  threads follow the task's CPUs.
"""

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

PANGOLIN_IMAGE = "community.wave.seqera.io/library/pangolin-data_pangolin_snakemake-minimal:638a1eb68adff9c7"
# Upstream: updatedata process_single, run process_medium.
UPDATEDATA_RESOURCES = flyte.Resources(cpu=1, memory="4Gi")
RUN_RESOURCES = flyte.Resources(cpu=2, memory="8Gi")

pangolin_updatedata_cmd = shell.create(
    name="pangolin_updatedata",
    image=PANGOLIN_IMAGE,
    resources=UPDATEDATA_RESOURCES,
    inputs={"args": str},
    defaults={"args": ""},
    outputs={"db": Dir},
    script=rf"""
        {TASK_LIMITS}
        ARGS={{inputs.args}}
        export XDG_CACHE_HOME=/tmp/.cache
        pangolin \
            $ARGS \
            --update-data \
            --threads $TASK_CPUS \
            --datadir {{outputs.db}} || exit 1
    """,
)

pangolin_run_cmd = shell.create(
    name="pangolin_run",
    image=PANGOLIN_IMAGE,
    resources=RUN_RESOURCES,
    inputs={"fasta": File, "db": Dir | None, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        FA=({{inputs.fasta}})
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        # An unset optional Dir is staged as a file, not a directory.
        DB={{inputs.db}}
        DATADIR=()
        if [ -d "$DB" ]; then DATADIR=(--datadir "$DB"); fi
        export XDG_CACHE_HOME=/tmp/.cache
        W=$(mktemp -d); cd "$W"
        ln -s "${{FA[0]}}" "$(basename "${{FA[0]}}")"
        pangolin \
            "$(basename "${{FA[0]}}")" \
            "${{DATADIR[@]}}" \
            --outfile "$P.pangolin.csv" \
            --threads $TASK_CPUS \
            $ARGS || exit 1
        mv "$P.pangolin.csv" "$OUT"/
    """,
)

env = flyte.TaskEnvironment.from_task("pangolin", pangolin_updatedata_cmd.as_task(), pangolin_run_cmd.as_task())


async def pangolin_updatedata(args: str = "") -> Dir:
    """A data directory holding the latest pangolin-data (empty when the image's is the latest)."""
    return await pangolin_updatedata_cmd(args=args)


async def pangolin_run(fasta: File, db: Dir | None = None, prefix: str = "sample", args: str = "") -> File:
    """Lineage assignments of the genomes in ``fasta``: ``<prefix>.pangolin.csv``."""
    results = await pangolin_run_cmd(fasta=fasta, db=db, prefix=prefix, args=args)
    report = await results.get_file(f"{prefix}.pangolin.csv")
    if report is None:
        raise FileNotFoundError(f"pangolin wrote no {prefix}.pangolin.csv")
    return report
