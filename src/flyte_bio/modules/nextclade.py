"""nextclade — clade assignment, mutation calling and QC of viral genomes.

Exposes:

- :func:`nextclade_datasetget` — downloads a Nextclade dataset (reference,
  annotation, tree and QC configuration) by name, e.g.
  ``nextstrain/sars-cov-2/wuhan-hu-1/orfs``, at ``tag`` (a release such as
  ``2024-01-16--20-31-02Z``; the latest when empty). Without a tag the
  dataset changes over time.
- :func:`nextclade_run` — runs Nextclade on a FASTA of sequences against a
  dataset and returns every output under ``<prefix>``: ``.csv``, ``.tsv``,
  ``.json``, ``.ndjson``, ``.auspice.json``, ``.aligned.fasta``, ``.nwk``,
  ``_translation.<gene>.fasta`` and, when there are any, ``.errors.csv`` and
  ``.insertions.csv``. Its jobs follow the task's CPUs. The results record
  the dataset's directory name (``datasetName``), so the dataset is staged as
  ``dataset_name``: upstream's runs see the last component of the dataset's
  name (e.g. ``orfs``).
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

NEXTCLADE_IMAGE = "community.wave.seqera.io/library/nextclade:3.21.2--d3538cbe586c0f6c"
# Upstream: process_low.
DEFAULT_RESOURCES = flyte.Resources(cpu=2, memory="4Gi")

nextclade_datasetget_cmd = shell.create(
    name="nextclade_datasetget",
    image=NEXTCLADE_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"dataset": str, "tag": str, "args": str},
    defaults={"tag": "", "args": ""},
    outputs={"dataset_dir": Dir},
    script=r"""
        DATASET={inputs.dataset}
        TAG={inputs.tag}
        ARGS={inputs.args}
        VERSION=()
        if [ -n "$TAG" ]; then VERSION=(--tag "$TAG"); fi
        nextclade \
            dataset \
            get \
            $ARGS \
            --name "$DATASET" \
            "${VERSION[@]}" \
            --output-dir {outputs.dataset_dir} || exit 1
    """,
)

nextclade_run_cmd = shell.create(
    name="nextclade_run",
    image=NEXTCLADE_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"fasta": File, "dataset": Dir, "dataset_name": str, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        FA=({{inputs.fasta}})
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        ln -s "${{FA[0]}}" "$(basename "${{FA[0]}}")"
        DATASET={{inputs.dataset_name}}
        ln -s {{inputs.dataset}} "$DATASET"
        nextclade \
            run \
            $ARGS \
            --jobs $TASK_CPUS \
            --input-dataset "$DATASET" \
            --output-all "$OUT"/ \
            --output-basename "$P" \
            "$(basename "${{FA[0]}}")" || exit 1
    """,
)

env = flyte.TaskEnvironment.from_task("nextclade", nextclade_datasetget_cmd.as_task(), nextclade_run_cmd.as_task())


async def nextclade_datasetget(dataset: str, tag: str = "", args: str = "") -> Dir:
    """The Nextclade dataset ``dataset`` at release ``tag`` (the latest when empty)."""
    return await nextclade_datasetget_cmd(dataset=dataset, tag=tag, args=args)


@dataclass
class NextcladeRunResult:
    csv: File | None  # <prefix>.csv
    tsv: File | None  # <prefix>.tsv
    json: File | None  # <prefix>.json
    json_auspice: File | None  # <prefix>.auspice.json
    ndjson: File | None  # <prefix>.ndjson
    fasta_aligned: File | None  # <prefix>.aligned.fasta
    nwk: File | None  # <prefix>.nwk
    csv_errors: File | None  # <prefix>.errors.csv
    csv_insertions: File | None  # <prefix>.insertions.csv
    fasta_translation: list[File]  # *_translation.*.fasta, by name
    results: Dir  # everything nextclade wrote


async def nextclade_run(
    fasta: File, dataset: Dir, dataset_name: str = "dataset", prefix: str = "sample", args: str = ""
) -> NextcladeRunResult:
    """Nextclade's analysis of the sequences in ``fasta`` against ``dataset``, staged as ``dataset_name``."""
    if not dataset_name or "/" in dataset_name:
        raise ValueError(f"dataset_name must be a plain directory name, got {dataset_name!r}")
    results = await nextclade_run_cmd(fasta=fasta, dataset=dataset, dataset_name=dataset_name, prefix=prefix, args=args)
    files = {f.name: f for f in await results.list_files()}
    return NextcladeRunResult(
        csv=files.get(f"{prefix}.csv"),
        tsv=files.get(f"{prefix}.tsv"),
        json=files.get(f"{prefix}.json"),
        json_auspice=files.get(f"{prefix}.auspice.json"),
        ndjson=files.get(f"{prefix}.ndjson"),
        fasta_aligned=files.get(f"{prefix}.aligned.fasta"),
        nwk=files.get(f"{prefix}.nwk"),
        csv_errors=files.get(f"{prefix}.errors.csv"),
        csv_insertions=files.get(f"{prefix}.insertions.csv"),
        fasta_translation=[
            f for name, f in sorted(files.items()) if "_translation." in name and name.endswith(".fasta")
        ],
        results=results,
    )
