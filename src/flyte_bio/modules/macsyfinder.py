"""macsyfinder — MacSyFinder, detection of macromolecular systems (secretion systems, pili, …) in proteomes.

Exposes:

- :func:`macsyfinder_download` — install a model package from the MacSyModels
  repository (``msf_data install``, e.g. ``TXSScan``) and return it as
  ``models.tar.gz``, an archive of the ``models`` folder. It's an archive
  rather than a Dir because a package holds hundreds of profiles, too many for
  a Dir to stage reliably.
- :func:`macsyfinder_search` — search a protein FASTA (plain or gzipped) for
  the systems of a model package (``macsyfinder``). Its output folder's files
  keep their names: ``all_systems.tsv`` (the summary), ``all_best_solutions*``
  (ordered replicons only), MacSyFinder's other tables, its standard output
  and error as ``macsyfinder.out`` and ``macsyfinder.err``, and the per-profile
  HMMER results, archived as ``hmmer_results.tar.gz``.

As upstream, the search gets the task's CPUs as workers (``--worker``).
"""

import os
from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

MACSYFINDER_IMAGE = "quay.io/biocontainers/macsyfinder:2.1.6--pyhdfd78af_0"
# Upstream's process_single and process_medium, scaled down.
DOWNLOAD_RESOURCES = flyte.Resources(cpu=1, memory="2Gi")
SEARCH_RESOURCES = flyte.Resources(cpu=4, memory="8Gi")

macsyfinder_download_cmd = shell.create(
    name="macsyfinder_download",
    image=MACSYFINDER_IMAGE,
    resources=DOWNLOAD_RESOURCES,
    inputs={"model_name": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        MODEL={inputs.model_name}
        ARGS={inputs.args}
        OUT={outputs.results}
        W=$(mktemp -d); cd "$W"
        export HOME="$W"
        mkdir -p models
        msf_data install --target models $ARGS "$MODEL" || exit 1
        tar -czf "$OUT/models.tar.gz" models
    """,
)

macsyfinder_search_cmd = shell.create(
    name="macsyfinder_search",
    image=MACSYFINDER_IMAGE,
    resources=SEARCH_RESOURCES,
    inputs={"proteins": File, "models": File, "model_names": str, "prefix": str, "args": str},
    defaults={"model_names": "", "args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        PROTEINS=({{inputs.proteins}}); MODELS=({{inputs.models}})
        MODEL_NAMES={{inputs.model_names}}
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        export HOME="$W"
        mkdir db
        tar -xzf "${{MODELS[0]}}" -C db
        MODEL_ARG=()
        if [ -n "$MODEL_NAMES" ]; then MODEL_ARG=(--models $MODEL_NAMES); fi
        ln -s "${{PROTEINS[0]}}" "$(basename "${{PROTEINS[0]}}")"
        macsyfinder --sequence-db "$(basename "${{PROTEINS[0]}}")" --models-dir db/models "${{MODEL_ARG[@]}}" \
            --out-dir "$P" --worker $TASK_CPUS $ARGS > macsyfinder.out 2> macsyfinder.err \
            || {{ tail -n 20 macsyfinder.err >&2; exit 1; }}
        mv macsyfinder.out macsyfinder.err "$P"/
        if [ -d "$P/hmmer_results" ]; then
            tar -czf "$OUT/hmmer_results.tar.gz" -C "$P" hmmer_results
            rm -r "$P/hmmer_results"
        fi
        find "$P" -mindepth 1 -maxdepth 1 -type f -exec mv {{}} "$OUT"/ \;
    """,
)

env = flyte.TaskEnvironment.from_task(
    "macsyfinder", macsyfinder_download_cmd.as_task(), macsyfinder_search_cmd.as_task()
)


async def macsyfinder_download(model_name: str, args: str = "") -> File:
    """The MacSyModels package ``model_name`` (e.g. ``TXSScan``): ``models.tar.gz``."""
    results = await macsyfinder_download_cmd(model_name=model_name, args=args)
    models = await results.get_file("models.tar.gz")
    if models is None:
        raise FileNotFoundError("msf_data wrote no models.tar.gz")
    return models


@dataclass
class MacsyfinderResult:
    summary: File | None  # all_systems.tsv
    best_solutions: list[File]  # all_best_solutions*
    stdout: File | None  # macsyfinder.out
    stderr: File | None  # macsyfinder.err
    hmmer: File | None  # hmmer_results.tar.gz
    results: Dir  # everything MacSyFinder wrote


async def macsyfinder_search(
    proteins: File, models: File, model_names: str = "", prefix: str = "sample", args: str = ""
) -> MacsyfinderResult:
    """Search ``proteins`` for the systems in ``models`` (from :func:`macsyfinder_download`)."""
    results = await macsyfinder_search_cmd(
        proteins=proteins, models=models, model_names=model_names, prefix=prefix, args=args
    )
    files = await results.list_files()
    return MacsyfinderResult(
        summary=await results.get_file("all_systems.tsv"),
        best_solutions=sorted(
            (f for f in files if os.path.basename(f.path).startswith("all_best_solutions")),
            key=lambda f: os.path.basename(f.path),
        ),
        stdout=await results.get_file("macsyfinder.out"),
        stderr=await results.get_file("macsyfinder.err"),
        hmmer=await results.get_file("hmmer_results.tar.gz"),
        results=results,
    )
