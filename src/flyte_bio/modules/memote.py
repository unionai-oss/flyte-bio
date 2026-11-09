"""memote — MEMOTE, standardised quality tests for genome-scale metabolic models.

Exposes:

- :func:`memote_run` — run MEMOTE's test suite on an SBML model
  (``memote run``): the results as ``<prefix>.json.gz``.
- :func:`memote_report` — MEMOTE's snapshot report of an SBML model
  (``memote report snapshot``): ``<prefix>.html``.

As upstream, both run with ``HOME`` in the working directory and the GLPK
solver (``COBRA_SOLVER=glpk_exact``, open source), so no solver licence is
needed; ``memote run`` ignores git (``--ignore-git``).
"""

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

MEMOTE_IMAGE = "quay.io/biocontainers/memote:0.17.0--pyhdfd78af_0"
# Upstream's process_low.
DEFAULT_RESOURCES = flyte.Resources(cpu=2, memory="8Gi")

memote_run_cmd = shell.create(
    name="memote_run",
    image=MEMOTE_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"model": File, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        MODEL=({inputs.model})
        P={inputs.prefix}
        ARGS={inputs.args}
        OUT={outputs.results}
        W=$(mktemp -d); cd "$W"
        export HOME="$W"
        export COBRA_SOLVER=glpk_exact
        memote run --filename "$P.json.gz" --ignore-git $ARGS "${MODEL[0]}" || exit 1
        mv "$P.json.gz" "$OUT"/
    """,
)

memote_report_cmd = shell.create(
    name="memote_report",
    image=MEMOTE_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"model": File, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        MODEL=({inputs.model})
        P={inputs.prefix}
        ARGS={inputs.args}
        OUT={outputs.results}
        W=$(mktemp -d); cd "$W"
        export HOME="$W"
        export COBRA_SOLVER=glpk_exact
        memote report snapshot --filename "$P.html" $ARGS "${MODEL[0]}" || exit 1
        mv "$P.html" "$OUT"/
    """,
)

env = flyte.TaskEnvironment.from_task("memote", memote_run_cmd.as_task(), memote_report_cmd.as_task())


async def pick(results: Dir, name: str) -> File:
    found = await results.get_file(name)
    if found is None:
        raise FileNotFoundError(f"memote wrote no {name}")
    return found


async def memote_run(model: File, prefix: str = "sample", args: str = "") -> File:
    """MEMOTE's test results for the SBML ``model``: ``<prefix>.json.gz``."""
    results = await memote_run_cmd(model=model, prefix=prefix, args=args)
    return await pick(results, f"{prefix}.json.gz")


async def memote_report(model: File, prefix: str = "sample", args: str = "") -> File:
    """MEMOTE's snapshot report of the SBML ``model``: ``<prefix>.html``."""
    results = await memote_report_cmd(model=model, prefix=prefix, args=args)
    return await pick(results, f"{prefix}.html")
