"""falint — fa-lint, validation of FASTA files.

Exposes :func:`falint`, which validates a FASTA and returns fa-lint's report:
``<prefix>.success.log`` when the FASTA is valid, otherwise
``<prefix>.error.log`` with fa-lint's errors. As upstream, an invalid FASTA
doesn't fail the task; check which log came back.

``-threads`` follows the task's CPUs. The FASTA is linked into the working
directory under its own name, as upstream stages it. Upstream tees fa-lint's
output into the logs through process substitutions; here fa-lint writes the
logs directly, which are then echoed to the task's stdout and stderr, so the
error check can't race the tee.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

FALINT_IMAGE = "quay.io/biocontainers/fa-lint:1.2.0--he881be0_0"
# Upstream: process_medium; a linter needs little.
DEFAULT_RESOURCES = flyte.Resources(cpu=2, memory="4Gi")

falint_cmd = shell.create(
    name="falint",
    image=FALINT_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"fasta": File, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        FA=({{inputs.fasta}})
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        W=$(mktemp -d); cd "$W"
        NAME=$(basename "${{FA[0]}}")
        ln -s "${{FA[0]}}" "$NAME"

        fa-lint \
            -threads $TASK_CPUS \
            $ARGS \
            -fasta "$NAME" \
            > "$P".success.log \
            2> "$P".error.log \
            || echo "Errors from fa-lint printed to $P.error.log"
        cat "$P".success.log
        cat "$P".error.log >&2

        if [ $(cat "$P".error.log | wc -l) -gt 0 ]; then
            echo "Validation failed..."
            rm "$P".success.log
            mv "$P".error.log {{outputs.results}}/
        else
            echo "Validation successful..."
            rm "$P".error.log
            mv "$P".success.log {{outputs.results}}/
        fi
    """,
)

env = flyte.TaskEnvironment.from_task("falint", falint_cmd.as_task())


@dataclass
class FalintResult:
    success_log: File | None  # <prefix>.success.log, when the FASTA is valid
    error_log: File | None  # <prefix>.error.log, when it isn't


async def falint(fasta: File, prefix: str = "sample", args: str = "") -> FalintResult:
    """Validate ``fasta``: exactly one of the result's logs is set."""
    results = await falint_cmd(fasta=fasta, prefix=prefix, args=args)
    return FalintResult(
        success_log=await results.get_file(f"{prefix}.success.log"),
        error_log=await results.get_file(f"{prefix}.error.log"),
    )
