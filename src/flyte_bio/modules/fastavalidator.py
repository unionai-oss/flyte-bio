"""fastavalidator — FASTA format validation with py_fasta_validator.

Exposes :func:`fastavalidator`, which checks that ``fasta`` is a well-formed
FASTA file. As upstream, a failed validation doesn't fail the task: it writes
the validator's messages to ``<prefix>.error.log``. A valid file gets
``<prefix>.success.log`` ("Validation successful...") instead.

Upstream has since deprecated this module (in favour of fa-lint, see
:mod:`flyte_bio.modules.falint`) and made it fail on purpose; this is its last
working version (modules commit 41dfa3f).
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

FASTAVALIDATOR_IMAGE = "quay.io/biocontainers/py_fasta_validator:0.6--py37h595c7a6_0"
# Upstream labels this process_single.
DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory="2Gi")

fastavalidator_cmd = shell.create(
    name="fastavalidator",
    image=FASTAVALIDATOR_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"fasta": File, "prefix": str},
    outputs={"results": Dir},
    script=r"""
        FA=({inputs.fasta})
        P={inputs.prefix}
        OUT={outputs.results}
        W=$(mktemp -d); cd "$W"

        # The validator's messages name the file, so it runs on the file's own name.
        NAME=$(basename "${FA[0]}")
        ln -s "${FA[0]}" "$NAME"
        py_fasta_validator \
            -f "$NAME" \
            2> "$P.error.log" \
            || echo "Errors from fasta_validate printed to $P.error.log"

        if [ $(cat "$P.error.log" | wc -l) -gt 0 ]; then
            echo "Validation failed..."

            cat \
                "$P.error.log"
            mv "$P.error.log" "$OUT"/
        else
            echo "Validation successful..."

            mv \
                "$P.error.log" \
                fasta_validate.stderr

            echo "Validation successful..." \
                > "$OUT/$P.success.log"
        fi
    """,
)

env = flyte.TaskEnvironment.from_task("fastavalidator", fastavalidator_cmd.as_task())


@dataclass
class FastaValidatorResult:
    success_log: File | None  # <prefix>.success.log, when the FASTA is valid
    error_log: File | None  # <prefix>.error.log, the validator's messages, when it isn't


async def fastavalidator(fasta: File, prefix: str = "sample") -> FastaValidatorResult:
    """Validate ``fasta``: a success log if it's valid, else the error log."""
    results = await fastavalidator_cmd(fasta=fasta, prefix=prefix)
    return FastaValidatorResult(
        success_log=await results.get_file(f"{prefix}.success.log"),
        error_log=await results.get_file(f"{prefix}.error.log"),
    )
