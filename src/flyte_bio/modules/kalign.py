"""kalign — multiple sequence alignment with Kalign.

Exposes :func:`kalign_align`, which aligns the sequences in ``fasta`` (may be
gzipped) and writes ``<prefix>.aln``, or ``<prefix>.aln.gz`` (gzipped with
pigz) with ``compress``. Kalign runs with ``-n`` set to the task's CPUs; its
other options go in ``args``.

As upstream, a Kalign crash with exit code 132 (an illegal instruction: the
build doesn't support some CPU types) is reported as such. Upstream writes the
compressed alignment through a process substitution; here Kalign writes the
file and pigz compresses it after, which is the same alignment.
"""

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

# Upstream's Docker image: Kalign 3.4.0 and pigz.
KALIGN_IMAGE = (
    "quay.io/biocontainers/mulled-v2-5cd0277547c6b33133225c8ce14c0cf2a4396ea2:"
    "0a70b6d89a3e06fbdc4a735461e8b98ff32ee5de-0"
)
# Upstream labels this process_medium; its tests run with 4 CPUs.
DEFAULT_RESOURCES = flyte.Resources(cpu=4, memory="8Gi")

kalign_align_cmd = shell.create(
    name="kalign_align",
    image=KALIGN_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"fasta": File, "prefix": str, "compress": bool, "args": str},
    defaults={"compress": False, "args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        set -o pipefail
        FA=({{inputs.fasta}})
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"

        unpigz -cdf "${{FA[0]}}" | \
        kalign \
            $ARGS \
            -n $TASK_CPUS \
            -o "$P.aln" || {{
            CODE=$?
            if [ $CODE -eq 132 ]; then
                echo "KALIGN failed because is incompatible with some CPU types," \
                    "see https://github.com/TimoLassmann/kalign/issues/46." >&2
            fi
            exit 1
        }}
        if [ {{inputs.compress}} = true ]; then
            pigz -p $TASK_CPUS "$P.aln"
        fi
        mv "$P".aln* "$OUT"/
    """,
)

env = flyte.TaskEnvironment.from_task("kalign", kalign_align_cmd.as_task())


async def kalign_align(fasta: File, prefix: str = "sample", compress: bool = False, args: str = "") -> File:
    """Align ``fasta``: ``<prefix>.aln``, or ``<prefix>.aln.gz`` with ``compress``."""
    results = await kalign_align_cmd(fasta=fasta, prefix=prefix, compress=compress, args=args)
    name = f"{prefix}.aln.gz" if compress else f"{prefix}.aln"
    found = await results.get_file(name)
    if found is None:
        raise FileNotFoundError(f"kalign wrote no {name}")
    return found
