"""msisensorpro — MSIsensor-pro microsatellite tools.

Exposes :func:`msisensorpro_scan`, which scans a reference FASTA for
homopolymers and microsatellites and writes the site list, as upstream, to
``<prefix>.msisensor_scan.list``.
"""

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

MSISENSORPRO_IMAGE = "quay.io/biocontainers/msisensor-pro:1.3.0--hfef96ef_0"
# Upstream labels scan process_low (2 CPUs, 12 GB); scan is single-threaded and
# a chromosome-sized FASTA needs far less memory.
SCAN_RESOURCES = flyte.Resources(cpu=1, memory="4Gi")

msisensorpro_scan_cmd = shell.create(
    name="msisensorpro_scan",
    image=MSISENSORPRO_IMAGE,
    resources=SCAN_RESOURCES,
    inputs={"fasta": File, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        FA=({inputs.fasta})
        P={inputs.prefix}
        ARGS={inputs.args}
        W=$(mktemp -d); cd "$W"
        FA_NAME=$(basename "${FA[0]}"); ln -s "${FA[0]}" "$FA_NAME"
        msisensor-pro scan -d "$FA_NAME" -o {outputs.results}/"$P".msisensor_scan.list $ARGS
    """,
)

env = flyte.TaskEnvironment.from_task("msisensorpro", msisensorpro_scan_cmd.as_task())


async def msisensorpro_scan(fasta: File, prefix: str, args: str = "") -> File:
    """``<prefix>.msisensor_scan.list``: the microsatellite sites MSIsensor-pro finds in ``fasta``."""
    results = await msisensorpro_scan_cmd(fasta=fasta, prefix=prefix, args=args)
    sites = await results.get_file(f"{prefix}.msisensor_scan.list")
    if sites is None:
        raise FileNotFoundError(f"msisensor-pro scan wrote no {prefix}.msisensor_scan.list")
    return sites
