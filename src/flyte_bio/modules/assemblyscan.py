"""assemblyscan — summary statistics of an assembly (assembly-scan).

Exposes :func:`assemblyscan`: contig counts, lengths, N50 and GC content of an
assembly FASTA, as ``<prefix>.tsv``, or ``<prefix>.json`` when ``args`` has
``--json`` (as upstream, which picks the extension the same way).

The assembly is linked into the working directory under its own name, which
the report records.
"""

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

ASSEMBLYSCAN_IMAGE = "community.wave.seqera.io/library/assembly-scan:1.0.0--fdab3a3c6ba9bbd0"
# Upstream's process_low.
ASSEMBLYSCAN_RESOURCES = flyte.Resources(cpu=1, memory="2Gi")

assemblyscan_cmd = shell.create(
    name="assemblyscan",
    image=ASSEMBLYSCAN_IMAGE,
    resources=ASSEMBLYSCAN_RESOURCES,
    inputs={"assembly": File, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        set -e
        ASM=({inputs.assembly})
        P={inputs.prefix}
        ARGS={inputs.args}
        FORMAT=tsv
        case " $ARGS " in *--json*) FORMAT=json ;; esac
        W=$(mktemp -d); cd "$W"
        NAME=$(basename "${ASM[0]}")
        ln -s "${ASM[0]}" "$NAME"
        assembly-scan $ARGS "$NAME" > "$P.$FORMAT" || exit 1
        mv "$P.$FORMAT" {outputs.results}/
    """,
)

env = flyte.TaskEnvironment.from_task("assemblyscan", assemblyscan_cmd.as_task())


async def assemblyscan(assembly: File, prefix: str = "sample", args: str = "") -> File:
    """assembly-scan report of ``assembly``: ``<prefix>.tsv``, or ``<prefix>.json`` with ``--json`` in ``args``."""
    name = f"{prefix}.{'json' if '--json' in args else 'tsv'}"
    results = await assemblyscan_cmd(assembly=assembly, prefix=prefix, args=args)
    report = await results.get_file(name)
    if report is None:
        raise FileNotFoundError(f"assembly-scan wrote no {name}")
    return report
