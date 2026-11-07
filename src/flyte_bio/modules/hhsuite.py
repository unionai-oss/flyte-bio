"""hhsuite — HH-suite utilities.

Exposes :func:`hhsuite_reformat`, which converts the multiple alignment
``aln`` (optionally gzipped, as ``.gz``; decompressed first, as upstream)
from ``informat`` to ``outformat`` with HH-suite's ``reformat.pl`` (formats
such as ``fas``, ``a2m``, ``a3m``, ``sto``, ``psi``, ``clu``). The result is
``<prefix>.<outformat>.gz``. ``reformat.pl``'s options (``-M first``,
``-r``) go in ``args``.
"""

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

HHSUITE_IMAGE = "quay.io/biocontainers/hhsuite:3.3.0--py39pl5321h0dd7abe_13"
# Upstream's process_single.
DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory="2Gi")

hhsuite_reformat_cmd = shell.create(
    name="hhsuite_reformat",
    image=HHSUITE_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"aln": File, "informat": str, "outformat": str, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        ALN=({inputs.aln})
        INFORMAT={inputs.informat}
        OUTFORMAT={inputs.outformat}
        P={inputs.prefix}
        ARGS={inputs.args}
        OUT={outputs.results}
        W=$(mktemp -d); cd "$W"

        NAME=$(basename "${ALN[0]}")
        if [ "${NAME%.gz}" != "$NAME" ]; then
            gzip -c -d "${ALN[0]}" > "${NAME%.gz}" || exit 1
            INPUT="${NAME%.gz}"
        else
            INPUT="${ALN[0]}"
        fi
        reformat.pl $ARGS "$INFORMAT" "$OUTFORMAT" "$INPUT" "$P.$OUTFORMAT" || exit 1
        gzip "$P.$OUTFORMAT"
        mv "$P.$OUTFORMAT.gz" "$OUT"/
    """,
)

env = flyte.TaskEnvironment.from_task("hhsuite", hhsuite_reformat_cmd.as_task())


async def hhsuite_reformat(aln: File, informat: str, outformat: str, prefix: str = "sample", args: str = "") -> File:
    """``aln`` converted from ``informat`` to ``outformat``: ``<prefix>.<outformat>.gz``."""
    results = await hhsuite_reformat_cmd(aln=aln, informat=informat, outformat=outformat, prefix=prefix, args=args)
    name = f"{prefix}.{outformat}.gz"
    found = await results.get_file(name)
    if found is None:
        raise FileNotFoundError(f"reformat.pl wrote no {name}")
    return found
