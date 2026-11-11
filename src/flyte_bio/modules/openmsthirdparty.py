"""openmsthirdparty — OpenMS adapters to third-party search engines.

Exposes :func:`openmsthirdparty_cometadapter`, which searches an mzML's
spectra against a protein FASTA with Comet through OpenMS' ``CometAdapter``
and returns the identifications, ``<prefix>.idXML``, and Percolator input,
``<prefix>.tsv``, when ``args`` ask for it (``-pin_out``).

As upstream, Comet gets the task's CPUs as ``-threads`` and reads its inputs
linked into the working directory under their own names, which the idXML
records. ``args`` is parsed with shell quoting, so values with spaces can be
quoted as upstream writes them (``-enzyme 'unspecific cleavage'``).
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.modules.openms import OPENMS_PREAMBLE

OPENMS_THIRDPARTY_IMAGE = "community.wave.seqera.io/library/openms-thirdparty:3.6.0--566f6f290d39203b"
# Upstream labels CometAdapter process_high.
HIGH_RESOURCES = flyte.Resources(cpu=8, memory="16Gi")

openmsthirdparty_cometadapter_cmd = shell.create(
    name="openmsthirdparty_cometadapter",
    image=OPENMS_THIRDPARTY_IMAGE,
    resources=HIGH_RESOURCES,
    inputs={"mzml": File, "fasta": File, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=OPENMS_PREAMBLE
    + r"""
        MZML=({inputs.mzml}); FA=({inputs.fasta}); stage "${MZML[0]}" "${FA[0]}"
        P={inputs.prefix}
        CometAdapter -in "$(basename "${MZML[0]}")" -database "$(basename "${FA[0]}")" -out "$P.idXML" \
            -threads $TASK_CPUS "${ARGS[@]}" || exit 1
        shopt -s nullglob
        mv ./*.idXML ./*.tsv "$OUT"/
    """,
)

env = flyte.TaskEnvironment.from_task("openmsthirdparty", openmsthirdparty_cometadapter_cmd.as_task())


@dataclass
class CometAdapterResult:
    idxml: File  # <prefix>.idXML
    pin: File | None  # <prefix>.tsv, Percolator input, when args ask for it
    results: Dir  # everything written


async def comet_result(results: Dir, prefix: str) -> CometAdapterResult:
    """Pick CometAdapter's outputs out of its ``results``."""
    idxml = await results.get_file(f"{prefix}.idXML")
    if idxml is None:
        raise FileNotFoundError(f"CometAdapter wrote no {prefix}.idXML")
    return CometAdapterResult(idxml=idxml, pin=await results.get_file(f"{prefix}.tsv"), results=results)


async def openmsthirdparty_cometadapter(
    mzml: File, fasta: File, prefix: str = "sample", args: str = ""
) -> CometAdapterResult:
    """Search ``mzml`` against ``fasta`` with Comet: ``<prefix>.idXML`` (and ``<prefix>.tsv``)."""
    results = await openmsthirdparty_cometadapter_cmd(mzml=mzml, fasta=fasta, prefix=prefix, args=args)
    return await comet_result(results, prefix)
