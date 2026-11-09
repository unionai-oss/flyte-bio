"""gmmdemux — demultiplex cell-hashing data with GMM-Demux.

Exposes :func:`gmmdemux`, which reads a 10x HTO count matrix directory and
classifies each droplet from the named hashtags (``hto_names``, comma
separated) with Gaussian mixture models. Outputs keep upstream's names:

- ``barcodes.tsv.gz``, ``features.tsv.gz``, ``matrix.mtx.gz``: the matrix of
  the droplets GMM-Demux keeps;
- ``GMM_*.csv`` and ``GMM_*.config``: the classification and its class key,
  full (``full_report``, GMM-Demux's ``-f``) or simplified (``-s``);
- with ``summary_report``, ``<prefix>_summary_report.txt`` (GMM-Demux's
  ``-r``; needs ``--summary <cells>`` in ``args``).

``args`` go to GMM-Demux (e.g. ``-t 0.8 --summary 11089 -a 0.05``);
``skip`` and ``examine`` are its ``--skip`` and ``--examine`` cell lists.
Upstream's output pattern for the summary report (``summary_report_*.txt``)
doesn't match the name it writes, so upstream never emits it; here it's
returned.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

GMMDEMUX_IMAGE = "quay.io/biocontainers/gmm-demux:0.2.2.3--pyh7cba7a3_0"
# Upstream: process_low.
DEFAULT_RESOURCES = flyte.Resources(cpu=2, memory="8Gi")

gmmdemux_cmd = shell.create(
    name="gmmdemux",
    image=GMMDEMUX_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "hto_matrix": Dir,
        "hto_names": str,
        "full_report": bool,
        "summary_report": bool,
        "skip": File | None,
        "examine": File | None,
        "prefix": str,
        "args": str,
    },
    defaults={"full_report": False, "summary_report": False, "args": ""},
    outputs={"results": Dir},
    script=r"""
        MATRIX={inputs.hto_matrix}
        shopt -s nullglob; SKIP=({inputs.skip}); EXAMINE=({inputs.examine}); shopt -u nullglob
        P={inputs.prefix}
        ARGS={inputs.args}
        OUT={outputs.results}
        W=$(mktemp -d); cd "$W"
        OPTS=()
        if [ {inputs.full_report} = true ]; then OPTS+=(-f .); else OPTS+=(-s .); fi
        if [ {inputs.summary_report} = true ]; then
            cat /dev/null > "${P}_summary_report.txt"
            OPTS+=(-r "${P}_summary_report.txt")
        fi
        if [ ${#SKIP[@]} -gt 0 ]; then OPTS+=(--skip "${SKIP[0]}"); fi
        if [ ${#EXAMINE[@]} -gt 0 ]; then OPTS+=(--examine "${EXAMINE[0]}"); fi
        GMM-demux $ARGS "${OPTS[@]}" "$MATRIX" {inputs.hto_names} -o . || exit 1
        mv barcodes.tsv.gz features.tsv.gz matrix.mtx.gz GMM_*.csv GMM_*.config "$OUT"/
        if [ -f "${P}_summary_report.txt" ]; then mv "${P}_summary_report.txt" "$OUT"/; fi
    """,
)

env = flyte.TaskEnvironment.from_task("gmmdemux", gmmdemux_cmd.as_task())


@dataclass
class GmmDemuxResult:
    barcodes: File  # barcodes.tsv.gz
    features: File  # features.tsv.gz
    matrix: File  # matrix.mtx.gz
    classification_report: File  # GMM_*.csv
    config_report: File  # GMM_*.config
    summary_report: File | None  # <prefix>_summary_report.txt
    results: Dir


async def gmmdemux(
    hto_matrix: Dir,
    hto_names: str,
    prefix: str,
    full_report: bool = False,
    summary_report: bool = False,
    skip: File | None = None,
    examine: File | None = None,
    args: str = "",
) -> GmmDemuxResult:
    """GMM-Demux classification of a 10x HTO matrix directory by the hashtags in ``hto_names`` (comma separated)."""
    results = await gmmdemux_cmd(
        hto_matrix=hto_matrix,
        hto_names=hto_names,
        full_report=full_report,
        summary_report=summary_report,
        skip=skip,
        examine=examine,
        prefix=prefix,
        args=args,
    )

    async def pick(name: str) -> File:
        found = await results.get_file(name)
        if found is None:
            raise FileNotFoundError(f"gmmdemux wrote no {name}")
        return found

    async def pick_report(suffix: str) -> File:
        found = [f async for f in results.walk() if f.name.startswith("GMM_") and f.name.endswith(suffix)]
        if len(found) != 1:
            raise FileNotFoundError(f"gmmdemux wrote {len(found)} GMM_*{suffix} files, expected one")
        return found[0]

    return GmmDemuxResult(
        barcodes=await pick("barcodes.tsv.gz"),
        features=await pick("features.tsv.gz"),
        matrix=await pick("matrix.mtx.gz"),
        classification_report=await pick_report(".csv"),
        config_report=await pick_report(".config"),
        summary_report=await results.get_file(f"{prefix}_summary_report.txt") if summary_report else None,
        results=results,
    )
