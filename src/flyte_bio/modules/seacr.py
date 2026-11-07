"""seacr — SEACR peak calling from CUT&RUN / CUT&Tag bedGraphs.

Exposes :func:`seacr_callpeak`, which calls enriched regions in a target
bedGraph either against a control (IgG) bedGraph or, without one, as the top
``threshold`` fraction of regions by signal. As upstream, the mode
(``norm``/``non`` and ``stringent``/``relaxed``) goes in ``args``, e.g.
``"norm stringent"``. The peaks are ``<prefix>.<stringent|relaxed>.bed``.
"""

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

# SEACR 1.3 with bedtools 2.30.0 and R 4.0.2.
SEACR_IMAGE = (
    "quay.io/biocontainers/mulled-v2-03bfeb32fe80910c231f630d4262b83677c8c0f4:"
    "f4bb19b68e66de27e4c64306f951d5ff11919931-0"
)
# Upstream: process_medium. SEACR itself is single-threaded.
DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory="6Gi")

seacr_callpeak_cmd = shell.create(
    name="seacr_callpeak",
    image=SEACR_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"bedgraph": File, "control": File | None, "threshold": float, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        BG=({inputs.bedgraph})
        shopt -s nullglob; CTRL=({inputs.control}); shopt -u nullglob
        P={inputs.prefix}
        ARGS={inputs.args}
        OUT={outputs.results}
        W=$(mktemp -d); cd "$W"

        # As upstream, the control bedGraph if given, else the threshold.
        SWITCH={inputs.threshold}
        if [ ${#CTRL[@]} -gt 0 ]; then SWITCH=${CTRL[0]}; fi
        SEACR_1.3.sh "${BG[0]}" "$SWITCH" $ARGS "$P" || exit 1
        BEDS=(*.bed)
        if [ ! -e "${BEDS[0]}" ]; then
            echo "SEACR wrote no peaks" >&2
            exit 1
        fi
        mv "${BEDS[@]}" "$OUT"/
    """,
)

env = flyte.TaskEnvironment.from_task("seacr", seacr_callpeak_cmd.as_task())


async def seacr_callpeak(
    bedgraph: File, control: File | None, threshold: float = 0.05, prefix: str = "sample", args: str = ""
) -> File:
    """Peaks of ``bedgraph`` against ``control``, or its top ``threshold`` without one: ``<prefix>.<mode>.bed``."""
    results = await seacr_callpeak_cmd(
        bedgraph=bedgraph, control=control, threshold=threshold, prefix=prefix, args=args
    )
    beds = [file async for file in results.walk() if file.name.endswith(".bed")]
    if len(beds) != 1:
        raise FileNotFoundError(f"seacr callpeak wrote {len(beds)} BED files, expected 1")
    return beds[0]
