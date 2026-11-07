"""macs2 — peak calling for ChIP-seq and ATAC-seq with MACS2.

Exposes :func:`macs2_callpeak`, which calls peaks on treatment reads
(optionally against control reads) and returns MACS2's outputs under their
own names:

- peaks: ``<prefix>_peaks.narrowPeak``, or ``<prefix>_peaks.broadPeak`` with
  ``--broad`` (which also writes ``<prefix>_peaks.gappedPeak``);
- the peak table ``<prefix>_peaks.xls``;
- summits ``<prefix>_summits.bed`` (narrow peaks only);
- with ``-B``/``--bdg``: the pileups ``<prefix>_treat_pileup.bdg`` and
  ``<prefix>_control_lambda.bdg``.

As upstream, the input format is ``BAMPE`` for paired-end reads and ``BAM``
for single-end, unless ``args`` sets ``--format``, which then wins. The reads
are linked into the working directory under their own names, in the caller's
order, so the command line MACS2 records in the ``.xls`` matches upstream's.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

MACS2_IMAGE = "quay.io/biocontainers/macs2:2.2.9.1--py39hff71179_1"
DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory="8Gi")

# A list[File] is staged under the files' names and globbed alphabetically, so
# each list's names ("/"-separated) carry the caller's order into the script.
macs2_callpeak_cmd = shell.create(
    name="macs2_callpeak",
    image=MACS2_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "treatment": list[File],
        "treatment_names": str,
        "control": list[File],
        "control_names": str,
        "gsize": str,
        "format": str,
        "prefix": str,
        "args": str,
    },
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        TREAT=({inputs.treatment})
        shopt -s nullglob; CTRL=({inputs.control}); shopt -u nullglob
        TREAT_NAMES={inputs.treatment_names}
        CTRL_NAMES={inputs.control_names}
        GSIZE={inputs.gsize}
        FORMAT={inputs.format}
        P={inputs.prefix}
        ARGS={inputs.args}
        OUT={outputs.results}
        W=$(mktemp -d); cd "$W"
        mkdir -p .tmp; export TMPDIR="$W/.tmp"

        link_in_order() {
            local DIR=$1 NAMES=$2 NAME
            local -a LIST
            IFS=/ read -ra LIST <<< "$NAMES"
            for NAME in "${LIST[@]}"; do
                ln -s "$DIR/$NAME" "$NAME"
                LINKED+=("$NAME")
            done
        }
        LINKED=(); link_in_order "$(dirname "${TREAT[0]}")" "$TREAT_NAMES"; TREATMENT=("${LINKED[@]}")
        CONTROL=()
        if [ ${#CTRL[@]} -gt 0 ]; then
            LINKED=(); link_in_order "$(dirname "${CTRL[0]}")" "$CTRL_NAMES"; CONTROL=(--control "${LINKED[@]}")
        fi
        INPUTS=("${TREATMENT[@]}" "${CONTROL[@]:1}")

        macs2 \
            callpeak \
            $ARGS \
            --gsize "$GSIZE" \
            --format "$FORMAT" \
            --name "$P" \
            --treatment "${TREATMENT[@]}" \
            "${CONTROL[@]}"
        rm -f "${INPUTS[@]}"  # inputs such as BED reads would otherwise be taken for outputs
        shopt -s nullglob
        mv ./*.narrowPeak ./*.broadPeak ./*.gappedPeak ./*.xls ./*.bed ./*.bdg "$OUT"/
    """,
)

env = flyte.TaskEnvironment.from_task("macs2", macs2_callpeak_cmd.as_task())


def file_name(file: File) -> str:
    return file.path.rstrip("/").rsplit("/", 1)[-1]


def ordered_names(files: list[File], label: str) -> str:
    names = [file_name(f) for f in files]
    if len(set(names)) != len(names):
        raise ValueError(f"macs2 {label} need distinct names, got {names}")
    return "/".join(names)


def split_format(args: str, single_end: bool) -> tuple[str, str]:
    """``args`` without its ``--format <value>``, and the format: that value, else BAM or BAMPE."""
    tokens = args.split()
    fmt = "BAM" if single_end else "BAMPE"
    if "--format" in tokens:
        i = tokens.index("--format")
        if i + 1 >= len(tokens):
            raise ValueError("macs2 args end with --format and no value")
        fmt = tokens[i + 1]
        del tokens[i : i + 2]
    return " ".join(tokens), fmt


@dataclass
class Macs2CallpeakResult:
    peak: File  # <prefix>_peaks.narrowPeak, or <prefix>_peaks.broadPeak with --broad
    xls: File  # <prefix>_peaks.xls
    gapped: File | None  # <prefix>_peaks.gappedPeak, with --broad
    bed: File | None  # <prefix>_summits.bed, for narrow peaks
    bdg: list[File]  # bedGraph pileups, with -B
    results: Dir  # everything MACS2 wrote


async def macs2_callpeak(
    treatment: list[File],
    gsize: str,
    control: list[File] | None = None,
    prefix: str = "sample",
    single_end: bool = False,
    args: str = "",
) -> Macs2CallpeakResult:
    """Call peaks on ``treatment`` reads against optional ``control`` reads.

    ``gsize`` is MACS2's effective genome size: a number, or ``hs``, ``mm``,
    ``ce`` or ``dm``.
    """
    if not treatment:
        raise ValueError("macs2_callpeak needs treatment reads")
    control = control or []
    rest, fmt = split_format(args, single_end)
    results = await macs2_callpeak_cmd(
        treatment=treatment,
        treatment_names=ordered_names(treatment, "treatment files"),
        control=control,
        control_names=ordered_names(control, "control files"),
        gsize=str(gsize),
        format=fmt,
        prefix=prefix,
        args=rest,
    )
    files = {file_name(f): f for f in await results.list_files()}
    peak = files.get(f"{prefix}_peaks.narrowPeak") or files.get(f"{prefix}_peaks.broadPeak")
    xls = files.get(f"{prefix}_peaks.xls")
    if peak is None or xls is None:
        raise FileNotFoundError(f"macs2 wrote no peaks or no peak table, only {sorted(files)}")
    return Macs2CallpeakResult(
        peak=peak,
        xls=xls,
        gapped=files.get(f"{prefix}_peaks.gappedPeak"),
        bed=files.get(f"{prefix}_summits.bed"),
        bdg=[f for name, f in sorted(files.items()) if name.endswith(".bdg")],
        results=results,
    )
