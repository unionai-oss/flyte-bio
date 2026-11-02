"""qcatch — cell-level quality control of an alevin-fry (simpleaf) quantification.

Exposes :func:`qcatch`, which reads a simpleaf ``af_quant`` directory, calls
cells and writes ``<prefix>_qcatch_report.html``, the matrix of the called
cells as ``<prefix>_filtered_quants.h5ad`` and ``<prefix>_metrics_summary.csv``.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

QCATCH_IMAGE = "community.wave.seqera.io/library/pip_qcatch:03b88593a5cca75b"

DEFAULT_RESOURCES = flyte.Resources(cpu=2, memory="8Gi")

# matplotlib and numba cache under $HOME by default, which may not be writable.
qcatch_cmd = shell.create(
    name="qcatch",
    image=QCATCH_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"quant": Dir, "chemistry": str, "prefix": str, "args": str},
    defaults={"chemistry": "", "args": ""},
    outputs={"results": Dir},
    script=r"""
        QUANT={inputs.quant}
        CHEM={inputs.chemistry}
        ARGS={inputs.args}
        W=$(mktemp -d); cd "$W"
        mkdir tmp
        export MPLCONFIGDIR="$W/tmp" NUMBA_CACHE_DIR="$W/tmp" NUMBA_DISABLE_CACHE=1
        CHEM_ARG=(); [ -n "$CHEM" ] && CHEM_ARG=(--chemistry "$CHEM")
        qcatch --input "$QUANT" --output {inputs.prefix} "${CHEM_ARG[@]}" \
            --save_filtered_h5ad --export_summary_table $ARGS
        mv {inputs.prefix}/QCatch_report.html {outputs.results}/{inputs.prefix}_qcatch_report.html
        mv {inputs.prefix}/filtered_quants.h5ad {outputs.results}/{inputs.prefix}_filtered_quants.h5ad
        mv {inputs.prefix}/summary_table.csv {outputs.results}/{inputs.prefix}_metrics_summary.csv
    """,
)


env = flyte.TaskEnvironment.from_task(
    "qcatch",
    qcatch_cmd.as_task(),
)


@dataclass
class QcatchResult:
    report: File  # <prefix>_qcatch_report.html
    filtered_h5ad: File  # <prefix>_filtered_quants.h5ad
    metrics_summary: File  # <prefix>_metrics_summary.csv


async def qcatch(
    quant: Dir,
    prefix: str,
    chemistry: str | None = None,
    remove_doublets: bool = False,
    n_partitions: int | None = None,
    args: str = "",
) -> QcatchResult:
    """QCatch over a simpleaf ``af_quant`` directory.

    ``chemistry`` is QCatch's (e.g. ``10X_3p_v3``); without it QCatch infers
    it from the quantification's metadata, or uses ``n_partitions``.
    ``remove_doublets`` also removes and plots doublets, as upstream.
    """
    extra = [args] if args else []
    if remove_doublets:
        extra.append("--remove_doublets --visualize_doublets")
    if n_partitions is not None:
        extra.append(f"--n_partitions {n_partitions}")
    results = await qcatch_cmd(quant=quant, chemistry=chemistry or "", prefix=prefix, args=" ".join(extra))

    async def pick(suffix: str) -> File:
        f = await results.get_file(f"{prefix}{suffix}")
        if f is None:
            raise FileNotFoundError(f"qcatch wrote no {prefix}{suffix}")
        return f

    return QcatchResult(
        report=await pick("_qcatch_report.html"),
        filtered_h5ad=await pick("_filtered_quants.h5ad"),
        metrics_summary=await pick("_metrics_summary.csv"),
    )
