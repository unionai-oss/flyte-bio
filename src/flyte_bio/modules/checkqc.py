"""checkqc — check an Illumina run folder against QC thresholds.

Exposes :func:`checkqc`, which runs CheckQC on a run folder (its
``RunInfo.xml``, ``RunParameters.xml`` and bcl2fastq/bclconvert statistics)
and returns its JSON report, ``checkqc_report.json``.

Run folders hold hundreds of files, so the folder comes in as a tar archive
(``<run>.tar.gz``) and is extracted inside the task into ``<run>/``. As
upstream, a single top-level directory in the archive is stripped. CheckQC's
own exit status (1 on fatal QC errors) is recorded in the report's
``exit_status``; the task fails only when CheckQC writes no report, e.g. when
the run folder lacks ``RunParameters.xml``.
"""

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

CHECKQC_IMAGE = "community.wave.seqera.io/library/sample-sheet_numpy_pandas_pip_pruned:0b9dc0869e46a949"

# Upstream's process_single: one CPU.
DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory="4Gi")

checkqc_cmd = shell.create(
    name="checkqc",
    image=CHECKQC_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"run_folder": File, "config": File | None, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        RUN=({inputs.run_folder})
        shopt -s nullglob; CONFIG=({inputs.config}); shopt -u nullglob
        ARGS={inputs.args}
        OPTS=()
        if [ ${#CONFIG[@]} -gt 0 ]; then OPTS+=(--config "${CONFIG[0]}"); fi
        NAME=$(basename "${RUN[0]}")
        NAME=${NAME%.gz}; NAME=${NAME%.tgz}; NAME=${NAME%.tar}
        W=$(mktemp -d); cd "$W"
        mkdir -p "$NAME"
        # As upstream: strip the top level only when every entry sits under one directory.
        tar -tf "${RUN[0]}" | sed 's|^\./||' | grep -v '^$' > entries.txt || exit 1
        TOPS=$(cut -d/ -f1 entries.txt | sort -u | wc -l)
        NESTED=$(grep -c / entries.txt || true)
        if [ "$TOPS" -eq 1 ] && [ "$NESTED" -gt 0 ]; then
            tar -C "$NAME" --strip-components 1 -xf "${RUN[0]}" || exit 1
        else
            tar -C "$NAME" -xf "${RUN[0]}" || exit 1
        fi
        checkqc $ARGS "${OPTS[@]}" --json "$NAME" > {outputs.results}/checkqc_report.json || true
        if [ ! -s {outputs.results}/checkqc_report.json ]; then
            echo "Error: CheckQC wrote no JSON report; most likely files are missing from the run folder." >&2
            exit 1
        fi
    """,
)

env = flyte.TaskEnvironment.from_task("checkqc", checkqc_cmd.as_task())


async def checkqc(run_folder: File, config: File | None = None, args: str = "") -> File:
    """CheckQC's ``checkqc_report.json`` for the run folder archived in ``run_folder`` (a ``.tar.gz``)."""
    results = await checkqc_cmd(run_folder=run_folder, config=config, args=args)
    report = await results.get_file("checkqc_report.json")
    if report is None:
        raise FileNotFoundError("checkqc wrote no checkqc_report.json")
    return report
