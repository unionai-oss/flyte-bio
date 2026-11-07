"""multiqcsav — MultiQC report with Illumina sequencing run (SAV) metrics.

Exposes :func:`multiqcsav`, which runs MultiQC with the MultiQC-SAV plugin
over an Illumina run folder's ``RunInfo.xml`` and ``InterOp/*.bin`` files,
plus any other QC files (e.g. FastQC zips), writing
``<prefix>.html`` (default ``multiqc_report.html``), ``*_data/`` and, with
``-p`` in ``args``, ``*_plots/`` into an output Dir.

Run folders hold hundreds of files, so the folder comes in as a tar archive
and the task picks out what upstream's module stages: the first
``RunInfo.xml`` and every ``.bin`` file directly under an ``InterOp``
directory. As upstream, each extra file and config file is staged in its own
numbered directory (``1/``, ``2/``, ...), and configs apply in list order.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.modules.multiqc import ordered

MULTIQCSAV_IMAGE = "community.wave.seqera.io/library/multiqc_multiqc_sav_pip_interop:9b10d606ce2f36b6"

DEFAULT_RESOURCES = flyte.Resources(cpu=2, memory="8Gi")

multiqcsav_cmd = shell.create(
    name="multiqcsav",
    image=MULTIQCSAV_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "run_folder": File,
        "extra": list[File],
        "configs": list[File],
        "logo": File | None,
        "replace_names": File | None,
        "sample_names": File | None,
        "prefix": str,
        "args": str,
    },
    defaults={"extra": [], "configs": [], "prefix": "", "args": ""},
    outputs={"results": Dir},
    script=r"""
        RUN=({inputs.run_folder})
        shopt -s nullglob
        EXTRA=({inputs.extra}); CONFIGS=({inputs.configs})
        LOGO=({inputs.logo}); REPLACE=({inputs.replace_names}); SAMPLES=({inputs.sample_names})
        shopt -u nullglob
        P={inputs.prefix}
        ARGS={inputs.args}
        X=$(mktemp -d)
        tar -C "$X" -xf "${RUN[0]}" || exit 1
        W=$(mktemp -d); cd "$W"
        XML=$(find "$X" -name RunInfo.xml | sort | head -n 1)
        if [ -z "$XML" ]; then echo "no RunInfo.xml in the run folder" >&2; exit 1; fi
        ln -s "$XML" RunInfo.xml
        mkdir InterOp
        while IFS= read -r BIN; do
            if [ "$(basename "$(dirname "$BIN")")" = InterOp ]; then ln -s "$BIN" InterOp/; fi
        done < <(find "$X" -type f -name '*.bin' | sort)
        OPTS=()
        i=1
        for f in "${EXTRA[@]}"; do mkdir "$i"; ln -s "$f" "$i"/; i=$((i + 1)); done
        for c in "${CONFIGS[@]}"; do
            mkdir "$i"; ln -s "$c" "$i"/; OPTS+=(--config "$i/$(basename "$c")"); i=$((i + 1))
        done
        if [ -n "$P" ]; then OPTS+=(--filename "$P.html"); fi
        if [ ${#LOGO[@]} -gt 0 ]; then
            ln -s "${LOGO[0]}" .
            OPTS+=(--cl-config "custom_logo: \"$(basename "${LOGO[0]}")\"")
        fi
        if [ ${#REPLACE[@]} -gt 0 ]; then OPTS+=(--replace-names "${REPLACE[0]}"); fi
        if [ ${#SAMPLES[@]} -gt 0 ]; then OPTS+=(--sample-names "${SAMPLES[0]}"); fi
        multiqc --force $ARGS "${OPTS[@]}" . || exit 1
        shopt -s nullglob
        mv ./*.html ./*_data ./*_plots {outputs.results}/
    """,
)

env = flyte.TaskEnvironment.from_task("multiqcsav", multiqcsav_cmd.as_task())


@dataclass
class MultiQCSavResult:
    report: File  # <prefix>.html
    results: Dir  # the report, *_data/ and (with -p) *_plots/


async def multiqcsav(
    run_folder: File,
    extra: list[File] | None = None,
    configs: list[File] | None = None,
    logo: File | None = None,
    replace_names: File | None = None,
    sample_names: File | None = None,
    prefix: str = "",
    args: str = "",
) -> MultiQCSavResult:
    """MultiQC with SAV metrics for the run folder archived in ``run_folder``, plus the ``extra`` QC files."""
    results = await multiqcsav_cmd(
        run_folder=run_folder,
        extra=extra or [],
        configs=await ordered(configs or []),
        logo=logo,
        replace_names=replace_names,
        sample_names=sample_names,
        prefix=prefix,
        args=args,
    )
    name = f"{prefix}.html" if prefix else "multiqc_report.html"
    report = await results.get_file(name)
    if report is None:
        raise FileNotFoundError(f"multiqc wrote no {name}")
    return MultiQCSavResult(report=report, results=results)
