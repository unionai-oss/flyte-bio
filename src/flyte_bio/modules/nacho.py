"""nacho — NanoString nCounter QC and normalisation with NACHO (R).

Exposes two functions, each running one of upstream's NACHO R scripts
(vendored under :mod:`flyte_bio.scripts`) on a set of RCC files and a sample
sheet. The sample sheet is a CSV whose ``RCC_FILE_NAME`` column names the RCC
files; as upstream, the RCC files are staged into ``input/``.

- :func:`nacho_qc`: QC reports and MultiQC custom content:
  ``NanoQC.html`` and ``NanoQC_with_outliers.html``; 16 plots
  ``*_mqc.png``; ``hk_detected_mqc.txt`` (the housekeeping genes detected)
  and ``normalized_qc_mqc.txt`` (per-sample QC metrics, PASS/FAIL).
- :func:`nacho_normalize`: normalised counts, ``<prefix>.tsv``, and the same
  without housekeeping-gene normalisation, ``<prefix>_wo_HKnorm.tsv``.
  ``args`` takes the script's options: ``--norm_method GLM`` (default) or
  ``GEO``.

Each also writes R's session info and a ``versions.yml``.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.scripts import path

NACHO_IMAGE = "community.wave.seqera.io/library/r-base_r-dplyr_r-fs_r-ggplot2_pruned:2bcf3b351adffd77"
DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory="4Gi")  # upstream: process_single

# The RCC files go in ./input and the sample sheet in the working directory
# under its own name, as upstream stages them.
STAGE_INPUTS = r"""
        SCRIPT=({inputs.script}); RCC=({inputs.rcc_files}); SHEET=({inputs.sample_sheet})
        OUT={outputs.results}
        W=$(mktemp -d); cd "$W"
        export HOME="$W"
        mkdir input
        for F in "${RCC[@]}"; do ln -s "$F" input/"$(basename "$F")"; done
        SHEET_NAME=$(basename "${SHEET[0]}")
        ln -s "${SHEET[0]}" "$SHEET_NAME"
"""

nacho_qc_cmd = shell.create(
    name="nacho_qc",
    image=NACHO_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"script": File, "rcc_files": list[File], "sample_sheet": File},
    outputs={"results": Dir},
    script=STAGE_INPUTS
    + r"""
        Rscript "${SCRIPT[0]}" "$SHEET_NAME" NACHO_QC || exit 1
        shopt -s nullglob
        mv ./*.html ./*_mqc.png ./*_mqc.txt ./R_sessionInfo.log ./versions.yml "$OUT"/
    """,
)

nacho_normalize_cmd = shell.create(
    name="nacho_normalize",
    image=NACHO_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"script": File, "rcc_files": list[File], "sample_sheet": File, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=STAGE_INPUTS
    + r"""
        P={inputs.prefix}
        Rscript "${SCRIPT[0]}" "$P" "$SHEET_NAME" {inputs.args} NACHO_NORMALIZE || exit 1
        shopt -s nullglob
        mv "$P.tsv" "${P}_wo_HKnorm.tsv" ./*.R_sessionInfo.log ./versions.yml "$OUT"/
    """,
)

env = flyte.TaskEnvironment.from_task("nacho", nacho_qc_cmd.as_task(), nacho_normalize_cmd.as_task())


@dataclass
class NachoQcResult:
    reports: list[File]  # NanoQC.html, NanoQC_with_outliers.html
    plots: list[File]  # *_mqc.png
    tables: list[File]  # hk_detected_mqc.txt, normalized_qc_mqc.txt
    results: Dir  # everything above, the R session info and versions.yml


@dataclass
class NachoNormalizeResult:
    normalized_counts: File  # <prefix>.tsv
    normalized_counts_wo_hk: File  # <prefix>_wo_HKnorm.tsv
    results: Dir  # both, the R session info and versions.yml


async def nacho_qc(rcc_files: list[File], sample_sheet: File) -> NachoQcResult:
    """NACHO QC reports and MultiQC content for ``rcc_files``, described by ``sample_sheet``."""
    script = await File.from_local(str(path("nacho_qc.R")))
    results = await nacho_qc_cmd(script=script, rcc_files=rcc_files, sample_sheet=sample_sheet)
    files = sorted(await results.list_files(), key=lambda f: f.name)
    reports = [f for f in files if f.name.endswith(".html")]
    if not reports:
        raise FileNotFoundError("nacho_qc wrote no HTML reports")
    return NachoQcResult(
        reports=reports,
        plots=[f for f in files if f.name.endswith("_mqc.png")],
        tables=[f for f in files if f.name.endswith("_mqc.txt")],
        results=results,
    )


async def nacho_normalize(
    rcc_files: list[File], sample_sheet: File, prefix: str = "sample", args: str = ""
) -> NachoNormalizeResult:
    """NACHO-normalised counts of ``rcc_files``; ``args`` e.g. ``--norm_method GEO``."""
    script = await File.from_local(str(path("nacho_norm.R")))
    results = await nacho_normalize_cmd(
        script=script, rcc_files=rcc_files, sample_sheet=sample_sheet, prefix=prefix, args=args
    )
    counts = await results.get_file(f"{prefix}.tsv")
    counts_wo_hk = await results.get_file(f"{prefix}_wo_HKnorm.tsv")
    if counts is None or counts_wo_hk is None:
        raise FileNotFoundError(f"nacho_normalize wrote no {prefix}.tsv or {prefix}_wo_HKnorm.tsv")
    return NachoNormalizeResult(normalized_counts=counts, normalized_counts_wo_hk=counts_wo_hk, results=results)
