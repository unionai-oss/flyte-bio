"""Tests for flyte_bio.modules.checkqc.

The cases and expected md5s are upstream's checkqc tests (modules commit
9a1cc9f, CheckQC 4.1.0). The "QC OK" case runs with upstream's test config
(``--downgrade-errors ReadsPerSampleHandler``); the "QC Error" case runs
without it. Upstream's run-folder tarballs are passed as they are.
"""

import json

from flyte.io import File

from flyte_bio.modules.checkqc import checkqc
from tests.framework import assert_md5, env, fixture

FLOWCELL = "genomics/homo_sapiens/illumina/interop_bcl2fastqstats/test_flowcell_stats.tar.gz"
MISEQ = (
    "https://github.com/nf-core/test-datasets/raw/demultiplex/testdata/MiSeq/220422_M11111_0222_000000000-K9H97.tar.gz"
)


async def read_json(file: File) -> dict:
    async with file.open("rb") as fh:
        return json.loads(bytes(await fh.read()).decode())


@env.task
async def test_checkqc_qc_ok() -> None:
    # upstream case: "homo_sapiens illumina test_interop_bcl2fastqstats_flowcell - QC OK"
    report = await checkqc(await fixture(FLOWCELL), args="--downgrade-errors ReadsPerSampleHandler")
    await assert_md5(report, "c684032b1172084ee1adbbf8461b9834", label="checkqc_report.json")
    data = await read_json(report)
    assert data["exit_status"] == 0, data["exit_status"]
    assert data["run_summary"]["instrument_and_reagent_type"] == "hiseq2500_rapidhighoutput_v4"


@env.task
async def test_checkqc_missing_run_parameters() -> None:
    # upstream case: "homo_sapiens illumina test_miseq_noRunParameters - Log Error: no File"
    try:
        await checkqc(await fixture(MISEQ))
    except Exception as e:  # the task fails: CheckQC writes no report
        print(f"checkqc failed as expected: {e}")
    else:
        raise AssertionError("checkqc succeeded on a run folder without RunParameters.xml")


@env.task
async def test_checkqc_qc_error() -> None:
    # upstream case: "homo_sapiens illumina test_interop_bcl2fastqstats_flowcell - QC Error"
    report = await checkqc(await fixture(FLOWCELL))
    await assert_md5(report, "ffb2efadbcfe9dfb8a3c665e80519073", label="checkqc_report.json")
    data = await read_json(report)
    assert data["exit_status"] == 1, data["exit_status"]
    handlers = {h["name"]: h for h in data["run_summary"]["handlers"]}
    assert handlers["UndeterminedPercentageHandler"]["error"] == 9, handlers["UndeterminedPercentageHandler"]


tests = [test_checkqc_qc_ok, test_checkqc_missing_run_parameters, test_checkqc_qc_error]
