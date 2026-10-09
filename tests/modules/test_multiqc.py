"""Tests for flyte_bio.modules.multiqc.

Mirrors upstream multiqc "sarscov2 single-end [fastqc]": the report and data
directory are produced (their content embeds the run date, so no md5) and
the FastQC sample was picked up.
"""

import tempfile
from pathlib import Path

from flyte.io import Dir

from flyte_bio.modules.multiqc import multiqc
from tests.framework import env, fixture


@env.task
async def test_multiqc_fastqc() -> None:
    # upstream case: multiqc "sarscov2 single-end [fastqc]"
    zipped = await fixture("genomics/sarscov2/illumina/fastqc/test_fastqc.zip")
    root = Path(tempfile.mkdtemp())
    await zipped.download(str(root / "test_fastqc.zip"))
    result = await multiqc(await Dir.from_local(root))
    stats = None
    async for f in result.results.walk():
        if f.path.endswith("multiqc_report_data/multiqc_general_stats.txt"):
            stats = f
    assert stats is not None, "no multiqc_report_data/multiqc_general_stats.txt"
    async with stats.open("rb") as fh:
        text = bytes(await fh.read()).decode()
    assert "test" in text, text[:300]


tests = [test_multiqc_fastqc]
