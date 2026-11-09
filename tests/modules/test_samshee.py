"""Tests for flyte_bio.modules.samshee.

The cases and expected md5s are upstream's samshee tests (modules commit
9a1cc9f, samshee 0.2.13), with their configs' schemas and options: a v1
MiSeq sheet checked against an inline schema and a schema file and written
sectioned, and a v2 NextSeq 2000 sheet checked against samshee's Illumina v2
schema.
"""

import os
import tempfile
from pathlib import Path

from flyte.io import File

from flyte_bio.modules.samshee import samshee
from tests.framework import assert_md5, env, fixture

TESTDATA = "https://raw.githubusercontent.com/nf-core/test-datasets/demultiplex/testdata/"
# The schema file upstream's v1 test writes.
SCHEMA_FILE = """{
                    "$schema": "https://json-schema.org/draft/2020-12/schema",
                    "required": ["Settings"]
                }"""


async def samshee_schema_file() -> File:
    path = Path(tempfile.mkdtemp()) / "schema.json"
    path.write_text(SCHEMA_FILE)
    return await File.from_local(str(path))


@env.task
async def test_samshee_samplesheet_v1() -> None:
    # upstream case: "test samplesheet_v1"
    sheet = await fixture(TESTDATA + "miseq_35147139/miseq_35147139_samplesheet.csv")
    formatted = await samshee(
        sheet, await samshee_schema_file(), schema='{"required": ["Data"]}', args="--output-format sectioned"
    )
    assert os.path.basename(formatted.path) == "miseq_35147139_samplesheet_formatted.csv", formatted.path
    await assert_md5(formatted, "2a844a232fcf9cd7131185a55c2602b0", label="v1 formatted sheet")


@env.task
async def test_samshee_samplesheet_v2() -> None:
    # upstream case: "test samplesheet_v2"
    sheet = await fixture(TESTDATA + "NextSeq2000/SampleSheet.csv")
    formatted = await samshee(sheet, schema='{"$ref": "urn:samshee:illuminav2/v1"}')
    assert os.path.basename(formatted.path) == "SampleSheet_formatted.csv", formatted.path
    await assert_md5(formatted, "779c2dda1994834a5795462d76a5c14b", label="v2 formatted sheet")


tests = [test_samshee_samplesheet_v1, test_samshee_samplesheet_v2]
