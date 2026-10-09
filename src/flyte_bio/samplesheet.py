"""samplesheet — the CSV a pipeline takes as its input, as upstream's ``--input``.

Each pipeline defines its own columns (rnaseq: ``sample, fastq_1, fastq_2,
strandedness``) and reads them straight off the rows: a samplesheet missing
one fails the run with an error naming the column. File paths in the sheet
must be URIs the cluster can read (``s3://``, ``gs://``, ``https://`` …).
"""

import csv
import io

from flyte.io import File


class Row(dict[str, str]):
    """One samplesheet row: header -> value (whitespace-stripped)."""

    def __init__(self, values: dict[str, str], line: int) -> None:
        super().__init__(values)
        self.line = line  # 1-based line number in the CSV, for error messages

    def __missing__(self, column: str) -> str:
        raise KeyError(f"samplesheet has no column {column!r} (line {self.line} has {sorted(self)})")


async def read_samplesheet(samplesheet: File) -> list[Row]:
    """The rows of a samplesheet CSV, in file order.

    The samplesheet is read in the calling task, so it should be in object
    storage or a local file passed to ``flyte run`` (which uploads it):
    reading ``https://`` files inside a task isn't reliable.
    """
    async with samplesheet.open("rb") as fh:
        text = bytes(await fh.read()).decode()
    reader = csv.DictReader(io.StringIO(text))
    return [
        Row({k.strip(): (v or "").strip() for k, v in row.items() if k is not None}, line)
        for line, row in enumerate(reader, start=2)
    ]
