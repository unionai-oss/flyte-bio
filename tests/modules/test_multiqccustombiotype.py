"""Tests for flyte_bio.modules.multiqccustombiotype.

Mirrors upstream custom/multiqccustombiotype, which builds its inputs inline:
the success case matches the snapshot md5s, and a table with more biotypes
than ``--max_biotypes`` (default 100) fails.
"""

import tempfile
from pathlib import Path

from flyte.io import File

from flyte_bio.modules.multiqccustombiotype import multiqc_custom_biotype
from tests.framework import assert_md5, env

HEADER = [
    "# id: 'biotype_counts'",
    "# section_name: 'Biotype Counts'",
    "# description: 'shows reads overlapping genomic features of different biotypes'",
    "# plot_type: 'bargraph'",
    "# anchor: 'featurecounts_biotype'",
]
COLUMNS = ["# Program:featureCounts v2.0.1", "Geneid\tChr\tStart\tEnd\tStrand\tLength\ttest"]


async def write(name: str, lines: list[str]) -> File:
    path = Path(tempfile.mkdtemp()) / name
    path.write_text("\n".join(lines) + "\n")
    return await File.from_local(path)


@env.task
async def test_biotype_counts() -> None:
    # upstream case: custom/multiqccustombiotype "featurecounts biotype counts"
    count = await write(
        "test.featureCounts.txt",
        [
            *COLUMNS,
            "protein_coding\tchr1\t1\t1000\t+\t1000\t5000",
            "rRNA\tchr1\t2000\t3000\t+\t1000\t500",
            "lincRNA\tchr1\t4000\t5000\t+\t1000\t200",
            "miRNA\tchr1\t6000\t7000\t+\t1000\t100",
        ],
    )
    header = await write("biotypes_header.txt", HEADER)
    result = await multiqc_custom_biotype(count, prefix="test", header=header)
    await assert_md5(result.counts_mqc, "6b122075c66e5565fd82804095dfd68a", label="biotype_counts_mqc")
    assert result.rrna_mqc is not None
    await assert_md5(result.rrna_mqc, "ac6ea1053b89872f745a73c4f61b5afe", label="biotype_counts_rrna_mqc")


@env.task
async def test_biotype_counts_too_many_fails() -> None:
    # upstream case: custom/multiqccustombiotype "featurecounts biotype counts - too many biotypes fails"
    rows = [f"biotype_{i}\tchr1\t{i * 1000}\t{i * 1000 + 500}\t+\t500\t{i}" for i in range(1, 151)]
    count = await write("too_many.featureCounts.txt", [*COLUMNS, *rows])
    header = await write("biotypes_header.txt", HEADER)
    try:
        await multiqc_custom_biotype(count, prefix="test", header=header)
    except Exception:
        return
    raise AssertionError("expected 150 biotypes (> --max_biotypes 100) to fail")


tests = [test_biotype_counts, test_biotype_counts_too_many_fails]
