"""Tests for flyte_bio.modules.untarfiles.

The cases and expected values are upstream's untarfiles tests at modules
commit f46e02b, the module's last working version (the snapshot is unchanged
at the pin, 9a1cc9f, where the deprecated module's tests expect failure).
Each runs with upstream's sample id, ``test``, as the prefix: the extracted
files' md5s, and for the flowcell archive the sorted names of its files
(from the archive's version upstream's snapshot was taken with).
"""

import os

from flyte_bio.modules.untarfiles import untarfiles
from tests.framework import assert_md5, env, fixture

# The flowcell archive as upstream's snapshot (2024) saw it: test-datasets added RunParameters.xml to it
# in 2026 (commit 709894c), so this pins the version before.
FLOWCELL = (
    "https://raw.githubusercontent.com/nf-core/test-datasets/68f1b662ef25d42c72907867bacae5573081015f/"
    "data/genomics/homo_sapiens/illumina/bcl/flowcell.tar.gz"
)

# The flowcell archive's files by name, as upstream's sorted snapshot counts them.
FLOWCELL_FILE_COUNTS = {
    "ControlMetricsOut.bin": 1,
    "CorrectedIntMetricsOut.bin": 1,
    "ErrorMetricsOut.bin": 1,
    "ExtractionMetricsOut.bin": 1,
    "IndexMetricsOut.bin": 1,
    "QMetricsOut.bin": 1,
    "RTAConfiguration.xml": 1,
    "RunInfo.xml": 1,
    "TileMetricsOut.bin": 1,
    "config.xml": 1,
    "s_1_1101.bcl": 110,
    "s_1_1101.control": 1,
    "s_1_1101.filter": 1,
    "s_1_1101.locs": 1,
    "s_1_1101.stats": 110,
}


@env.task
async def test_untarfiles_kraken2() -> None:
    # upstream case: untarfiles "test-untarfiles"
    result = await untarfiles(await fixture("genomics/sarscov2/genome/db/kraken2.tar.gz"), prefix="test")
    names = [os.path.basename(f.path) for f in result.files]
    assert names == ["hash.k2d", "opts.k2d", "taxo.k2d"], names
    for file, md5 in zip(
        result.files,
        ["8b8598468f54a7087c203ad0190555d9", "a033d00cf6759407010b21700938f543", "094d5891cdccf2f1468088855c214b2c"],
    ):
        await assert_md5(file, md5, label=os.path.basename(file.path))


@env.task
async def test_untarfiles_subdirectories() -> None:
    # upstream case: untarfiles "test-untarfiles-subdirectories"
    result = await untarfiles(await fixture(FLOWCELL), prefix="test")
    names = sorted(os.path.basename(f.path) for f in result.files)
    expected = sorted(name for name, count in FLOWCELL_FILE_COUNTS.items() for _ in range(count))
    assert names == expected, f"{len(names)} files, expected {len(expected)}: {names[:20]}"


@env.task
async def test_untarfiles_onlyfiles() -> None:
    # upstream case: untarfiles "test-untarfiles-onlyfiles"
    result = await untarfiles(await fixture("generic/tar/hello.tar.gz"), prefix="test")
    names = [os.path.basename(f.path) for f in result.files]
    assert names == ["hello.txt"], names
    await assert_md5(result.files[0], "e59ff97941044f85df5297e1c302d260", label="hello.txt")


tests = [test_untarfiles_kraken2, test_untarfiles_subdirectories, test_untarfiles_onlyfiles]
