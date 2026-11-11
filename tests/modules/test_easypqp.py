"""Tests for flyte_bio.modules.easypqp.

The cases and expected values are upstream's easypqp/convert and
easypqp/library tests (modules commit 9a1cc9f, easypqp 0.1.60), with the
library test config's ``ext.args``. As upstream, the PSM pickle isn't
byte-stable between runs, so only its name and size are checked; the peak
pickle is md5ed, and the library's first 5 lines are compared.
"""

from flyte.io import File

from flyte_bio.modules.easypqp import EasypqpConvertResult, easypqp_convert, easypqp_library
from tests.framework import assert_md5, assert_nonempty, env, fixture

IDXML = "proteomics/openms/HepG2_rep1_small.idXML"
MZML = "proteomics/msspectra/HepG2_rep1_small.mzML"

LIBRARY_HEAD = [
    (
        "PrecursorMz\tProductMz\tAnnotation\tProteinId\tGeneName\tPeptideSequence\tModifiedPeptideSequence"
        "\tPrecursorCharge\tLibraryIntensity\tNormalizedRetentionTime\tPrecursorIonMobility"
    ),
    *(
        f"271.410568\t{mz}\t{ion}\tsp|Q96EN8|MOCOS_HUMAN\t-\tKAAGVLEGALGP\tKAAGVLEGALGP\t4\t{intensity}"
        "\t93.0079411195042\t"
        for mz, ion, intensity in [
            ("110.070828", "b4^3", "10000.0"),
            ("116.070606", "y1^1", "398.89114"),
            ("129.10224", "b1^1", "3139.388"),
            ("271.176468", "b3^1", "6836.362"),
        ]
    ),
]


def easypqp_name(file: File) -> str:
    return file.path.rsplit("/", 1)[-1]


async def easypqp_test_convert() -> EasypqpConvertResult:
    return await easypqp_convert(await fixture(IDXML), await fixture(MZML))


@env.task
async def test_easypqp_convert() -> None:
    # upstream case: easypqp/convert "proteomics - idxml - mzml"
    result = await easypqp_test_convert()
    assert easypqp_name(result.psmpkl) == "HepG2_rep1_small.psmpkl", result.psmpkl.path
    await assert_nonempty(result.psmpkl, label="psmpkl")
    assert easypqp_name(result.peakpkl) == "HepG2_rep1_small.peakpkl", result.peakpkl.path
    await assert_md5(result.peakpkl, "0785c76c6d5e085831ce6d136ce9c0ba", label="peakpkl")


@env.task
async def test_easypqp_library() -> None:
    # upstream case: easypqp/library "proteomics - psmpkl - peakpkl"
    converted = await easypqp_test_convert()
    result = await easypqp_library(
        [converted.psmpkl],
        [converted.peakpkl],
        prefix="test",
        args="--perform_rt_calibration False --perform_im_calibration False --nofdr",
    )
    async with result.tsv.open("rb") as fh:
        lines = bytes(await fh.read()).decode().splitlines()[:5]
    assert lines == LIBRARY_HEAD, f"library lines 1-5 differ: {lines}"
    names = [easypqp_name(f) for f in result.run_peaks]
    assert names == ["HepG2_rep1_small_run_peaks.tsv"], names


tests = [test_easypqp_convert, test_easypqp_library]
