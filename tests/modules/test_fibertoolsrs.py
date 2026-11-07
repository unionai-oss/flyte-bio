"""Tests for flyte_bio.modules.fibertoolsrs.

The cases and expected md5s are upstream's fibertoolsrs/predictm6a,
fibertoolsrs/extract and fibertoolsrs/addnucleosomes tests (modules commit
9a1cc9f, fibertools-rs 0.8.2), on the module's 2 CPUs: the output BAMs'
``@PG`` command lines record ``--threads``, and upstream's snapshots are of
2-thread runs. extract's input is predict-m6a's output, as upstream's test
setup. nf-test hashes the BAMs as they are and the BEDs decompressed.
"""

from flyte_bio.modules.fibertoolsrs import fibertoolsrs_addnucleosomes, fibertoolsrs_extract, fibertoolsrs_predictm6a
from tests.framework import assert_gunzipped_md5, assert_md5, env, fixture

PACBIO_BAM = "genomics/homo_sapiens/pacbio/bam/NA03697B2_downsampled.pbmm2.repeats.bam"
NANOPORE_BAM = "genomics/homo_sapiens/nanopore/bam/test.sorted.bam"


@env.task
async def test_fibertoolsrs_predictm6a() -> None:
    # upstream case: fibertoolsrs/predictm6a "homo_sapiens - bam"
    bam = await fibertoolsrs_predictm6a(await fixture(PACBIO_BAM), prefix="test_m6a")
    assert bam.name == "test_m6a.bam", bam.name
    await assert_md5(bam, "e5054c244fa78c7b504ae0d2ca7e57ca", label="test_m6a.bam")


@env.task
async def test_fibertoolsrs_extract() -> None:
    # upstream cases: fibertoolsrs/extract "homo_sapiens - bam" (all), "... - cpg" and "... - invalid_extract_type"
    m6a = await fibertoolsrs_predictm6a(await fixture(PACBIO_BAM), prefix="test_m6a")
    for extract_type, md5 in [("all", "3bdde73d10c5369f509cfac1e59b5f68"), ("cpg", "ed4f4ea245c88cb2c40c0866eb779fb0")]:
        bed = await fibertoolsrs_extract(m6a, extract_type, prefix="test")
        assert bed.name == "test.bed.gz", bed.name
        await assert_gunzipped_md5(bed, md5, label=f"extract --{extract_type}")
    try:
        await fibertoolsrs_extract(m6a, "m5c", prefix="test")
    except ValueError as e:
        assert "Invalid extract_type" in str(e), e
    else:
        raise AssertionError("fibertoolsrs_extract accepted extract_type 'm5c'")


@env.task
async def test_fibertoolsrs_addnucleosomes() -> None:
    # upstream case: fibertoolsrs/addnucleosomes "homo_sapiens - bam"
    bam = await fibertoolsrs_addnucleosomes(await fixture(NANOPORE_BAM), prefix="test_nucleosomes")
    assert bam.name == "test_nucleosomes.bam", bam.name
    await assert_md5(bam, "895e6f7118274676037e3b6609d52752", label="test_nucleosomes.bam")


tests = [test_fibertoolsrs_predictm6a, test_fibertoolsrs_extract, test_fibertoolsrs_addnucleosomes]
