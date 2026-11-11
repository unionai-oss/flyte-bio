"""Tests for flyte_bio.modules.fcs.

The cases and md5s are upstream's fcs/fcsadaptor tests (modules commit
9a1cc9f, FCS-adaptor 0.5.0). As upstream, the logs (and, in the adapter case,
the report) are checked by content, not md5. The adapter case's input is the
SARS-CoV-2 genome with two PacBio adapter sequences inserted, built the way
upstream's test builds it.
"""

import tempfile
from pathlib import Path

from flyte.io import File

from flyte_bio.modules.fcs import fcs_fcsadaptor
from tests.framework import assert_gunzipped_md5, assert_md5, env, fixture

LOG_SUCCESS = "[step all_cleaned_fasta] completed success"
PIPELINE_ARGS_MD5 = "4a748f51e82fa0639cb06242bc3df318"


async def text(file: File) -> str:
    async with file.open("rb") as fh:
        return bytes(await fh.read()).decode()


async def genome_with_adapters() -> File:
    """Upstream's ``sarscov2_with_adapters.fasta``: an adapter mid-genome and one at the end."""
    genome = await text(await fixture("genomics/sarscov2/genome/genome.fasta"))
    lines = genome.split("\n")
    while lines and lines[-1] == "":  # as Java's split, which drops trailing empty strings
        lines.pop()
    half = len(lines) // 2
    combined = [
        *lines[:half],
        "ATCTCTCTCTTTTCCTCCTCCTCCGTTGTTGTTGTTGAGAGAGAT",
        *lines[half:],
        "AAAAAAAAAAAAAAAAAATTAACGGAGGAGGAGGA",
    ]
    path = Path(tempfile.mkdtemp()) / "sarscov2_with_adapters.fasta"
    path.write_text("\n".join(combined))
    return await File.from_local(str(path))


@env.task
async def test_fcs_fcsadaptor_bacteroides() -> None:
    # upstream case: fcs/fcsadaptor "bacteroides_fragilis-genome_fna_gz"
    genome = await fixture("genomics/prokaryotes/bacteroides_fragilis/genome/genome.fna.gz")
    result = await fcs_fcsadaptor(genome, prefix="test")
    await assert_md5(result.adaptor_report, "27375be4671e01d2d2674ddd1d44414a", label="adaptor report")
    await assert_md5(result.pipeline_args, PIPELINE_ARGS_MD5, label="pipeline_args.yaml")
    assert result.cleaned_assembly is None, "a clean genome gave a cleaned assembly"
    assert LOG_SUCCESS in await text(result.log), "FCS-adaptor log lacks its success line"


@env.task
async def test_fcs_fcsadaptor_pacbio_adapters() -> None:
    # upstream case: fcs/fcsadaptor "sarscov2-genome_fasta-pacbio_adapters"
    result = await fcs_fcsadaptor(await genome_with_adapters(), prefix="test")
    assert result.cleaned_assembly is not None, "no cleaned assembly despite adapters"
    await assert_gunzipped_md5(result.cleaned_assembly, "9a922a4dbd72360cdf97fedf4538983b", label="cleaned")
    await assert_md5(result.pipeline_args, PIPELINE_ARGS_MD5, label="pipeline_args.yaml")
    await assert_md5(result.skipped_trims, "cbcdf071a7f61519e34d6b60c1aa9f51", label="skipped_trims.jsonl")
    assert LOG_SUCCESS in await text(result.log), "FCS-adaptor log lacks its success line"
    report = await text(result.adaptor_report)
    for adapter in ("Pacific Biosciences Blunt Adapter", "Pacific Biosciences C2 Primer"):
        assert adapter in report, f"report lacks {adapter!r}"


tests = [test_fcs_fcsadaptor_bacteroides, test_fcs_fcsadaptor_pacbio_adapters]
