"""Tests for flyte_bio.modules.repeatmasker.

The cases and md5s are upstream's repeatmasker/repeatmasker tests (modules
commit 9a1cc9f, RepeatMasker 4.1.5), with their test config's ``args``
(``-no_is -gff -engine hmmer``) and HMM library. As upstream, the summary
table (which records the run's date) is checked by content.
"""

from flyte.io import File

from flyte_bio.modules.repeatmasker import RepeatmaskerResult, repeatmasker
from tests.framework import assert_md5, env, fixture

GENOME = "genomics/homo_sapiens/genome/"
ARGS = "-no_is -gff -engine hmmer"


async def check(result: RepeatmaskerResult) -> None:
    await assert_md5(result.masked, "788f27945020adc342db9d56ff29a840", label="test.masked")
    await assert_md5(result.out, "71ec2969baaa8ea492808c3ae6e7012e", label="test.out")
    assert result.gff is not None, "RepeatMasker wrote no GFF"
    await assert_md5(result.gff, "5cd19732692c496fed8e7ad3c2f52f17", label="test.gff")
    assert "run with nhmmscan" in await text(result.tbl), "summary table lacks 'run with nhmmscan'"


async def text(file: File) -> str:
    async with file.open("rb") as fh:
        return bytes(await fh.read()).decode()


@env.task
async def test_repeatmasker_fasta_hmm_lib() -> None:
    # upstream case: repeatmasker/repeatmasker "homo_sapiens - genome - fasta - hmm lib"
    lib = await fixture(GENOME + "test_starfusion_dfam.hmm")
    await check(await repeatmasker(await fixture(GENOME + "genome.fasta"), lib, prefix="test", args=ARGS))


@env.task
async def test_repeatmasker_fasta_gz_hmm_lib() -> None:
    # upstream case: repeatmasker/repeatmasker "homo_sapiens - genome - fasta.gz - hmm lib"
    lib = await fixture(GENOME + "test_starfusion_dfam.hmm")
    await check(await repeatmasker(await fixture(GENOME + "genome.fasta.gz"), lib, prefix="test", args=ARGS))


tests = [test_repeatmasker_fasta_hmm_lib, test_repeatmasker_fasta_gz_hmm_lib]
