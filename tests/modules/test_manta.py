"""Tests for flyte_bio.modules.manta.

The cases are upstream's manta/somatic tests (modules commit 9a1cc9f, Manta
1.6.0): a CRAM tumour/normal pair with a configuration file, and the same
pair within a BED's regions. As upstream, each VCF is checked for its
``##fileformat=VCFv4.1`` header line.
"""

import asyncio
import gzip
import tempfile
from pathlib import Path

from flyte.io import File

from flyte_bio.modules.manta import MantaSomaticResult, manta_somatic
from tests.framework import env, fixture

HS = "genomics/homo_sapiens/"
CRAM = HS + "illumina/cram/{}.paired_end.recalibrated.sorted.cram"
CHR21 = HS + "genome/chr21/sequence/"


async def tumour_normal_inputs() -> dict[str, File]:
    normal, normal_index, tumor, tumor_index, fasta, fai = await asyncio.gather(
        fixture(CRAM.format("test")),
        fixture(CRAM.format("test") + ".crai"),
        fixture(CRAM.format("test2")),
        fixture(CRAM.format("test2") + ".crai"),
        fixture(CHR21 + "genome.fasta"),
        fixture(CHR21 + "genome.fasta.fai"),
    )
    return {
        "normal": normal,
        "normal_index": normal_index,
        "tumor": tumor,
        "tumor_index": tumor_index,
        "fasta": fasta,
        "fai": fai,
    }


async def check_vcfs(result: MantaSomaticResult, prefix: str) -> None:
    for kind in ("candidate_small_indels", "candidate_sv", "diploid_sv", "somatic_sv"):
        vcf: File = getattr(result, f"{kind}_vcf")
        tbi: File = getattr(result, f"{kind}_vcf_tbi")
        assert vcf.name == f"{prefix}.{kind}.vcf.gz", vcf.name
        assert tbi.name == f"{prefix}.{kind}.vcf.gz.tbi", tbi.name
        async with vcf.open("rb") as fh:
            lines = gzip.decompress(bytes(await fh.read())).decode().splitlines()
        assert "##fileformat=VCFv4.1" in lines, f"{vcf.name} lacks ##fileformat=VCFv4.1"


@env.task
async def test_manta_somatic_cram() -> None:
    # upstream case: manta/somatic "human - cram" (with a manta_options.ini configuration)
    config = Path(tempfile.mkdtemp()) / "manta_options.ini"
    config.write_text("[manta]\nenableRemoteReadRetrievalForInsertionsInGermlineCallingModes = 0\n")
    inputs = await tumour_normal_inputs()
    result = await manta_somatic(**inputs, prefix="test", config=await File.from_local(str(config)))
    await check_vcfs(result, "test")


@env.task
async def test_manta_somatic_cram_bed() -> None:
    # upstream case: manta/somatic "human - cram - bed"
    inputs, bed, bed_tbi = await asyncio.gather(
        tumour_normal_inputs(),
        fixture(CHR21 + "multi_intervals.bed.gz"),
        fixture(CHR21 + "multi_intervals.bed.gz.tbi"),
    )
    result = await manta_somatic(**inputs, prefix="bed_test", target_bed=bed, target_bed_tbi=bed_tbi)
    await check_vcfs(result, "bed_test")


tests = [test_manta_somatic_cram, test_manta_somatic_cram_bed]
