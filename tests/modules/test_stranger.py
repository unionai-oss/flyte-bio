"""Tests for flyte_bio.modules.stranger.

The case and expected md5s are upstream's stranger test (modules commit
9a1cc9f, stranger 0.10.2, htslib 1.23.1), with its test config's prefix
``stranger``. As upstream, its input is ExpansionHunter's genotyping of the
human test BAM against the test repeat catalogue; ExpansionHunter isn't
ported, so its upstream command runs here as a test setup task. The VCF is
compared decompressed (as nf-test hashes it), the index byte for byte.
"""

import asyncio

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS
from flyte_bio.modules.stranger import stranger
from tests.framework import assert_gunzipped_md5, assert_md5, env, fixture

HUMAN = "genomics/homo_sapiens/"
CATALOG = HUMAN + "genome/loci/repeat_expansions.json"

# Upstream's expansionhunter command, for the stranger test's setup.
stranger_setup_expansionhunter_cmd = shell.create(
    name="test_stranger_setup_expansionhunter",
    image="community.wave.seqera.io/library/expansionhunter:5.0.0--389ada7e191a4fba",
    resources=flyte.Resources(cpu=2, memory="4Gi"),
    inputs={"bam": File, "bai": File, "fasta": File, "fai": File, "variant_catalog": File, "prefix": str},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        BAM=({{inputs.bam}}); BAI=({{inputs.bai}}); FA=({{inputs.fasta}}); FAI=({{inputs.fai}})
        CATALOG=({{inputs.variant_catalog}})
        P={{inputs.prefix}}
        W=$(mktemp -d); cd "$W"
        for F in "${{BAM[0]}}" "${{BAI[0]}}" "${{FA[0]}}" "${{FAI[0]}}" "${{CATALOG[0]}}"; do
            ln -s "$F" "$(basename "$F")"
        done
        ExpansionHunter \
            --threads $TASK_CPUS \
            --reads "$(basename "${{BAM[0]}}")" \
            --output-prefix "$P" \
            --reference "$(basename "${{FA[0]}}")" \
            --variant-catalog "$(basename "${{CATALOG[0]}}")" || exit 1
        bgzip --threads $TASK_CPUS "$P.vcf" || exit 1
        mv "$P.vcf.gz" {{outputs.results}}/
    """,
)
env.add_dependency(flyte.TaskEnvironment.from_task("test_stranger_setup", stranger_setup_expansionhunter_cmd.as_task()))


@env.task
async def test_stranger_repeat_catalogue() -> None:
    # upstream case: stranger "homo_sapiens - vcf, repeat_catalogue"
    bam, bai, fasta, fai, catalog = await asyncio.gather(
        fixture(HUMAN + "illumina/bam/test.paired_end.sorted.bam"),
        fixture(HUMAN + "illumina/bam/test.paired_end.sorted.bam.bai"),
        fixture(HUMAN + "genome/genome.fasta"),
        fixture(HUMAN + "genome/genome.fasta.fai"),
        fixture(CATALOG),
    )
    called = await stranger_setup_expansionhunter_cmd(
        bam=bam, bai=bai, fasta=fasta, fai=fai, variant_catalog=catalog, prefix="test"
    )
    vcf = await called.get_file("test.vcf.gz")
    assert vcf is not None, "expansionhunter wrote no test.vcf.gz"
    result = await stranger(vcf, catalog, prefix="stranger")
    assert result.vcf.path.endswith("/stranger.vcf.gz"), result.vcf.path
    await assert_gunzipped_md5(result.vcf, "7bccff6da8e34606dccad76c44434e76", label="stranger.vcf.gz")
    await assert_md5(result.tbi, "587822ec191e1840e306146b7e5dafcd", label="stranger.vcf.gz.tbi")


tests = [test_stranger_repeat_catalogue]
