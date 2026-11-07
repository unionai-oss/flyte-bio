"""Tests for flyte_bio.modules.snpeff.

The cases and expected md5s are upstream's snpeff/download and snpeff/snpeff
tests (modules commit 9a1cc9f, snpEff 5.4c). Upstream downloads the
C. elegans database ``WBcel235.99`` and snapshots its
``snpEffectPredictor.bin``.

Upstream annotates sarscov2's ``test.vcf`` with ``WBcel235.99`` (its test
config's ``snpeff_db``), passing as cache a bucket copy of ``WBcel235.105``;
that cache doesn't hold ``WBcel235.99``, so snpEff downloads it itself. Here
the annotation reads :func:`snpeff_download`'s cached ``WBcel235.99`` instead.
Its snapshot records the genes table's md5, the other outputs' names and
nft-vcf's ``variantsMD5`` of the annotated VCF (whose header records snpEff's
command line).
"""

import hashlib

from flyte.io import Dir, File

from flyte_bio.modules.snpeff import snpeff_download, snpeff_snpeff
from tests.framework import assert_md5, env, fixture

DB = "WBcel235.99"
VCF = "genomics/sarscov2/illumina/vcf/test.vcf"


async def snpeff_cache() -> Dir:
    return await snpeff_download(DB)


async def text_lines(file: File) -> list[str]:
    async with file.open("rb") as fh:
        return bytes(await fh.read()).decode().splitlines()


async def variants_md5(file: File) -> str:
    """nft-vcf's ``variantsMD5``: md5 of the record lines joined without newlines."""
    body = "".join(line for line in await text_lines(file) if not line.startswith("#"))
    return hashlib.md5(body.encode()).hexdigest()


@env.task
async def test_snpeff_download() -> None:
    cache = await snpeff_cache()
    files = {f.path.split("/" + DB + "/", 1)[-1]: f async for f in cache.walk() if "/" + DB + "/" in f.path}
    assert sorted(files) == ["snpEffectPredictor.bin"], f"cache files: {sorted(files)}"
    await assert_md5(
        files["snpEffectPredictor.bin"], "07eca39b2a3176a0292d11d6256c216d", label="snpEffectPredictor.bin"
    )


@env.task
async def test_snpeff_snpeff_annotate() -> None:
    result = await snpeff_snpeff(await fixture(VCF), DB, cache=await snpeff_cache(), prefix="test")
    assert result.genes_txt is not None and result.genes_txt.name == "test.genes.txt", result.genes_txt
    await assert_md5(result.genes_txt, "130536bf0237d7f3f746d32aaa32840a", label="test.genes.txt")
    assert result.report is not None and result.report.name == "test.csv", result.report
    assert result.summary_html is not None and result.summary_html.name == "snpEff_summary.html", result.summary_html
    assert result.vcf.name == "test.ann.vcf", result.vcf
    md5 = await variants_md5(result.vcf)
    assert md5 == "fd266e623365c7d1837bc9eb088fe70c", f"variantsMD5 {md5}"


tests = [test_snpeff_download, test_snpeff_snpeff_annotate]
