"""Tests for flyte_bio.modules.tximport.

Expected md5s are the upstream tximeta/tximport snapshot values; the output
tables are a pure function of the quant files and the tx2gene table.
"""

from flyte.io import File

from flyte_bio.modules.tx2gene import tx2gene
from flyte_bio.modules.tximport import TximportResult, tximport
from tests.framework import assert_md5, env, fixture, fixture_dir

YEAST = "genomics/eukaryotes/saccharomyces_cerevisiae/"


async def assert_tximport_md5s(result: TximportResult, expected: dict[str, str], label: str) -> None:
    for field, md5 in expected.items():
        output: File = getattr(result, field)
        await assert_md5(output, md5, label=f"{label} {field}")


@env.task
async def test_tximport_kallisto() -> None:
    # upstream case: tximeta/tximport "saccharomyces_cerevisiae - kallisto - gtf"
    gtf = await fixture(YEAST + "genome_gfp.gtf")
    quants = await fixture_dir(YEAST + "kallisto_results.tar.gz", "kallisto_results")
    tx2gene_tsv = await tx2gene(gtf, quants, quant_type="kallisto", gene_id="gene_id", extra="gene_name")
    result = await tximport(quants, tx2gene_tsv, quant_type="kallisto", prefix="test")
    await assert_tximport_md5s(
        result,
        {
            "counts_gene": "e89c28692ea214396b2d4cb702a804c3",
            "counts_gene_length_scaled": "4944841ac711124d29673b6b6ed16ef3",
            "counts_gene_scaled": "39d14e361434978b3cadae901a26a028",
            "counts_transcript": "42e0106e75fa97c1c684c6d9060f1724",
            "lengths_gene": "db6becdf807fd164a9c63dd1dd916d9c",
            "lengths_transcript": "f974b52840431a5dae57bcb615badbf1",
            "tpm_gene": "85d108269769ae0d841247b9b9ed922d",
            "tpm_transcript": "65862ed9d4a05abfab952e680dc0e49d",
            "tx2gene_augmented": "0e2418a69d2eba45097ebffc2f700bfe",
        },
        label="kallisto",
    )


@env.task
async def test_tximport_salmon() -> None:
    # upstream case: tximeta/tximport "saccharomyces_cerevisiae - salmon - gtf"
    gtf = await fixture(YEAST + "genome_gfp.gtf")
    quants = await fixture_dir(YEAST + "salmon_results.tar.gz")
    tx2gene_tsv = await tx2gene(gtf, quants, quant_type="salmon", gene_id="gene_id", extra="gene_name")
    result = await tximport(quants, tx2gene_tsv, quant_type="salmon", prefix="test")
    await assert_tximport_md5s(
        result,
        {
            "counts_gene": "c14cab7e15cfac73ec0602dc2c404551",
            "counts_gene_length_scaled": "5f92a6784f6edc5e3b336c71c3ee7daf",
            "counts_gene_scaled": "fdfb3d23aaf5d4316d81247ec4664ca0",
            "counts_transcript": "ff0f5be09ca7a322672c0074ba35da17",
            "lengths_gene": "1691ea2677612805cd699265c83024d7",
            "lengths_transcript": "db6d8ab9f8e1123d5984fd534b4347dc",
            "tpm_gene": "6076364cc78741a4f8bc8935a045d13d",
            "tpm_transcript": "7a334b565e1e865efb1caf615f194ef7",
            "tx2gene_augmented": "0e2418a69d2eba45097ebffc2f700bfe",
        },
        label="salmon",
    )


tests = [test_tximport_kallisto, test_tximport_salmon]
