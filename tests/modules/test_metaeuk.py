"""Tests for flyte_bio.modules.metaeuk.

The cases are upstream's metaeuk/easypredict tests (modules commit 9a1cc9f,
MetaEuk 6.a5d39d9) on the human test genome: against the yeast UPS protein
FASTA, checking upstream's md5s, and against the Kalamari MMseqs2 database
(downloaded by :func:`mmseqs_databases`, as upstream's setup does), checking
the output names only, as upstream marks those outputs unstable.
"""

from flyte_bio.modules.metaeuk import MetaeukResult, metaeuk_easypredict
from flyte_bio.modules.mmseqs import mmseqs_databases
from tests.framework import assert_md5, env, fixture

GENOME = "genomics/homo_sapiens/genome/genome.fasta"


def check_names(result: MetaeukResult) -> None:
    names = [result.faa.name, result.codon.name, result.tsv.name, result.gff.name]
    assert names == ["test.fas", "test.codon.fas", "test.headersMap.tsv", "test.gff"], f"outputs: {names}"


@env.task
async def test_metaeuk_easypredict_fasta() -> None:
    result = await metaeuk_easypredict(
        await fixture(GENOME), await fixture("proteomics/database/yeast_UPS.fasta"), prefix="test"
    )
    check_names(result)
    await assert_md5(result.codon, "6eed1f7651d7ba82171cb02014884e5f", label="test.codon.fas")
    await assert_md5(result.faa, "99e9062723af0b474c21de91a4f5dab0", label="test.fas")
    await assert_md5(result.gff, "be171e6f1e1e546909ad540ff0a9ffe4", label="test.gff")
    await assert_md5(result.tsv, "b3628c76d7385d326bd8a7931d974cff", label="test.headersMap.tsv")


@env.task
async def test_metaeuk_easypredict_mmseqs_db() -> None:
    database = await mmseqs_databases("Kalamari")
    result = await metaeuk_easypredict(await fixture(GENOME), database, prefix="test")
    check_names(result)


tests = [
    test_metaeuk_easypredict_fasta,
    test_metaeuk_easypredict_mmseqs_db,
]
