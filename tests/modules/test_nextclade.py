"""Tests for flyte_bio.modules.nextclade.

The cases are upstream's nextclade/datasetget and nextclade/run tests
(modules commit 9a1cc9f, Nextclade 3.21.2), on the SARS-CoV-2 dataset
pinned to release ``2024-01-16--20-31-02Z``: every dataset file's md5, and
every Nextclade output's md5 except the JSON, of which lines 6-12 are
checked, as upstream. Upstream has no errors or insertions CSV for this
genome; neither does this run. The dataset is staged as ``orfs``, the name
upstream's run sees (it is recorded in the NDJSON).
"""

from flyte.io import Dir

from flyte_bio.modules.nextclade import nextclade_datasetget, nextclade_run
from tests.framework import assert_md5, env, fixture

DATASET = "nextstrain/sars-cov-2/wuhan-hu-1/orfs"
TAG = "2024-01-16--20-31-02Z"
DATASET_MD5S = {
    "CHANGELOG.md": "ebbe8be5a3c378ed903c1afb4d8c441d",
    "README.md": "c69387d632361334f0d7c9b66065f947",
    "genome_annotation.gff3": "4dff84d2d6ada820e0e3a8bc6798d402",
    "pathogen.json": "db5bbec52359c1e168ffc5e6dc0ea32a",
    "reference.fasta": "c7ce05f28e4ec0322c96f24e064ef55c",
    "sequences.fasta": "c2a4d6cbb837dce22d81f9c36dd0629e",
    "tree.json": "e180607cd34a6cb6bab101d295f6cedf",
}
RUN_MD5S = {
    "test.csv": "547159941c3ac2a5b7525c90d19b1226",
    "test.aligned.fasta": "1bf54662837b0df37f1857c7fa631225",
    "test.auspice.json": "192080dfa0ca85b193db68769116aa33",
    "test.ndjson": "72c53dae4348f21974ba3a57d9c26032",
    "test.nwk": "fae09e1063d673a4db750f605b7de587",
    "test.tsv": "085b4e0400368749b49a26a6caa441f4",
}
TRANSLATION_MD5S = {
    "test.cds_translation.E.fasta": "1a6d93bd7abfeb193476a86950f07202",
    "test.cds_translation.M.fasta": "6f79ab0742c078fcd9d2a474518c7022",
    "test.cds_translation.N.fasta": "bef9912f101777bfff36225d7e5c3c1f",
    "test.cds_translation.ORF1a.fasta": "4dc0e2eb1f2a61939ba49500ee3fa41a",
    "test.cds_translation.ORF1b.fasta": "9c2e83d26161b5c887ff51cd64bd15bb",
    "test.cds_translation.ORF3a.fasta": "adcbede4ebc2cac7af755a9d29d28ea3",
    "test.cds_translation.ORF6.fasta": "3785b34cce978c95256f83db6ee82af0",
    "test.cds_translation.ORF7a.fasta": "6ede1acb9e75afc84aa30bbc40551d37",
    "test.cds_translation.ORF7b.fasta": "460e4cbc5f8c632c2bc9a8aedad5cf43",
    "test.cds_translation.ORF8.fasta": "c733c88e61b29542664368fbf6dd4c76",
    "test.cds_translation.ORF9b.fasta": "0aa13afc6cbf445fc92caa2e6c0a7548",
    "test.cds_translation.S.fasta": "77740927a3f00b7e5bfac392fa6d264c",
}
JSON_LINES = [
    '  "cladeNodeAttrKeys": [',
    "    {",
    '      "name": "Nextclade_pango",',
    '      "displayName": "Pango lineage (Nextclade)",',
    (
        '      "description": "Pango lineage as inferred by Nextclade from the nearest neighbour in the reference '
        "tree. 98% accurate for recent sequences, for higher accuracy use dedicated pangolin software in UShER or "
        "pangoLEARN mode. Recombinants may get (wrongly) assigned to a designated recombinant lineage if they have "
        'similar breakpoints.",'
    ),
    '      "hideInWeb": false,',
    '      "skipAsReference": false',
]


async def check_md5s(directory: Dir, md5s: dict[str, str]) -> None:
    for name, md5 in md5s.items():
        found = await directory.get_file(name)
        assert found is not None, f"no {name}"
        await assert_md5(found, md5, label=name)


@env.task
async def test_nextclade_datasetget() -> None:
    # upstream case: nextclade/datasetget "species-tag"
    dataset = await nextclade_datasetget(DATASET, TAG)
    names = sorted(f.name for f in await dataset.list_files())
    assert names == sorted(DATASET_MD5S), f"dataset files {names}"
    await check_md5s(dataset, DATASET_MD5S)


@env.task
async def test_nextclade_run() -> None:
    # upstream case: nextclade/run "sarscov2 default"
    dataset = await nextclade_datasetget(DATASET, TAG)
    result = await nextclade_run(
        await fixture("genomics/sarscov2/genome/genome.fasta"), dataset, dataset_name="orfs", prefix="test"
    )
    assert result.csv_errors is None and result.csv_insertions is None
    await check_md5s(result.results, RUN_MD5S)
    assert [f.name for f in result.fasta_translation] == sorted(TRANSLATION_MD5S)
    await check_md5s(result.results, TRANSLATION_MD5S)
    assert result.json is not None
    async with result.json.open("rb") as fh:
        lines = bytes(await fh.read()).decode().splitlines()[5:12]
    assert lines == JSON_LINES, f"test.json lines 6-12: {lines}"


tests = [test_nextclade_datasetget, test_nextclade_run]
