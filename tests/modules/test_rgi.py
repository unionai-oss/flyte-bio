"""Tests for flyte_bio.modules.rgi.

The cases are upstream's rgi/cardannotation, rgi/bwt and rgi/main tests
(modules commit 9a1cc9f, RGI 6.0.8):

- cardannotation and main use the CARD data upstream's snapshot was made
  with. Upstream downloads CARD's ``latest`` release, which has moved on since;
  the snapshot records version 4.0.2, so these tests fetch CARD's versioned
  4.0.2 download (whose files have the snapshot's md5s). As upstream,
  cardannotation checks every file of the processed database by md5, and main
  checks its JSON and table by md5.
- bwt uses upstream's pinned CARD test data (version 4.0.1) and the minigut
  reads from upstream's mag test data. As upstream, it checks the output file
  names (the JSON and tables vary between runs), the intermediate FASTA by
  md5, and the gene mapping table's lines, sorted.

CARD's data is free for academic use only (see flyte_bio.modules.rgi).
"""

from flyte.io import File

from flyte_bio.modules.rgi import RgiDatabase, rgi_bwt, rgi_cardannotation, rgi_main
from tests.framework import assert_md5, env, fixture

CARD_4_0_2 = "https://card.mcmaster.ca/download/0/broadstreet-v4.0.2.tar.bz2"
CARD_TEST_DATA = "genomics/prokaryotes/metagenome/rgi/card-data.tar.bz2"
MINIGUT = "https://raw.githubusercontent.com/nf-core/test-datasets/mag/test_data/test_minigut_"

CARDANNOTATION_MD5S = {
    "CARD-Download-README.txt": "b57364cc1a5fab0541e13a01fabe5dea",
    "PMID.tsv": "6aac0b0e6fc3ee1e1e53cf9251522af3",
    "aro_categories.tsv": "1a45f02fdc6a8b0b64b84ca7cda8dff3",
    "aro_categories_index.tsv": "9a65b210f7c446d7a4e6dc8aa97fd51c",
    "aro_index.tsv": "57163f33db066882971eb4e4df4a062b",
    "card.json": "b73e9f02206e13f02af3809a1ebdf5b5",
    "card_database_v4.0.2.fasta": "500f594b42c6ac97fdbe17a2f788ffd8",
    "card_database_v4.0.2_all.fasta": "aa8bb51771f095348d37c72b59b181e9",
    "nucleotide_fasta_protein_homolog_model.fasta": "e2645d3cfee5bda65d9b7e0956df6463",
    "nucleotide_fasta_protein_knockout_model.fasta": "d8224ce019ad69a041d83ccf7877b4e0",
    "nucleotide_fasta_protein_overexpression_model.fasta": "7d81c25204b97213d7cc7aaf16aba0ba",
    "nucleotide_fasta_protein_variant_model.fasta": "7d2df2efed96e790c4c0d429025a0361",
    "nucleotide_fasta_rRNA_gene_variant_model.fasta": "1321456b509ffa30f3ba69fe631b643a",
    "protein_fasta_protein_homolog_model.fasta": "ba91b130bd7dd8a5eba31e600f7782cc",
    "protein_fasta_protein_knockout_model.fasta": "f16667df70d1a5f910d3ae45f5250c3d",
    "protein_fasta_protein_overexpression_model.fasta": "711e3e9acc87a3c1379e0ddbf6652c92",
    "protein_fasta_protein_variant_model.fasta": "88e30f42525980a63667336a24945900",
    "shortname_antibiotics.tsv": "eab9bb0394aef8c7d04e0453532d2581",
    "shortname_pathogens.tsv": "76cd65af74163663f6fa8dc92ab30e6b",
    "snps.txt": "1f4148e4a8c66f88c99652197bc73023",
}

BWT_TABLES = [
    "test.allele_mapping_data.txt",
    "test.artifacts_mapping_stats.txt",
    "test.coverage.temp.txt",
    "test.coverage_all_positions.summary.temp.txt",
    "test.coverage_all_positions.temp.txt",
    "test.gene_mapping_data.txt",
    "test.overall_mapping_stats.txt",
    "test.reference_mapping_stats.txt",
    "test.seqs.temp.txt",
    "test.temp.txt",
]

EFFLUX_DRUGS = (
    "fluoroquinolone antibiotic; cephalosporin; glycylcycline; penicillin beta-lactam; tetracycline antibiotic; "
    "rifamycin antibiotic; phenicol antibiotic; disinfecting agents and antiseptics"
)
GENE_MAPPING_ROWS = [
    [
        "ARO Term",
        "ARO Accession",
        "Reference Model Type",
        "Reference DB",
        "Alleles with Mapped Reads",
        "Reference Allele(s) Identity to CARD Reference Protein (%)",
        "Resistomes & Variants: Observed in Genome(s)",
        "Resistomes & Variants: Observed in Plasmid(s)",
        "Resistomes & Variants: Observed Pathogen(s)",
        "Completely Mapped Reads",
        "Mapped Reads with Flanking Sequence",
        "All Mapped Reads",
        "Average Percent Coverage",
        "Average Length Coverage (bp)",
        "Average MAPQ (Completely Mapped Reads)",
        "Number of Mapped Baits",
        "Number of Mapped Baits with Reads",
        "Average Number of reads per Bait",
        "Number of reads per Bait Coefficient of Variation (%)",
        "Number of reads mapping to baits and mapping to complete gene",
        "Number of reads mapping to baits and mapping to complete gene (%)",
        "Mate Pair Linkage (# reads)",
        "Reference Length",
        "AMR Gene Family",
        "Drug Class",
        "Resistance Mechanism",
    ],
    [
        "Escherichia coli acrA",
        "3004043",
        "protein homolog model",
        "CARD",
        "1",
        "100.0",
        "no data",
        "no data",
        "Escherichia coli",
        "51.00",
        "0.00",
        "51.00",
        "100.00",
        "1194.00",
        "188.71",
        "0",
        "0",
        "0",
        "0",
        "N/A",
        "N/A",
        "",
        "1194",
        "resistance-nodulation-cell division (RND) antibiotic efflux pump",
        EFFLUX_DRUGS,
        "antibiotic efflux",
    ],
    [
        "Escherichia coli emrE",
        "3004039",
        "protein homolog model",
        "CARD",
        "1",
        "100.0",
        "no data",
        "no data",
        "Escherichia coli",
        "23.00",
        "0.00",
        "23.00",
        "100.00",
        "333.00",
        "175.96",
        "0",
        "0",
        "0",
        "0",
        "N/A",
        "N/A",
        "",
        "333",
        "small multidrug resistance (SMR) antibiotic efflux pump",
        "macrolide antibiotic",
        "antibiotic efflux",
    ],
    [
        "LptD",
        "3005059",
        "protein homolog model",
        "CARD",
        "1",
        "100.0",
        "no data",
        "no data",
        "Klebsiella pneumoniae",
        "38.00",
        "0.00",
        "38.00",
        "44.57",
        "1047.00",
        "131.34",
        "0",
        "0",
        "0",
        "0",
        "N/A",
        "N/A",
        "",
        "2349",
        "ATP-binding cassette (ABC) antibiotic efflux pump",
        "carbapenem; peptide antibiotic; aminocoumarin antibiotic; rifamycin antibiotic",
        "antibiotic efflux",
    ],
    [
        "acrB",
        "3000216",
        "protein homolog model",
        "CARD",
        "1",
        "100.0",
        "no data",
        "no data",
        "Escherichia coli",
        "133.00",
        "0.00",
        "133.00",
        "100.00",
        "3150.00",
        "191.36",
        "0",
        "0",
        "0",
        "0",
        "N/A",
        "N/A",
        "",
        "3150",
        "resistance-nodulation-cell division (RND) antibiotic efflux pump",
        EFFLUX_DRUGS,
        "antibiotic efflux",
    ],
    [
        "kdpE",
        "3003841",
        "protein homolog model",
        "CARD",
        "1",
        "100.0",
        "no data",
        "no data",
        "Escherichia coli",
        "45.00",
        "0.00",
        "45.00",
        "100.00",
        "678.00",
        "175.80",
        "0",
        "0",
        "0",
        "0",
        "N/A",
        "N/A",
        "",
        "678",
        "kdpDE",
        "aminoglycoside antibiotic",
        "antibiotic efflux",
    ],
    [
        "leuO",
        "3003843",
        "protein homolog model",
        "CARD",
        "1",
        "100.0",
        "no data",
        "no data",
        "Escherichia coli",
        "55.00",
        "0.00",
        "55.00",
        "100.00",
        "945.00",
        "173.55",
        "0",
        "0",
        "0",
        "0",
        "N/A",
        "N/A",
        "",
        "945",
        "major facilitator superfamily (MFS) antibiotic efflux pump",
        "nucleoside antibiotic; disinfecting agents and antiseptics",
        "antibiotic efflux",
    ],
    [
        "rmtH",
        "3003198",
        "protein homolog model",
        "CARD",
        "1",
        "100.0",
        "no data",
        "no data",
        "Klebsiella pneumoniae",
        "4.00",
        "0.00",
        "4.00",
        "6.59",
        "50.00",
        "79.00",
        "0",
        "0",
        "0",
        "0",
        "N/A",
        "N/A",
        "",
        "759",
        "16S rRNA methyltransferase (G1405)",
        "aminoglycoside antibiotic",
        "antibiotic target alteration",
    ],
    [
        "tet(Q)",
        "3000191",
        "protein homolog model",
        "CARD",
        "1",
        "100.0",
        "no data",
        "no data",
        "Bacteroides fragilis",
        "138.00",
        "0.00",
        "138.00",
        "100.00",
        "1974.00",
        "179.62",
        "0",
        "0",
        "0",
        "0",
        "N/A",
        "N/A",
        "",
        "1974",
        "tetracycline-resistant ribosomal protection protein",
        "tetracycline antibiotic",
        "antibiotic target protection",
    ],
]
GENE_MAPPING_LINES = ["\t".join(row) for row in GENE_MAPPING_ROWS]


async def text(file: File) -> str:
    async with file.open("rb") as fh:
        return bytes(await fh.read()).decode()


def names(files: list[File]) -> list[str]:
    return [f.name for f in files]


async def card_4_0_2() -> RgiDatabase:
    return await rgi_cardannotation(await fixture(CARD_4_0_2))


@env.task
async def test_rgi_cardannotation() -> None:
    database = await card_4_0_2()
    assert database.db_version == "4.0.2", database.db_version
    files = {f.name: f for f in await database.db.list_files()}
    assert sorted(files) == sorted(CARDANNOTATION_MD5S), f"database files: {sorted(files)}"
    for name, md5 in CARDANNOTATION_MD5S.items():
        await assert_md5(files[name], md5, label=name)


@env.task
async def test_rgi_bwt_minigut_fastq_gz() -> None:
    database = await rgi_cardannotation(await fixture(CARD_TEST_DATA))
    assert database.db_version == "4.0.1", database.db_version
    reads = [await fixture(MINIGUT + "R1.fastq.gz"), await fixture(MINIGUT + "R2.fastq.gz")]
    result = await rgi_bwt(reads, database.db, prefix="test")
    assert names(result.json) == ["test.allele_mapping_data.json"], names(result.json)
    assert names(result.tsv) == BWT_TABLES, names(result.tsv)
    assert names(result.temp) == ["test.temp.sam.temp.fsa"], names(result.temp)
    await assert_md5(result.temp[0], "bd47f4f8e6eefb9f9e70deb5cd42e7f9", label="temp FASTA")
    gene_mapping = next(f for f in result.tsv if f.name == "test.gene_mapping_data.txt")
    lines = sorted((await text(gene_mapping)).splitlines())
    assert lines == GENE_MAPPING_LINES, f"gene mapping lines: {lines}"


@env.task
async def test_rgi_main_haemophilus_influenzae_genome_fna_gz() -> None:
    database = await card_4_0_2()
    genome = await fixture("genomics/prokaryotes/haemophilus_influenzae/genome/genome.fna.gz")
    result = await rgi_main(genome, database.db, prefix="test")
    await assert_md5(result.json, "74b89120ae465892ef013014e44a8ab4", label="json")
    await assert_md5(result.tsv, "9dd0afe6946da9c6608e775f954813a1", label="txt")


tests = [test_rgi_cardannotation, test_rgi_bwt_minigut_fastq_gz, test_rgi_main_haemophilus_influenzae_genome_fna_gz]
