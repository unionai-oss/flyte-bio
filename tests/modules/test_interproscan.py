"""Tests for flyte_bio.modules.interproscan.

The cases and expected values are upstream's interproscan tests (modules
commit 9a1cc9f, InterProScan 5.59-91.0) with their configs' ``args``. As
upstream, the test-database case downloads upstream's small InterProScan data
archive with aria2 and runs TIGRFAM and SFLD against it; the sarscov2 cases
run Coils on the container's sample data. Each checks the first TSV row, the
XML's and JSON's md5s, and the GFF3 without its sixth line (a match,
which records the run's date), as upstream does. The stop-codon case appends
a ``*`` to every protein and checks the TSV's first row is unchanged.
"""

import tempfile
from pathlib import Path

from flyte.io import File

from flyte_bio.modules.aria2 import aria2
from flyte_bio.modules.interproscan import InterproscanResult, interproscan
from tests.framework import TEST_DATA_BASE, assert_md5, env, fixture

PROTEOME = "genomics/sarscov2/genome/proteome.fasta"
DATABASE_ARGS = "--applications TIGRFAM,sfld --disable-precalc"
COILS_ARGS = "-appl Coils"
DATABASE_TSV_ROW = (
    "GI|225038609|EFDID|719595|FULL\t079fff43a0270e432d339ea71b6f0acf\t350\tSFLD\tSFLDS00057\t"
    "Glutaminase/Asparaginase\t17\t347\t0.0\tT"
)
COILS_TSV_ROW = "ENSSASP00005000004.1\t4c35f09aac2f7be4f3cffd30c6aecac8\t1273\tCoils\tCoil\tCoil\t1176\t1203\t-\tT"
DATABASE_GFF3 = [
    "##gff-version 3",
    "##feature-ontology http://song.cvs.sourceforge.net/viewvc/song/ontology/sofa.obo?revision=1.269",
    "##interproscan-version 5.59-91.0",
    "##sequence-region GI|225038609|EFDID|719595|FULL 1 350",
    "GI|225038609|EFDID|719595|FULL\t.\tpolypeptide\t1\t350\t.\t+\t.\tID=GI|225038609|EFDID|719595|FULL;md5=079fff43a0270e432d339ea71b6f0acf",
    "##FASTA",
    ">GI|225038609|EFDID|719595|FULL",
    "MDKLLWNKKEISVSERISMKRIVMINTGGTFSSQRGENGLSPKLTGTQIRNFLGEFEEDL",
    "ELSTEDYCALDSSNITPEDWVQLADKISQIIYSCDGVVIIHGTDTMAYTASMLSFMLQNL",
    "PIPVVLTGSQLPIGVPMSDAVNNCRCAVQMAASGLGGVYVAFDHKLMLGCRTSKVRTVSF",
    "NAFESINYPYVGEVNALGMQLYPTRLSKPTGEFQLQTAYSDKIAVLKLFPGMRPDLFSFL",
    "QEKGYEGIYIEGFGLGGVPFVKNDITEEISKASKAGIPILVGSQCSYEGSNLGIYETGLR",
    "VLESGGIPVHDMTQEAIVTKLMWCLGQTKDREKIHQLFHTNLIQEVTLPY",
    ">match$1_17_347",
    "ISMKRIVMINTGGTFSSQRGENGLSPKLTGTQIRNFLGEFEEDLELSTEDYCALDSSNIT",
    "PEDWVQLADKISQIIYSCDGVVIIHGTDTMAYTASMLSFMLQNLPIPVVLTGSQLPIGVP",
    "MSDAVNNCRCAVQMAASGLGGVYVAFDHKLMLGCRTSKVRTVSFNAFESINYPYVGEVNA",
    "LGMQLYPTRLSKPTGEFQLQTAYSDKIAVLKLFPGMRPDLFSFLQEKGYEGIYIEGFGLG",
    "GVPFVKNDITEEISKASKAGIPILVGSQCSYEGSNLGIYETGLRVLESGGIPVHDMTQEA",
    "IVTKLMWCLGQTKDREKIHQLFHTNLIQEVT",
]
COILS_GFF3 = [
    "##gff-version 3",
    "##feature-ontology http://song.cvs.sourceforge.net/viewvc/song/ontology/sofa.obo?revision=1.269",
    "##interproscan-version 5.59-91.0",
    "##sequence-region ENSSASP00005000004.1 1 1273",
    "ENSSASP00005000004.1\t.\tpolypeptide\t1\t1273\t.\t+\t.\tID=ENSSASP00005000004.1;md5=4c35f09aac2f7be4f3cffd30c6aecac8",
    "##FASTA",
    ">ENSSASP00005000004.1",
    "MFVFLVLLPLVSSQCVNLTTRTQLPPAYTNSFTRGVYYPDKVFRSSVLHSTQDLFLPFFS",
    "NVTWFHAIHVSGTNGTKRFDNPVLPFNDGVYFASTEKSNIIRGWIFGTTLDSKTQSLLIV",
    "NNATNVVIKVCEFQFCNDPFLGVYYHKNNKSWMESEFRVYSSANNCTFEYVSQPFLMDLE",
    "GKQGNFKNLREFVFKNIDGYFKIYSKHTPINLVRDLPQGFSALEPLVDLPIGINITRFQT",
    "LLALHRSYLTPGDSSSGWTAGAAAYYVGYLQPRTFLLKYNENGTITDAVDCALDPLSETK",
    "CTLKSFTVEKGIYQTSNFRVQPTESIVRFPNITNLCPFGEVFNATRFASVYAWNRKRISN",
    "CVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIAD",
    "YNYKLPDDFTGCVIAWNSNNLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQAGSTPC",
    "NGVEGFNCYFPLQSYGFQPTNGVGYQPYRVVVLSFELLHAPATVCGPKKSTNLVKNKCVN",
    "FNFNGLTGTGVLTESNKKFLPFQQFGRDIADTTDAVRDPQTLEILDITPCSFGGVSVITP",
    "GTNTSNQVAVLYQDVNCTEVPVAIHADQLTPTWRVYSTGSNVFQTRAGCLIGAEHVNNSY",
    "ECDIPIGAGICASYQTQTNSPRRARSVASQSIIAYTMSLGAENSVAYSNNSIAIPTNFTI",
    "SVTTEILPVSMTKTSVDCTMYICGDSTECSNLLLQYGSFCTQLNRALTGIAVEQDKNTQE",
    "VFAQVKQIYKTPPIKDFGGFNFSQILPDPSKPSKRSFIEDLLFNKVTLADAGFIKQYGDC",
    "LGDIAARDLICAQKFNGLTVLPPLLTDEMIAQYTSALLAGTITSGWTFGAGAALQIPFAM",
    "QMAYRFNGIGVTQNVLYENQKLIANQFNSAIGKIQDSLSSTASALGKLQDVVNQNAQALN",
    "TLVKQLSSNFGAISSVLNDILSRLDKVEAEVQIDRLITGRLQSLQTYVTQQLIRAAEIRA",
    "SANLAATKMSECVLGQSKRVDFCGKGYHLMSFPQSAPHGVVFLHVTYVPAQEKNFTTAPA",
    "ICHDGKAHFPREGVFVSNGTHWFVTQRNFYEPQIITTDNTFVSGNCDVVIGIVNNTVYDP",
    "LQPELDSFKEELDKYFKNHTSPDVDLGDISGINASVVNIQKEIDRLNEVAKNLNESLIDL",
    "QELGKYEQYIKWPWYIWLGFIAGLIAIVMVTIMLCCMTSCCSCLKGCCSCGSCCKFDEDD",
    "SEPVLKGVKLHYT",
    ">match$1_1176_1203",
    "VVNIQKEIDRLNEVAKNLNESLIDLQEL",
]


async def text_lines(file: File | None, label: str) -> list[str]:
    assert file is not None, f"interproscan wrote no {label}"
    async with file.open("rb") as fh:
        return bytes(await fh.read()).decode().splitlines()


async def check_result(result: InterproscanResult, tsv_row: str, xml_md5: str, json_md5: str, gff3: list[str]) -> None:
    tsv = await text_lines(result.tsv, "tsv")
    assert tsv_row in tsv[0], f"first TSV row is {tsv[0]!r}"
    assert result.xml is not None and result.json is not None, "interproscan wrote no xml or json"
    await assert_md5(result.xml, xml_md5, label="xml")
    await assert_md5(result.json, json_md5, label="json")
    actual = await text_lines(result.gff3, "gff3")
    actual = actual[0:5] + actual[6:]
    assert actual == gff3, f"gff3 differs: {actual}"


@env.task
async def test_interproscan_l_asparaginase_test_database() -> None:
    data = await aria2(TEST_DATA_BASE + "proteomics/interproscan/interproscan_test.tar.gz")
    faa = await fixture("proteomics/interproscan/l_arginase.faa")
    result = await interproscan(faa, data=data, prefix="test", args=DATABASE_ARGS)
    await check_result(
        result, DATABASE_TSV_ROW, "7248992d9c1618cf7baa7515ae79ce32", "e0d127dd8a952cbd798999851d1338e6", DATABASE_GFF3
    )


@env.task
async def test_interproscan_sarscov2_proteome() -> None:
    result = await interproscan(await fixture(PROTEOME), prefix="test", args=COILS_ARGS)
    await check_result(
        result, COILS_TSV_ROW, "7a211c1a4761e2b9b8700e6e9abbb15f", "b05cffc28b7bfeb3dabe43c2927b2024", COILS_GFF3
    )


@env.task
async def test_interproscan_sarscov2_proteome_gz() -> None:
    result = await interproscan(await fixture(PROTEOME + ".gz"), prefix="test", args=COILS_ARGS)
    await check_result(
        result, COILS_TSV_ROW, "7a211c1a4761e2b9b8700e6e9abbb15f", "b05cffc28b7bfeb3dabe43c2927b2024", COILS_GFF3
    )


def with_stop_codons(fasta: str) -> str:
    """``fasta`` with each sequence on one line and ending in ``*``, as upstream's test builds it."""
    records = []
    for chunk in Path(fasta).read_text().split(">")[1:]:
        header, *sequence = chunk.splitlines()
        records.append(f">{header}\n{''.join(sequence)}*\n")
    out = Path(tempfile.mkdtemp()) / "proteome_with_stops.fasta"
    out.write_text("".join(records))
    return str(out)


@env.task
async def test_interproscan_sarscov2_proteome_stop_codons() -> None:
    proteome = await (await fixture(PROTEOME)).download()
    fasta = await File.from_local(with_stop_codons(str(proteome)))
    result = await interproscan(fasta, prefix="test", args=COILS_ARGS)
    tsv = await text_lines(result.tsv, "tsv")
    assert COILS_TSV_ROW in tsv[0], f"first TSV row is {tsv[0]!r}"


tests = [
    test_interproscan_l_asparaginase_test_database,
    test_interproscan_sarscov2_proteome,
    test_interproscan_sarscov2_proteome_gz,
    test_interproscan_sarscov2_proteome_stop_codons,
]
