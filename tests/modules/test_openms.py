"""Tests for flyte_bio.modules.openms.

The cases and expected values are upstream's openms/* tests (modules commit
9a1cc9f, OpenMS 3.6.0), with their test configs' ``ext.args`` and prefixes.
Where upstream md5s an output, so do these tests. Where the output records
its run (dates, work paths, run-to-run differences), they check what upstream
snapshots instead:

- FileFilter: the mzML's first 15 lines.
- FeatureFinderIdentification: the feature count (``<featureList count>``)
  and the number of ``<feature id=`` elements, with 2 CPUs as upstream's
  config pins (salting the cache, which ignores resources).
- IDConflictResolver: the number of peptide identifications in the
  consensus elements, one hit each.
- IDMerger: lines the merged idXML must hold.
- PeakPickerHiRes: the centroided spectra's "mass spectrum" term. Upstream's
  test config sets ``-algorithm:ms_levels 2`` under a process name that
  doesn't match the module, so its test runs without args; so does this one.
- PSMFeatureExtractor: the md5 of the first 92 lines (lines joined with
  newlines, no trailing one), which stop before the timestamp.
- TextExporter: the md5 of the ``#`` header lines, joined likewise.

The FeatureFinderIdentification, IDConflictResolver and PSMFeatureExtractor
cases start, as upstream, from a Comet search of a Thermo RAW file against the
human proteome (see tests.modules.test_openmsthirdparty).
"""

import asyncio
import hashlib
import re

import flyte
from flyte.io import File

from flyte_bio.modules.openms import (
    openms_decoydatabase,
    openms_featurefinderidentification_cmd,
    openms_fileconverter,
    openms_filefilter,
    openms_idconflictresolver,
    openms_idfilter,
    openms_idmerger,
    openms_idripper,
    openms_idscoreswitcher,
    openms_mapaligneridentification,
    openms_maprttransformer,
    openms_mztabexporter,
    openms_peakpickerhires,
    openms_peptideindexer,
    openms_pick,
    openms_psmfeatureextractor,
    openms_textexporter,
)
from flyte_bio.modules.thermorawfileparser import thermorawfileparser
from tests.framework import assert_md5, env, fixture
from tests.modules.test_openmsthirdparty import HUMAN_PROTEOME, RAW, openmsthirdparty_comet_search

REP1 = "proteomics/openms/HepG2_rep1_small.idXML"
REP2 = "proteomics/openms/HepG2_rep2_small.idXML"
PEAKPICKER_MZML = "proteomics/msspectra/peakpicker_tutorial_1.mzML"
# Upstream's config pins FeatureFinderIdentification to 2 CPUs; the salt keeps its cache entries apart.
TWO_CPUS = {"resources": flyte.Resources(cpu=2, memory="8Gi"), "cache": flyte.Cache(behavior="auto", salt="cpu=2")}

FILEFILTER_HEAD = [
    '<?xml version="1.0" encoding="ISO-8859-1"?>',
    (
        '<indexedmzML xmlns="http://psi.hupo.org/ms/mzml" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"'
        ' xsi:schemaLocation="http://psi.hupo.org/ms/mzml http://psidev.info/files/ms/mzML/xsd/mzML1.1.0_idx.xsd">'
    ),
    (
        '<mzML xmlns="http://psi.hupo.org/ms/mzml" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"'
        ' xsi:schemaLocation="http://psi.hupo.org/ms/mzml http://psidev.info/files/ms/mzML/xsd/mzML1.1.0.xsd"'
        ' accession="" version="1.1.0">'
    ),
    '\t<cvList count="5">',
    (
        '\t\t<cv id="MS" fullName="Proteomics Standards Initiative Mass Spectrometry Ontology"'
        ' URI="http://psidev.cvs.sourceforge.net/*checkout*/psidev/psi/psi-ms/mzML/controlledVocabulary/psi-ms.obo"/>'
    ),
    '\t\t<cv id="UO" fullName="Unit Ontology" URI="http://obo.cvs.sourceforge.net/obo/obo/ontology/phenotype/unit.obo"/>',
    (
        '\t\t<cv id="BTO" fullName="BrendaTissue545" version="unknown"'
        ' URI="http://www.brenda-enzymes.info/ontology/tissue/tree/update/update_files/BrendaTissueOBO"/>'
    ),
    (
        '\t\t<cv id="GO" fullName="Gene Ontology - Slim Versions" version="unknown"'
        ' URI="http://www.geneontology.org/GO_slims/goslim_goa.obo"/>'
    ),
    (
        '\t\t<cv id="PATO" fullName="Quality ontology" version="unknown"'
        ' URI="http://obo.cvs.sourceforge.net/*checkout*/obo/obo/ontology/phenotype/quality.obo"/>'
    ),
    "\t</cvList>",
    "\t<fileDescription>",
    "\t\t<fileContent>",
    '\t\t\t<cvParam cvRef="MS" accession="MS:1000579" name="MS1 spectrum" />',
    '\t\t\t<cvParam cvRef="MS" accession="MS:1000580" name="MSn spectrum" />',
    "\t\t</fileContent>",
]


async def openms_text(file: File) -> str:
    async with file.open("rb") as fh:
        return bytes(await fh.read()).decode()


def openms_lines_md5(lines: list[str]) -> str:
    return hashlib.md5("\n".join(lines).encode()).hexdigest()


async def openms_test_features() -> File:
    """Upstream's FeatureFinderIdentification setup: features of the RAW file at Comet's identifications."""
    mzml, comet = await openmsthirdparty_comet_search()
    task = openms_featurefinderidentification_cmd.override(**TWO_CPUS)
    results = await task(mzml=mzml, id_int=comet.idxml, id_ext=None, out_type="featureXML", prefix="test", args="")
    return await openms_pick(results, "test.featureXML")


@env.task
async def test_openms_decoydatabase() -> None:
    # upstream case: openms/decoydatabase "test_openms_decoydatabase"
    decoys = await openms_decoydatabase(await fixture("proteomics/database/yeast_UPS.fasta"))
    assert decoys.path.endswith("/yeast_UPS_decoy.fasta"), decoys.path
    await assert_md5(decoys, "96c0295f69e0cb92b705ca4c33e6f047")


@env.task
async def test_openms_featurefinderidentification() -> None:
    # upstream case: openms/featurefinderidentification "proteomics - mzml - idxml"
    features = await openms_test_features()
    text = await openms_text(features)
    count = re.search(r'<featureList count="(\d+)"', text)
    assert count is not None, "no <featureList count>"
    assert count.group(1) == "3949", f"featureList count {count.group(1)}"
    elements = len(re.findall(r'<feature id="', text))
    assert elements == 11847, f"{elements} <feature id= elements"


@env.task
async def test_openms_fileconverter() -> None:
    # upstream case: openms/fileconverter "proteomics - mzML - to mgf"
    mgf = await openms_fileconverter(await fixture(PEAKPICKER_MZML), "mgf", prefix="test")
    assert mgf.path.endswith("/test.mgf"), mgf.path
    await assert_md5(mgf, "2ab0ae45c51f98d1eca2ba862ad74be8")


@env.task
async def test_openms_filefilter() -> None:
    # upstream case: openms/filefilter "filter - mzml"
    mzml = await thermorawfileparser(await fixture(RAW), prefix="test")
    filtered = await openms_filefilter(mzml, prefix="test_filtered", args="-mz 1000:2000")
    assert filtered.path.endswith("/test_filtered.mzML"), filtered.path
    async with filtered.open("rb") as fh:
        head = bytes(await fh.read(1 << 16)).decode(errors="replace")
    lines = head.splitlines()[:15]
    assert lines == FILEFILTER_HEAD, f"mzML lines 1-15 differ: {lines}"


@env.task
async def test_openms_idconflictresolver() -> None:
    # upstream case: openms/idconflictresolver "proteomics - consensusxml"
    consensus = await openms_fileconverter(await openms_test_features(), "consensusXML", prefix="test")
    resolved = await openms_idconflictresolver(consensus, prefix="test_resolved")
    assert resolved.path.endswith("/test_resolved.consensusXML"), resolved.path
    text = await openms_text(resolved)
    section = text[text.index("<consensusElementList>") :]
    elements = section.count("<consensusElement ")
    peptide_ids = section.count("<PeptideIdentification ")
    peptide_hits = section.count("<PeptideHit ")
    assert peptide_hits == peptide_ids, f"{peptide_hits} hits for {peptide_ids} identifications"
    assert peptide_ids <= elements, f"{peptide_ids} identifications for {elements} consensus elements"
    assert peptide_ids == 3949, f"{peptide_ids} peptide identifications"


@env.task
async def test_openms_idfilter() -> None:
    # upstream case: openms/idfilter "proteomics - filter - idxml"
    filtered = await openms_idfilter(await fixture(REP1), await fixture(REP2), prefix="test")
    assert filtered.path.endswith("/test.idXML"), filtered.path
    await assert_md5(filtered, "336992f7531351328ab293ff99d16682")


@env.task
async def test_openms_idmerger() -> None:
    # upstream case: openms/idmerger "proteomics - openms - idxml"
    merged = await openms_idmerger([await fixture(REP1), await fixture(REP2)], prefix="test")
    assert merged.path.endswith("/test.idXML"), merged.path
    lines = (await openms_text(merged)).splitlines()
    for expected in (
        "</PeptideHit>",
        '<VariableModification name="Oxidation (M)" />',
        'accession="sp|P48507|GSH0_HUMAN"',
        'accession="DECOY_sp|Q96SQ7|ATOH8_HUMAN"',
    ):
        assert any(expected in line for line in lines), f"merged idXML lacks {expected!r}"


@env.task
async def test_openms_idripper() -> None:
    # upstream case: openms/idripper "proteomics - split - idxml"
    merged = await openms_idmerger(
        [await fixture(REP1), await fixture(REP2)],
        prefix="test",
        args="-annotate_file_origin true -merge_proteins_add_PSMs",
    )
    ripped = await openms_idripper(merged)
    names = [f.path.rsplit("/", 1)[-1] for f in ripped]
    assert names == ["HepG2_rep1_small.idXML", "HepG2_rep2_small.idXML"], names
    await assert_md5(ripped[0], "983f169dfe2a6f47e2a672f11678875d", label=names[0])
    await assert_md5(ripped[1], "9ae046dd7eb06f87d041c65f3f0c3eda", label=names[1])


@env.task
async def test_openms_idscoreswitcher() -> None:
    # upstream case: openms/idscoreswitcher "proteomics - switch_score - idxml"
    switched = await openms_idscoreswitcher(
        await fixture(REP1), prefix="test", args="-new_score COMET:xcorr -new_score_orientation higher_better"
    )
    await assert_md5(switched, "c08b05d0511454c74560f5a84645d29a")


@env.task
async def test_openms_mapaligneridentification() -> None:
    # upstream case: openms/mapaligneridentification "proteomics - idxml"
    trafos = await openms_mapaligneridentification([await fixture(REP1), await fixture(REP2)])
    names = [f.path.rsplit("/", 1)[-1] for f in trafos]
    assert names == ["HepG2_rep1_small.trafoXML", "HepG2_rep2_small.trafoXML"], names
    await assert_md5(trafos[0], "b70147203181749bb27975929bf6bd33", label=names[0])
    await assert_md5(trafos[1], "c129fc3f79bd40bfdd899689599bc25d", label=names[1])


@env.task
async def test_openms_maprttransformer() -> None:
    # upstream case: openms/maprttransformer "proteomics - idxml - trafoxml"
    rep1, rep2 = await asyncio.gather(fixture(REP1), fixture(REP2))
    trafos = await openms_mapaligneridentification([rep1, rep2])
    aligned = await openms_maprttransformer(rep1, trafos[0], prefix="test_aligned")
    assert aligned.path.endswith("/test_aligned.idXML"), aligned.path
    await assert_md5(aligned, "658a7c43d713ede36f65d2c9dd6523ef")


@env.task
async def test_openms_mztabexporter() -> None:
    # upstream case: openms/mztabexporter "proteomics - idxml"
    mztab = await openms_mztabexporter(await fixture(REP1), prefix="test")
    assert mztab.path.endswith("/test.mzTab"), mztab.path
    await assert_md5(mztab, "5e45450fb6ba4288711299f3753fc4ee")


@env.task
async def test_openms_peakpickerhires() -> None:
    # upstream case: openms/peakpickerhires "proteomics - mzML"
    picked = await openms_peakpickerhires(await fixture(PEAKPICKER_MZML), prefix="test")
    assert picked.path.endswith("/test.mzML"), picked.path
    expected = '\t\t\t<cvParam cvRef="MS" accession="MS:1000294" name="mass spectrum" />'
    assert expected in (await openms_text(picked)).splitlines(), "no 'mass spectrum' term"


@env.task
async def test_openms_peptideindexer() -> None:
    # upstream case: openms/peptideindexer "proteomics - refresh_idxml"
    decoys = await openms_decoydatabase(await fixture(HUMAN_PROTEOME))
    indexed = await openms_peptideindexer(
        await fixture(REP1),
        decoys,
        prefix="test_indexed",
        args="-missing_decoy_action 'silent' -unmatched_action 'remove'",
    )
    await assert_md5(indexed, "0ea90f7f688df1c302f08a1586b7b75d")


@env.task
async def test_openms_psmfeatureextractor() -> None:
    # upstream case: openms/psmfeatureextractor "proteomics - psmfeatureextractor - idxml"
    comet = (await openmsthirdparty_comet_search())[1]
    features = await openms_psmfeatureextractor(comet.idxml, prefix="test_features")
    lines = (await openms_text(features)).splitlines()
    assert any("PeptideHit" in line for line in lines), "no PeptideHit"
    actual = openms_lines_md5(lines[:92])
    assert actual == "0d75087ea4ca67400f9a8f9a49be8444", f"lines 1-92 md5 {actual}"


@env.task
async def test_openms_textexporter() -> None:
    # upstream case: openms/textexporter "proteomics - export - idxml"
    tsv = await openms_textexporter(await fixture(REP1), prefix="test")
    lines = (await openms_text(tsv)).splitlines()
    assert any("PROTEIN" in line for line in lines), "no PROTEIN rows"
    assert any("PEPTIDE" in line for line in lines), "no PEPTIDE rows"
    actual = openms_lines_md5([line for line in lines if line.startswith("#")])
    assert actual == "ed2ef59b2d596b75b2eccd9ebc735bff", f"header lines md5 {actual}"


tests = [
    test_openms_decoydatabase,
    test_openms_featurefinderidentification,
    test_openms_fileconverter,
    test_openms_filefilter,
    test_openms_idconflictresolver,
    test_openms_idfilter,
    test_openms_idmerger,
    test_openms_idripper,
    test_openms_idscoreswitcher,
    test_openms_mapaligneridentification,
    test_openms_maprttransformer,
    test_openms_mztabexporter,
    test_openms_peakpickerhires,
    test_openms_peptideindexer,
    test_openms_psmfeatureextractor,
    test_openms_textexporter,
]
