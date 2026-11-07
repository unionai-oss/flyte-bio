"""Tests for flyte_bio.modules.thermorawfileparser.

The case is upstream's thermorawfileparser test (modules commit 9a1cc9f,
ThermoRawFileParser 1.4.5): an Orbitrap Elite RAW file converted to mzML
with no ``args``. The mzML records the input's path near its top and a
checksum at its end, so, as upstream, the test compares lines 22 to 51 (the
instrument, software and processing header) rather than an md5.
"""

from flyte_bio.modules.thermorawfileparser import thermorawfileparser
from tests.framework import env, fixture

RAW = "proteomics/msspectra/PXD012083_e005640_II.raw"

EXPECTED_HEADER = [
    '    <referenceableParamGroupList count="1">',
    '      <referenceableParamGroup id="commonInstrumentParams">',
    '        <cvParam cvRef="MS" accession="MS:1001910" value="" name="LTQ Orbitrap Elite" />',
    '        <cvParam cvRef="MS" accession="MS:1000529" value="SN03142B" name="instrument serial number" />',
    "      </referenceableParamGroup>",
    "    </referenceableParamGroupList>",
    '    <softwareList count="1">',
    '      <software id="ThermoRawFileParser" version="1.4.5">',
    '        <cvParam cvRef="MS" accession="MS:1003145" name="ThermoRawFileParser" />',
    "      </software>",
    "    </softwareList>",
    '    <instrumentConfigurationList count="1">',
    '      <instrumentConfiguration id="IC1">',
    '        <referenceableParamGroupRef ref="commonInstrumentParams" />',
    '        <componentList count="3">',
    '          <source order="1">',
    '            <cvParam cvRef="MS" accession="MS:1000073" value="" name="electrospray ionization" />',
    "          </source>",
    '          <analyzer order="2">',
    '            <cvParam cvRef="MS" accession="MS:1000484" value="" name="orbitrap" />',
    "          </analyzer>",
    '          <detector order="3">',
    '            <cvParam cvRef="MS" accession="MS:1000624" value="" name="inductive detector" />',
    "          </detector>",
    "        </componentList>",
    "      </instrumentConfiguration>",
    "    </instrumentConfigurationList>",
    '    <dataProcessingList count="1">',
    '      <dataProcessing id="ThermoRawFileParserProcessing">',
    '        <processingMethod order="0" softwareRef="ThermoRawFileParser">',
]


@env.task
async def test_thermorawfileparser_parse_raw() -> None:
    spectra = await thermorawfileparser(await fixture(RAW), prefix="test")
    assert spectra.name == "test.mzML", f"spectra file is {spectra.name}"
    async with spectra.open("rb") as fh:
        head = bytes(await fh.read(1 << 16)).decode(errors="replace")
    lines = head.splitlines()[21:51]
    assert lines == EXPECTED_HEADER, f"mzML lines 22-51 differ: {lines}"


tests = [
    test_thermorawfileparser_parse_raw,
]
