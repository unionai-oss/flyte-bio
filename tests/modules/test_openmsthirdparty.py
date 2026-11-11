"""Tests for flyte_bio.modules.openmsthirdparty.

The case is upstream's openmsthirdparty/cometadapter test (modules commit
9a1cc9f, OpenMS 3.6.0, Comet 2024.01 rev. 1), with its test config's
``ext.args``. As upstream, a Thermo RAW file is converted to mzML
(ThermoRawFileParser) and searched against the human proteome with decoys
added (DecoyDatabase). The idXML records a timestamp, so, as upstream, the
test compares its first 31 lines (the search parameters) and checks that it
holds protein and peptide hits.

Upstream's tests run with at most 4 CPUs, so the search pins 4 CPUs, salting
the cache (which ignores resources). :func:`openmsthirdparty_comet_search` is
shared with the OpenMS tests that start from Comet's identifications.
"""

import flyte
from flyte.io import File

from flyte_bio.modules.openms import openms_decoydatabase
from flyte_bio.modules.openmsthirdparty import (
    CometAdapterResult,
    comet_result,
    openmsthirdparty_cometadapter_cmd,
)
from flyte_bio.modules.thermorawfileparser import thermorawfileparser
from tests.framework import env, fixture

RAW = "proteomics/msspectra/PXD012083_e005640_II.raw"
HUMAN_PROTEOME = "proteomics/database/UP000005640_9606.fasta"
COMET_ARGS = (
    "-instrument low_res -fragment_bin_offset 0.4 -precursor_mass_tolerance 5 -precursor_error_units 'ppm' "
    "-fragment_mass_tolerance 0.50025 -digest_mass_range '800:5000' -max_variable_mods_in_peptide 1 "
    "-precursor_charge '2:5' -activation_method 'CID' -variable_modifications 'Oxidation (M)' "
    "-enzyme 'unspecific cleavage' -spectrum_batch_size 0"
)
FOUR_CPUS = {"resources": flyte.Resources(cpu=4, memory="16Gi"), "cache": flyte.Cache(behavior="auto", salt="cpu=4")}

EXPECTED_HEAD = [
    '<?xml version="1.0" encoding="UTF-8"?>',
    '<?xml-stylesheet type="text/xsl" href="https://www.openms.de/xml-stylesheet/IdXML.xsl" ?>',
    (
        '<IdXML version="1.5" xsi:noNamespaceSchemaLocation="https://www.openms.de/xml-schema/IdXML_1_5.xsd"'
        ' xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">'
    ),
    (
        '\t<SearchParameters id="SP_0" db="UP000005640_9606_decoy.fasta" db_version="" taxonomy=""'
        ' mass_type="monoisotopic" charges="2:5" enzyme="unspecific cleavage" missed_cleavages="1"'
        ' precursor_peak_tolerance="0" precursor_peak_tolerance_ppm="true" peak_mass_tolerance="0.50025"'
        ' peak_mass_tolerance_ppm="false" >'
    ),
    '\t\t<FixedModification name="Carbamidomethyl (C)" />',
    '\t\t<VariableModification name="Oxidation (M)" />',
    *(
        f'\t\t\t\t<UserParam type="{kind}" name="CometAdapter:1:{name}" value="{value}"/>'
        for kind, name, value in [
            ("string", "in", "test.mzML"),
            ("string", "out", "test.idXML"),
            ("string", "database", "UP000005640_9606_decoy.fasta"),
            ("string", "comet_executable", "comet.exe"),
            ("string", "pin_out", ""),
            ("string", "default_params_file", ""),
            ("float", "precursor_mass_tolerance", "5.0"),
            ("string", "precursor_error_units", "ppm"),
            ("string", "isotope_error", "off"),
            ("float", "fragment_mass_tolerance", "0.50025"),
            ("string", "fragment_error_units", "Da"),
            ("float", "fragment_bin_offset", "0.4"),
            ("string", "instrument", "low_res"),
            ("string", "use_A_ions", "false"),
            ("string", "use_B_ions", "true"),
            ("string", "use_C_ions", "false"),
            ("string", "use_X_ions", "false"),
            ("string", "use_Y_ions", "true"),
            ("string", "use_Z_ions", "false"),
            ("string", "use_NL_ions", "false"),
            ("string", "enzyme", "unspecific cleavage"),
            ("string", "second_enzyme", ""),
            ("string", "num_enzyme_termini", "fully"),
            ("int", "missed_cleavages", "1"),
            ("int", "min_peptide_length", "5"),
        ]
    ),
]


async def openmsthirdparty_text(file: File) -> str:
    async with file.open("rb") as fh:
        return bytes(await fh.read()).decode()


async def openmsthirdparty_comet_search() -> tuple[File, CometAdapterResult]:
    """Upstream's setup: the RAW file as ``test.mzML``, searched with Comet against the decoyed human proteome."""
    mzml = await thermorawfileparser(await fixture(RAW), prefix="test")
    decoys = await openms_decoydatabase(await fixture(HUMAN_PROTEOME))
    task = openmsthirdparty_cometadapter_cmd.override(**FOUR_CPUS)
    results = await task(mzml=mzml, fasta=decoys, prefix="test", args=COMET_ARGS)
    return mzml, await comet_result(results, "test")


@env.task
async def test_openmsthirdparty_cometadapter() -> None:
    # upstream case: openmsthirdparty/cometadapter "proteomics - comet"
    result = (await openmsthirdparty_comet_search())[1]
    assert result.idxml.path.endswith("/test.idXML"), result.idxml.path
    assert result.pin is None, "Comet wrote Percolator input without -pin_out"
    lines = (await openmsthirdparty_text(result.idxml)).splitlines()
    assert lines[:31] == EXPECTED_HEAD, f"idXML lines 1-31 differ: {lines[:31]}"
    assert any("ProteinHit" in line for line in lines), "no ProteinHit"
    assert any("PeptideHit" in line for line in lines), "no PeptideHit"


tests = [test_openmsthirdparty_cometadapter]
