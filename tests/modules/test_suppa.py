"""Tests for flyte_bio.modules.suppa.

The cases and expected values are upstream's suppa/generateevents,
suppa/psiperevent, suppa/psiperisoform, suppa/diffsplice and
suppa/clusterevents tests (modules commit 9a1cc9f, SUPPA 2.4), chained as
upstream's setup steps chain them: events from the chr21 GTF, PSIs of upstream's
mock expression tables (written here byte for byte as upstream's tests
write them), then diffSplice and clusterEvents.

As upstream, generateEvents is checked by the event files' names and the event
GTFs' md5s, and the empirical diffSplice runs on ``.ioe`` events (which
sample at random) by their outputs' names. Upstream's "poolgenes"
generateEvents cases pass ``pool_genes = false``; they are ported as written.
Upstream's psiPerEvent "A5" case takes the first ``SS`` event file, which is
``A3``'s, and its "A3" case the second; ported as written.
"""

import tempfile
from pathlib import Path

from flyte.io import File

from flyte_bio.modules.suppa import (
    SuppaDiffSplice,
    SuppaEvents,
    suppa_clusterevents,
    suppa_diffsplice,
    suppa_generateevents,
    suppa_psiperevent,
    suppa_psiperisoform,
)
from tests.framework import assert_md5, env, fixture

GTF = "genomics/homo_sapiens/genome/chr21/sequence/chr21_gencode.gtf"
ALL_EVENTS = "SE SS MX RI FL"
EVENT_TYPES = ["A3", "A5", "AF", "AL", "MX", "RI", "SE"]
STRICT_GTFS = {
    "A3": "b913fbfe02db79b24383bccf3dc03f72",
    "A5": "1d6aebfaf31d2cd46e6f7525bd20e5a5",
    "AF": "9809083cd1a3a63a4e525be37fc2bd3d",
    "AL": "d35ff041853a3ddc182ed8ec5e2dc70a",
    "MX": "184e9fdac966818f8c0892974cce3317",
    "RI": "ee08abdbba6714252fd24b752317eed6",
    "SE": "6e4842f5c02f1321c83e758f49f3f6f8",
}
VARIABLE_GTFS = {
    **STRICT_GTFS,
    "A3": "fd7c24dac526acb2a9639438d27fac01",
    "A5": "149f0fa75a37efff30e06727c0c83f54",
    "RI": "cfb22e58c83a1b07170c2c179d9e0121",
}
# Upstream's psiPerEvent / psiPerIsoform mock expression table.
MOCK_EXPRESSION = "sample1\nENST00000625012.1\t10.5\nENSE00003760353.1\t95.1\nENSE00003754925.1\t0.0"
# The chr21 transcripts of upstream's diffSplice / clusterEvents mock expression tables.
TRANSCRIPTS = (
    "ENST00000614229.4 ENST00000400566.6 ENST00000658863.1 ENST00000635225.1 ENST00000635621.1 ENST00000449602.2 "
    "ENST00000654045.1 ENST00000669763.1 ENST00000660967.1 ENST00000339775.10 ENST00000284881.9 ENST00000671660.1 "
    "ENST00000656151.1 ENST00000656404.1 ENST00000666665.1 ENST00000400094.5 ENST00000348990.9 ENST00000346798.8 "
    "ENST00000357903.7 ENST00000439274.6 ENST00000440126.7 ENST00000303775.10 ENST00000286788.9 ENST00000626972.2 "
    "ENST00000399913.5 ENST00000399907.5 ENST00000622113.4 ENST00000618007.4 ENST00000285681.6 ENST00000400183.6 "
    "ENST00000400166.5 ENST00000400165.5 ENST00000400128.5 ENST00000400127.5 ENST00000299295.7 ENST00000543733.5 "
    "ENST00000400135.5 ENST00000480456.6 ENST00000400532.5 ENST00000596385.5 ENST00000660773.1 ENST00000399935.6 "
    "ENST00000399934.5 ENST00000339024.8 ENST00000344577.6 ENST00000536861.5 ENST00000619120.4 ENST00000634644.1 "
    "ENST00000635525.2 ENST00000348354.7 ENST00000400558.7 ENST00000663293.1 ENST00000654371.1 ENST00000669643.1 "
    "ENST00000284971.7 ENST00000354192.7 ENST00000359726.7 ENST00000351429.7 ENST00000540844.5 ENST00000389125.7 "
    "ENST00000427445.6 ENST00000285679.10 ENST00000356275.10 ENST00000400131.5 ENST00000312957.9 ENST00000600590.1 "
    "ENST00000668267.1 ENST00000399947.6 ENST00000341618.8"
).split()


async def local_file(name: str, content: str) -> File:
    path = Path(tempfile.mkdtemp()) / name
    path.write_text(content)
    return await File.from_local(str(path))


def tenths(n: int, value_in_tenths: int) -> str:
    """``n * value`` printed as Groovy prints a one-decimal BigDecimal product (``(i + 1) * 10.5``)."""
    total = n * value_in_tenths
    return f"{total // 10}.{total % 10}"


async def diffsplice_tables() -> tuple[File, File]:
    """Upstream's diffSplice mock expression tables of conditions 1 and 2."""
    rows_1 = ["sample1\tsample2"]
    rows_2 = ["sample1\tsample2"]
    for i, tx in enumerate(TRANSCRIPTS):
        rows_1.append(f"{tx}\t{tenths(i + 1, 105)}\t{tenths(i + 1, 21)}")
        rows_2.append(f"{tx}\t{tenths(i + 1, 15)}\t{tenths(i + 1, 153)}")
    return await local_file("mock_tpm_1.txt", "\n".join(rows_1)), await local_file("mock_tpm_2.txt", "\n".join(rows_2))


async def clusterevents_tables() -> tuple[File, File]:
    """Upstream's clusterEvents mock tables: high/low expression alternating, inverted in condition 2."""
    rows_1 = ["sample1\tsample2"]
    rows_2 = ["sample3\tsample4"]
    for i, tx in enumerate(TRANSCRIPTS):
        high = i % 2 == 0
        value_1, value_2 = ("100.0", "5.0") if high else ("5.0", "100.0")
        rows_1.append(f"{tx}\t{value_1}\t{value_1}")
        rows_2.append(f"{tx}\t{value_2}\t{value_2}")
    return await local_file("mock_tpm_1.txt", "\n".join(rows_1)), await local_file("mock_tpm_2.txt", "\n".join(rows_2))


async def events(format: str, event_type: str = "", boundary: str = "") -> SuppaEvents:
    return await suppa_generateevents(
        await fixture(GTF), format=format, prefix="human", event_type=event_type, boundary=boundary
    )


async def check_events(result: SuppaEvents, suffix: str, gtfs: dict[str, str]) -> None:
    names = [f.name for f in result.events]
    expected = [f"human_{event}_{suffix}.ioe" for event in EVENT_TYPES]
    assert names == expected, f"event files {names}, expected {expected}"
    gtf_names = [f.name for f in result.gtfs]
    expected = [f"human_{event}_{suffix}.gtf" for event in EVENT_TYPES]
    assert gtf_names == expected, f"event GTFs {gtf_names}, expected {expected}"
    for event, file in zip(EVENT_TYPES, result.gtfs, strict=True):
        await assert_md5(file, gtfs[event], label=file.name)


async def check_ioi(result: SuppaEvents) -> None:
    assert [f.name for f in result.events] == ["human.ioi"], f"event files {[f.name for f in result.events]}"


@env.task
async def test_suppa_generateevents_ioe_strict() -> None:
    # upstream case: suppa/generateevents "human - ioe strict - gtf"
    await check_events(await events("ioe", ALL_EVENTS, "S"), "strict", STRICT_GTFS)


@env.task
async def test_suppa_generateevents_ioe_variable() -> None:
    # upstream case: suppa/generateevents "human - ioe variable - gtf"
    await check_events(await events("ioe", ALL_EVENTS, "V"), "variable_10", VARIABLE_GTFS)


@env.task
async def test_suppa_generateevents_ioi() -> None:
    # upstream case: suppa/generateevents "human - ioi - gtf"
    await check_ioi(await events("ioi"))


@env.task
async def test_suppa_generateevents_ioe_strict_poolgenes() -> None:
    # upstream case: suppa/generateevents "human - ioe strict poolgenes - gtf" (upstream passes pool_genes = false)
    await check_events(await events("ioe", ALL_EVENTS, "S"), "strict", STRICT_GTFS)


@env.task
async def test_suppa_generateevents_ioe_variable_poolgenes() -> None:
    # upstream case: suppa/generateevents "human - ioe variable poolgenes - gtf" (upstream passes pool_genes = false)
    await check_events(await events("ioe", ALL_EVENTS, "V"), "variable_10", VARIABLE_GTFS)


@env.task
async def test_suppa_generateevents_ioi_poolgenes() -> None:
    # upstream case: suppa/generateevents "human - ioi poolgenes - gtf" (upstream passes pool_genes = false)
    await check_ioi(await events("ioi"))


async def check_psi(event_type: str, index: int, expected: str) -> None:
    result = await events("ioe", event_type, "S")
    expression = await local_file("mock_expression_matrix.txt", MOCK_EXPRESSION)
    psi = await suppa_psiperevent(expression, result.events[index], prefix="sample1")
    await assert_md5(psi, expected, label=f"{result.events[index].name} psi")


@env.task
async def test_suppa_psiperevent_se() -> None:
    # upstream case: suppa/psiperevent "human - skipping exon - ioe"
    await check_psi("SE", 0, "eebc549dee0797dd5bf62d65d1296dae")


@env.task
async def test_suppa_psiperevent_a5() -> None:
    # upstream case: suppa/psiperevent "human - alternative splice site A5 - ioe" (the first SS file, A3's)
    await check_psi("SS", 0, "04daef324165f4f77f872124844ff657")


@env.task
async def test_suppa_psiperevent_a3() -> None:
    # upstream case: suppa/psiperevent "human - alternative splice site A3 - ioe" (the second SS file, A5's)
    await check_psi("SS", 1, "237d9e5dc3fb79e8149a5527ed51a3c1")


@env.task
async def test_suppa_psiperevent_mx() -> None:
    # upstream case: suppa/psiperevent "human - mutually exclusive exons - ioe"
    await check_psi("MX", 0, "9fb5a764439d9b907bf50944248bb8cd")


@env.task
async def test_suppa_psiperevent_ri() -> None:
    # upstream case: suppa/psiperevent "human - retained intron - ioe"
    await check_psi("RI", 0, "85917088c83d53c6c578bfdee28539e6")


@env.task
async def test_suppa_psiperevent_af() -> None:
    # upstream case: suppa/psiperevent "human - alternative first exon - ioe"
    await check_psi("FL", 0, "40a096c2d820965a912ec3907f94bb0c")


@env.task
async def test_suppa_psiperevent_al() -> None:
    # upstream case: suppa/psiperevent "human - alternative last exons - ioe"
    await check_psi("FL", 1, "9f665c73ab259b5d475bc10a5e5f467b")


@env.task
async def test_suppa_psiperisoform() -> None:
    # upstream case: suppa/psiperisoform "human - tpm"
    expression = await local_file("mock_expression_matrix.txt", MOCK_EXPRESSION)
    psi = await suppa_psiperisoform(expression, await fixture(GTF), prefix="sample1")
    assert psi.name == "sample1.psi_isoform.psi", f"unexpected name {psi.name}"
    await assert_md5(psi, "997759ca63af7047f47d3cb2355dc97f", label="psi")


async def diffsplice(
    tables: tuple[File, File], format: str, prefix: str, method: str, **settings: object
) -> SuppaDiffSplice:
    """Upstream's setup chain: events (SE for ioe), each condition's PSIs, then diffSplice."""
    result = await events(format, "SE", "S") if format == "ioe" else await events("ioi")
    (event_file,) = result.events
    tpm_1, tpm_2 = tables
    psi_1 = await suppa_psiperevent(tpm_1, event_file, prefix="samples1")
    psi_2 = await suppa_psiperevent(tpm_2, event_file, prefix="samples2")
    return await suppa_diffsplice(
        event_file,
        psi_1,
        tpm_1,
        psi_2,
        tpm_2,
        method=method,
        prefix=prefix,
        local_area=1000,
        multi_testing_correction_alpha=0.05,
        **settings,
    )


def check_names(result: SuppaDiffSplice, prefix: str, tpm_events: bool = False) -> None:
    assert result.dpsi.name == f"{prefix}.dpsi" and result.psivec.name == f"{prefix}.psivec"
    if tpm_events:
        assert result.tpm_events is not None, "no tpm events table"
        assert result.tpm_events.name == f"{prefix}_avglogtpm.tab", f"unexpected {result.tpm_events.name}"
    else:
        assert result.tpm_events is None, "unexpected tpm events table"


@env.task
async def test_suppa_diffsplice_ioe() -> None:
    # upstream case: suppa/diffsplice "human - ioe"
    result = await diffsplice(await diffsplice_tables(), "ioe", "test_SE", "empirical")
    check_names(result, "test_SE")


@env.task
async def test_suppa_diffsplice_ioe_lower_bound() -> None:
    # upstream case: suppa/diffsplice "human - ioe - custom lower bound delta psi"
    result = await diffsplice(await diffsplice_tables(), "ioe", "test_SE", "empirical", lower_bound_delta_psi=2)
    check_names(result, "test_SE")


@env.task
async def test_suppa_diffsplice_ioe_save_tpm_events() -> None:
    # upstream case: suppa/diffsplice "human - ioe - save tpm events"
    result = await diffsplice(await diffsplice_tables(), "ioe", "test_SE", "empirical", save_tpm_events=True)
    check_names(result, "test_SE", tpm_events=True)


@env.task
async def test_suppa_diffsplice_ioi() -> None:
    # upstream case: suppa/diffsplice "human - ioi"
    result = await diffsplice(await diffsplice_tables(), "ioi", "test", "empirical")
    check_names(result, "test")
    await assert_md5(result.dpsi, "8419e12b53abb23ea61a2c2ce56cb338", label="dpsi")
    await assert_md5(result.psivec, "9e57ba7f3062b28cd4dce83ab66d15e3", label="psivec")


async def clustered(**settings: object) -> File:
    result = await diffsplice(await clusterevents_tables(), "ioe", "test_SE", "classical")
    return await suppa_clusterevents(
        result.dpsi,
        result.psivec,
        groups="1-2,3-4",
        prefix="test_SE",
        significance_threshold=1.0,
        dpsi_threshold=0.0,
        maximum_distance=0.0,
        metric="euclidean",
        **settings,
    )


@env.task
async def test_suppa_clusterevents_dbscan() -> None:
    # upstream case: suppa/clusterevents "human - dbscan"
    clustvec = await clustered(min_number_events=1, clustering_method="DBSCAN")
    await assert_md5(clustvec, "b4ca48e5306183aaecc9b890fe87e13c", label="clustvec")


@env.task
async def test_suppa_clusterevents_optics() -> None:
    # upstream case: suppa/clusterevents "human - optics"
    clustvec = await clustered(separation=0.9, min_number_events=2, clustering_method="OPTICS")
    await assert_md5(clustvec, "dea9bf23050c80a91a5abc357168ec9a", label="clustvec")


tests = [
    test_suppa_generateevents_ioe_strict,
    test_suppa_generateevents_ioe_variable,
    test_suppa_generateevents_ioi,
    test_suppa_generateevents_ioe_strict_poolgenes,
    test_suppa_generateevents_ioe_variable_poolgenes,
    test_suppa_generateevents_ioi_poolgenes,
    test_suppa_psiperevent_se,
    test_suppa_psiperevent_a5,
    test_suppa_psiperevent_a3,
    test_suppa_psiperevent_mx,
    test_suppa_psiperevent_ri,
    test_suppa_psiperevent_af,
    test_suppa_psiperevent_al,
    test_suppa_psiperisoform,
    test_suppa_diffsplice_ioe,
    test_suppa_diffsplice_ioe_lower_bound,
    test_suppa_diffsplice_ioe_save_tpm_events,
    test_suppa_diffsplice_ioi,
    test_suppa_clusterevents_dbscan,
    test_suppa_clusterevents_optics,
]
