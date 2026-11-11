"""Tests for flyte_bio.modules.freyja.

The cases are upstream's freyja/update, freyja/variants, freyja/demix and
freyja/boot tests (modules commit 9a1cc9f, freyja 2.0.6), with upstream's
checks: update's file names and that its pathogen config mentions MPX (it
downloads live data, as upstream); variants' md5s; demix's output containing
``summarized``; boot's outputs containing ``0.025,``. demix and boot run
with upstream's test ``--depthcutoff 1``.

Upstream's demix and boot cases take their barcodes and lineage data from a
live ``freyja update``. Here they use a pinned copy instead: the
2026-10-06 barcodes, hierarchy and metadata from the Freyja-data archive,
the same files bundled in the freyja 2.0.6 image. Without ``--lineageyml``
demix falls back on the image's bundled hierarchy, which newer live barcodes
outrun: with the 2026-11-09 barcodes the "no meta/topology" case fails
(freyja's "Lineage hierarchy file is likely behind the selected barcode
file"), so live data makes that case depend on the day it runs.
"""

from flyte_bio.modules.freyja import (
    FreyjaVariantsResult,
    freyja_boot,
    freyja_demix,
    freyja_update,
    freyja_variants,
)
from tests.framework import assert_md5, env, fixture

# Freyja-data's archive of the 2026-10-06 data, as bundled in the freyja 2.0.6 image.
FREYJA_DATA = "https://raw.githubusercontent.com/andersen-lab/Freyja-data/1a14468ec194e67e2e0072ec8d0dfdf72e07b8f0/"
PINNED_BARCODES = FREYJA_DATA + "history_barcodes/usher_barcodes10_06_2026-06-39.feather"
PINNED_TOPOLOGY = FREYJA_DATA + "history_lineage_hierarchy/lineages10_06_2026-06-39.yml"
PINNED_META = FREYJA_DATA + "history_metadata/curated_lineages10_06_2026-06-39.json"


async def text(file) -> str:
    async with file.open("rb") as fh:
        return bytes(await fh.read()).decode()


async def freyja_test_variants() -> FreyjaVariantsResult:
    return await freyja_variants(
        await fixture("genomics/sarscov2/illumina/bam/test.paired_end.sorted.bam"),
        await fixture("genomics/sarscov2/genome/genome.fasta"),
        prefix="test",
    )


@env.task
async def test_freyja_update() -> None:
    # upstream case: freyja/update "sarscov2 - value db name "
    result = await freyja_update()
    assert "MPX" in await text(result.config), "pathogen config lacks MPX"
    names = (result.barcodes.name, result.lineages_meta.name if result.lineages_meta else None)
    assert names == ("usher_barcodes.feather", "curated_lineages.json"), names
    assert result.lineages_topology.name == "lineages.yml", result.lineages_topology.name


@env.task
async def test_freyja_update_mpx() -> None:
    # upstream case: freyja/update "MPX - value db name "
    result = await freyja_update(args="--pathogen MPX")
    assert "MPX" in await text(result.config), "pathogen config lacks MPX"
    names = (result.barcodes.name, result.lineages_topology.name)
    assert names == ("MPX_barcodes.csv", "MPX_lineages.yml"), names


@env.task
async def test_freyja_variants() -> None:
    # upstream case: freyja/variants "sarscov2 - illumina - test_paired_end_sorted - bam"
    result = await freyja_test_variants()
    await assert_md5(result.variants, "966450bae4d9abae278572927b821983", label="test.variants.tsv")
    await assert_md5(result.depths, "27f79b28a365a8af915895b484d1153e", label="test.depth.tsv")


@env.task
async def test_freyja_demix() -> None:
    # upstream case: freyja/demix "sarscov2 - illumina - test_paired_end_sorted - bam"
    variants = await freyja_test_variants()
    demix = await freyja_demix(
        variants.variants,
        variants.depths,
        await fixture(PINNED_BARCODES),
        lineages_meta=await fixture(PINNED_META),
        lineages_topology=await fixture(PINNED_TOPOLOGY),
        prefix="test",
        args="--depthcutoff 1",
    )
    assert "summarized" in await text(demix), "demix output lacks 'summarized'"


@env.task
async def test_freyja_demix_no_lineages() -> None:
    # upstream case: freyja/demix "... - no meta/topology"
    variants = await freyja_test_variants()
    demix = await freyja_demix(
        variants.variants, variants.depths, await fixture(PINNED_BARCODES), prefix="test", args="--depthcutoff 1"
    )
    assert "summarized" in await text(demix), "demix output lacks 'summarized'"


@env.task
async def test_freyja_boot() -> None:
    # upstream case: freyja/boot "sarscov2 - illumina - test_paired_end_sorted - bam"
    variants = await freyja_test_variants()
    boot = await freyja_boot(
        variants.variants,
        variants.depths,
        3,
        await fixture(PINNED_BARCODES),
        lineages_meta=await fixture(PINNED_META),
        lineages_topology=await fixture(PINNED_TOPOLOGY),
        prefix="test",
        args="--depthcutoff 1",
    )
    assert "0.025," in await text(boot.summarized), "summarized lacks '0.025,'"
    assert "0.025," in await text(boot.lineages), "lineages lacks '0.025,'"


tests = [
    test_freyja_update,
    test_freyja_update_mpx,
    test_freyja_variants,
    test_freyja_demix,
    test_freyja_demix_no_lineages,
    test_freyja_boot,
]
