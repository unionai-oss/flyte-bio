"""Tests for flyte_bio.modules.bandage.

The cases are upstream's bandage/image tests (modules commit 9a1cc9f,
Bandage 0.9.0): the PNG's name and lines 4-8 of the SVG, for a plain GFA and
for one bgzipped first with htslib (as upstream's setup does).
"""

from flyte.io import File

from flyte_bio.modules.bandage import bandage_image
from flyte_bio.modules.htslib import htslib_bgziptabix
from tests.framework import env, fixture

GFA = "genomics/sarscov2/illumina/gfa/assembly.gfa"
SVG_LINES = [
    ' xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink"  version="1.2" baseProfile="tiny">',
    "<title>Qt SVG Document</title>",
    "<desc>Generated with Qt</desc>",
    "<defs>",
    "</defs>",
]


async def check_bandage_images(gfa: File) -> None:
    result = await bandage_image(gfa, prefix="B-3106")
    assert result.png.name == "B-3106.png", result.png.name
    async with result.svg.open("rb") as fh:
        lines = bytes(await fh.read()).decode().splitlines()[3:8]
    assert lines == SVG_LINES, f"B-3106.svg lines 4-8: {lines}"


@env.task
async def test_bandage_image() -> None:
    # upstream case: bandage/image "test-bandage-image"
    await check_bandage_images(await fixture(GFA))


@env.task
async def test_bandage_image_gzip() -> None:
    # upstream case: bandage/image "test-bandage-image - gzip"
    gfa_gz, _ = await htslib_bgziptabix(await fixture(GFA), prefix="B-3106", out_ext="gfa", make_index=False)
    await check_bandage_images(gfa_gz)


tests = [test_bandage_image, test_bandage_image_gzip]
