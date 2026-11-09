"""Tests for flyte_bio.modules.mcquant.

The case and expected md5 are upstream's mcquant test (modules commit 9a1cc9f,
quantification 1.5.4). Upstream ships a test config adding
``--intensity_props intensity_sum``, but its test doesn't load it, so the
snapshot (and this test) use mcquant's defaults. As upstream, the md5 is over
the table's lines joined by newlines, after rounding every number with seven
or more decimals to six (Java's ``%.6f``: half-up on the number's shortest
decimal form), so trailing floating-point digits don't matter.
"""

import hashlib
import re
from decimal import ROUND_HALF_UP, Decimal

from flyte.io import File

from flyte_bio.modules.mcquant import mcquant
from tests.framework import env, fixture

DATA = "imaging/quantification/"
LONG_DECIMAL = re.compile(r"\d+\.\d{7,}")


def mcquant_round(match: re.Match) -> str:
    value = Decimal(repr(float(match.group())))
    return str(value.quantize(Decimal("0.000001"), rounding=ROUND_HALF_UP))


async def rounded_md5(table: File) -> str:
    async with table.open("rb") as fh:
        lines = bytes(await fh.read()).decode().splitlines()
    text = "\n".join(LONG_DECIMAL.sub(mcquant_round, line) for line in lines)
    return hashlib.md5(text.encode()).hexdigest()


@env.task
async def test_mcquant_image_mask_markers() -> None:
    # upstream case: "imaging - image,mask,marker"
    table = await mcquant(
        await fixture(DATA + "cycif_tonsil_registered.ome.tif"),
        await fixture(DATA + "cell.ome.tif"),
        await fixture(DATA + "markers.csv"),
    )
    assert table.name == "cycif_tonsil_registered_cell.csv", f"got {table.name}"
    md5 = await rounded_md5(table)
    assert md5 == "e1559b2be1611815a3335bf83b48d22d", f"rounded table md5 {md5}"


tests = [test_mcquant_image_mask_markers]
