"""Tests for flyte_bio.modules.scimap.

The case and expected values are upstream's scimap/mcmicro test (modules
commit 9a1cc9f, scimap 2.1.3). As upstream, the clustering outputs are checked
by name and by the table's header, since the clusters themselves aren't
reproducible.
"""

from flyte_bio.modules.scimap import scimap_mcmicro
from tests.framework import env, fixture

HEADER = (
    "CellID,Na/K ATPase,CD3,CD45RO,Antigen Ki67,Pan-cytokeratin,Aortic smooth muscle actin,X_centroid,Y_centroid,"
    "Area,MajorAxisLength,MinorAxisLength,Eccentricity,Solidity,Extent,Orientation,imageid,kmeans,leiden,"
    "phenograph,spatial_kmeans"
)


@env.task
async def test_scimap_mcmicro() -> None:
    # upstream case: "test-scimap-mcmicro"
    result = await scimap_mcmicro(await fixture("imaging/downstream/cycif_tonsil_cell.csv"))
    assert result.csv is not None and result.h5ad is not None, "scimap wrote no .h5ad.csv or .h5ad"
    assert result.csv.name == "cycif_tonsil_cell.h5ad.csv", f"got {result.csv.name}"
    assert result.h5ad.name == "cycif_tonsil_cell.h5ad", f"got {result.h5ad.name}"
    async with result.csv.open("rb") as fh:
        first = bytes(await fh.read(1 << 16)).decode(errors="replace").splitlines()[0]
    assert first == HEADER, f"unexpected header {first!r}"


tests = [test_scimap_mcmicro]
