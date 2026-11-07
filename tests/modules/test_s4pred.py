"""Tests for flyte_bio.modules.s4pred.

The case and expected md5s are upstream's s4pred/runmodel test (modules commit
9a1cc9f, S4PRED 1.2.1): one ``.ss2`` per protein of the sarscov2 proteome.
"""

import hashlib

from flyte_bio.modules.s4pred import s4pred_runmodel
from tests.framework import env, fixture

PREDICTIONS = {
    "ENSSASP00005000002.1.ss2": "ca1ecb796927219ab1bd2408e17d94d6",
    "ENSSASP00005000003.1.ss2": "ad37a8ec4cef27c91c9958bcb652903f",
    "ENSSASP00005000004.1.ss2": "dae6e305856d12c20c97b86492c937cd",
    "ENSSASP00005000005.1.ss2": "f3e72765229ac89063633f6d4e40c4fe",
    "ENSSASP00005000006.1.ss2": "731586e1e6176ada900c181edb17e3c9",
    "ENSSASP00005000007.1.ss2": "2f7459999558834d339b781c96da1482",
    "ENSSASP00005000008.1.ss2": "c620bbff5f073089df9554ca439933eb",
    "ENSSASP00005000009.1.ss2": "19f317c1b3d468ec2e876ecce4d4b60a",
    "ENSSASP00005000010.1.ss2": "020604f0e037a104983e7ba325821c0a",
    "ENSSASP00005000011.1.ss2": "ee9dc69931be3ec7eeed88e42f7b90e7",
    "ENSSASP00005000012.1.ss2": "f08bffd8ce8031fb7bb658f8e2408626",
    "ENSSASP00005000013.1.ss2": "445e88597dc8fc0f40e21d95d148ced1",
}


@env.task
async def test_s4pred_sarscov2_proteome() -> None:
    preds = await s4pred_runmodel(await fixture("genomics/sarscov2/genome/proteome.fasta"))
    actual = {}
    async for file in preds.walk():
        async with file.open("rb") as fh:
            actual[file.name] = hashlib.md5(bytes(await fh.read())).hexdigest()
    assert sorted(actual) == sorted(PREDICTIONS), f"files {sorted(actual)}"
    for name, md5 in PREDICTIONS.items():
        assert actual[name] == md5, f"{name} md5 {actual[name]}, expected {md5}"


tests = [
    test_s4pred_sarscov2_proteome,
]
