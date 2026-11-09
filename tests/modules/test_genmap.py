"""Tests for flyte_bio.modules.genmap.

The case and expected md5s are upstream's genmap/index test (modules commit
9a1cc9f, GenMap 1.3.0): every file of the index of the SARS-CoV-2 genome.
"""

from flyte_bio.modules.genmap import genmap_index
from tests.framework import assert_md5, env, fixture

INDEX_MD5S = {
    "index.ids.concat": "da6caa25f62c5407ccdfbcce1fa92408",
    "index.ids.limits": "f82636c5da188aec131d3a809473eff1",
    "index.info.concat": "8ba5273aa9e58722bf45b9cc39fc6bfe",
    "index.info.limits": "3522f2811f4ddf04598809fc84a1459e",
    "index.lf.drp": "dd85d6a23af2c7adf2695658e3056c08",
    "index.lf.drp.sbl": "f1d3ff8443297732862df21dc4e57262",
    "index.lf.drs": "93b885adfe0da089cdf634904fd59f71",
    "index.lf.drv": "e06b605496bd91b32afa3c4f56d934ac",
    "index.lf.drv.sbl": "8dd6bb7329a71449b0a1b292b5999164",
    "index.lf.pst": "e8daba34298e99e42942435286f9b3f0",
    "index.rev.lf.drp": "5d9107e3aeec0721553dd661d4365fef",
    "index.rev.lf.drp.sbl": "f1d3ff8443297732862df21dc4e57262",
    "index.rev.lf.drs": "93b885adfe0da089cdf634904fd59f71",
    "index.rev.lf.drv": "df7e795edc0a034577a9d2599fe8cfeb",
    "index.rev.lf.drv.sbl": "8dd6bb7329a71449b0a1b292b5999164",
    "index.rev.lf.pst": "e8daba34298e99e42942435286f9b3f0",
    "index.sa.ind": "e21e5c7ce887cc8e3d0fa44ab1019cab",
    "index.sa.len": "5dfc20cfe8ed9892451461a8d402f51c",
    "index.sa.val": "400ee7f2fe93b2000ae3a5da5e509730",
    "index.txt.concat": "b4303962e0c176107945f3405370e6ae",
    "index.txt.limits": "4480a068db603e4c9a27bc4fa9ceaf14",
}


@env.task
async def test_genmap_index() -> None:
    index = await genmap_index(await fixture("genomics/sarscov2/genome/genome.fasta"), prefix="test")
    files = {file.path.rstrip("/").rsplit("/", 1)[-1]: file for file in await index.list_files()}
    assert sorted(files) == sorted(INDEX_MD5S), f"index files: {sorted(files)}"
    for name, md5 in INDEX_MD5S.items():
        await assert_md5(files[name], md5, label=name)


tests = [test_genmap_index]
