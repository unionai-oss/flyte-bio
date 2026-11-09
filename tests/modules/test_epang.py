"""Tests for flyte_bio.modules.epang.

The cases are upstream's epang/place and epang/split tests (modules commit
9a1cc9f, EPA-ng 0.3.8).

- epang_place: three query sequences placed on the reference tree with the LG
  model (upstream's config passes ``--model <meta.model>``), once with the tree
  as is and once gzipped (built in the test, as upstream builds it). Upstream
  snapshots only the output names, as the placements' likelihoods and branch
  lengths are printed at full precision; the tests also check that the jplace
  places each of the three queries on the reference tree.
- epang_split: upstream's setup adds three sequences to the PF14720 seed
  alignment with MAFFT; the full alignment is then split back into query and
  reference sequences, by md5 (decompressed, as nf-test hashes ``.gz``
  outputs).
"""

import gzip
import json
import tempfile
from pathlib import Path

from flyte.io import File

from flyte_bio.modules.epang import epang_place, epang_split
from flyte_bio.modules.mafft import mafft_align
from tests.framework import assert_gunzipped_md5, env, fixture

EPANG = "https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/epang/"
PHYLOPLACE = "https://raw.githubusercontent.com/nf-core/test-datasets/phyloplace/testdata/"
# The query alignment's headers, as epa-ng names its placements.
QUERIES = {"A0A1B7XJK0_9DELT (A0A1B7XJK0)", "A0A1J4NSI0_9ACTN (A0A1J4NSI0)", "B4F2Z1_PROMH (B4F2Z1)"}


async def file_bytes(file: File) -> bytes:
    async with file.open("rb") as fh:
        return bytes(await fh.read())


async def check_place(reference_tree: File) -> None:
    result = await epang_place(
        query=await fixture(EPANG + "query.alnfaa.gz"),
        reference_alignment=await fixture(EPANG + "reference.alnfaa.gz"),
        reference_tree=reference_tree,
        prefix="test",
        args="--model LG",
    )
    assert result.jplace is not None, "epa-ng wrote no jplace"
    assert result.jplace.name == "test.epa_result.jplace.gz", result.jplace.name
    assert result.log.name == "test.epa_info.log", result.log.name
    jplace = json.loads(gzip.decompress(await file_bytes(result.jplace)))
    placed = {name for placement in jplace["placements"] for name in placement["n"]}
    assert placed == QUERIES, f"placed {sorted(placed)}"
    assert jplace["tree"].startswith("((((((Q0P8L5_CAMJE/235-318:"), jplace["tree"][:40]


@env.task
async def test_epang_place() -> None:
    await check_place(await fixture(EPANG + "reference.newick"))


@env.task
async def test_epang_place_gzipped_reference_tree() -> None:
    tree = await file_bytes(await fixture(EPANG + "reference.newick"))
    path = Path(tempfile.mkdtemp()) / "reference.newick.gz"
    path.write_bytes(gzip.compress(tree))
    await check_place(await File.from_local(str(path)))


@env.task
async def test_epang_split_pf14720() -> None:
    seed = await fixture(PHYLOPLACE + "PF14720_seed.alnfaa")
    full = await mafft_align(seed, add=await fixture(PHYLOPLACE + "PF14720_3_sequences.faa"), prefix="test")
    result = await epang_split(seed, full, prefix="test")
    assert result.query.name == "test.query.fasta.gz", result.query.name
    assert result.reference.name == "test.reference.fasta.gz", result.reference.name
    await assert_gunzipped_md5(result.query, "ff1159ff8a5d2587fb20309a81da3c29", label="test.query.fasta.gz")
    await assert_gunzipped_md5(result.reference, "e7753871101c48c42fc417ce0e9b398a", label="test.reference.fasta.gz")


tests = [test_epang_place, test_epang_place_gzipped_reference_tree, test_epang_split_pf14720]
