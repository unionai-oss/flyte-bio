"""Tests for flyte_bio.modules.xeniumranger.

The cases are upstream's xeniumranger/relabel, xeniumranger/resegment and
xeniumranger/importsegmentation tests (modules commit 9a1cc9f, Xenium Ranger
4.0.1.1):

- relabel and resegment run on a tiny human kidney protein bundle, which
  upstream unzips; here it's repacked as a ``.tar.gz`` in the test pod.
  relabel takes the bundle's own ``gene_panel.json``.
- import-segmentation runs on upstream's test bundle (already a ``.tar.gz``),
  importing the bundle's own ``cells.zarr.zip`` as nuclei with an expansion
  distance of 5, in 2 CPUs and 4 GB as upstream's test config sets.

Xenium Ranger's outputs aren't deterministic, so upstream snapshots none of
their md5s: it checks the output folder isn't empty (and, for resegment, that
its main files exist), and snapshots the tool version (``xeniumranger -V``). The
output checks are made here on the output archive. The version isn't checked:
the module doesn't report it, and ``experiment.xenium``'s
``analysis_sw_version`` keeps the input bundle's version.
"""

import asyncio
import tarfile
import tempfile
import zipfile
from pathlib import Path

import flyte
from flyte.io import File

from flyte_bio.modules.xeniumranger import (
    xeniumranger_importsegmentation_cmd,
    xeniumranger_outs,
    xeniumranger_relabel,
    xeniumranger_resegment,
)
from tests.framework import env, fixture

KIDNEY = (
    "https://raw.githubusercontent.com/nf-core/test-datasets/spatialxe/Xenium_V1_Protein_Human_Kidney_tiny_outs.zip"
)
BUNDLE = "https://raw.githubusercontent.com/nf-core/test-datasets/spatialaxe/xenium_bundle.tar.gz"
RESEGMENT_FILES = [
    "analysis_summary.html",
    "cells.csv.gz",
    "cells.parquet",
    "cells.zarr.zip",
    "transcripts.parquet",
    "transcripts.zarr.zip",
    "analysis.zarr.zip",
    "cell_feature_matrix.zarr.zip",
]


def xeniumranger_repack_zip(zip_path: str, name: str) -> tuple[str, str]:
    """The zip at ``zip_path`` repacked as ``<name>.tar.gz`` (files at its top level), and its gene panel."""
    root = Path(tempfile.mkdtemp())
    with zipfile.ZipFile(zip_path) as zf:
        zf.extractall(root / "files")
    archive = root / f"{name}.tar.gz"
    with tarfile.open(archive, "w:gz") as tar:
        for path in sorted((root / "files").iterdir()):
            tar.add(path, arcname=path.name)
    return str(archive), str(root / "files" / "gene_panel.json")


async def xeniumranger_kidney_bundle() -> tuple[File, File]:
    zipped = await fixture(KIDNEY)
    archive, panel = await asyncio.to_thread(xeniumranger_repack_zip, str(await zipped.download()), "kidney_bundle")
    return await File.from_local(archive), await File.from_local(panel)


def xeniumranger_extract_member(archive: str, member: str) -> str:
    """Extract the member of ``archive`` whose path ends in ``member``."""
    root = Path(tempfile.mkdtemp())
    with tarfile.open(archive) as tar:
        found = next(m for m in tar.getmembers() if m.name == member or m.name.endswith("/" + member))
        tar.extract(found, root, filter="data")
    return str(root / found.name)


def xeniumranger_outs_names(archive: str, prefix: str) -> list[str]:
    """The file names under ``<prefix>/`` in ``archive``."""
    names: list[str] = []
    with tarfile.open(archive) as tar:
        for member in tar.getmembers():
            if member.isfile():
                assert member.name.startswith(prefix + "/"), f"{member.name} is outside {prefix}/"
                names.append(member.name.removeprefix(prefix + "/"))
    return names


async def check_outs(archive: File, prefix: str, expected: list[str]) -> None:
    names = await asyncio.to_thread(xeniumranger_outs_names, str(await archive.download()), prefix)
    assert names, f"{prefix}: empty outs"
    missing = [name for name in expected if name not in names]
    assert not missing, f"{prefix}: outs lack {missing}"


@env.task
async def test_xeniumranger_relabel() -> None:
    # upstream case: "xeniumranger relabel"
    bundle, panel = await xeniumranger_kidney_bundle()
    outs = await xeniumranger_relabel(bundle, panel, prefix="test_xeniumranger_relabel")
    await check_outs(outs, "test_xeniumranger_relabel", [])


@env.task
async def test_xeniumranger_resegment() -> None:
    # upstream case: "xeniumranger resegment"
    bundle, _ = await xeniumranger_kidney_bundle()
    outs = await xeniumranger_resegment(bundle, prefix="test_xeniumranger_resegment")
    await check_outs(outs, "test_xeniumranger_resegment", RESEGMENT_FILES)


@env.task
async def test_xeniumranger_importsegmentation_nuclei() -> None:
    # upstream case: "xeniumranger importsegmentation nuclei", in its test config's 2 CPUs and 4 GB
    bundle = await fixture(BUNDLE)
    nuclei = await asyncio.to_thread(xeniumranger_extract_member, str(await bundle.download()), "cells.zarr.zip")
    prefix = "test_xeniumranger_importsegmentation"
    task = xeniumranger_importsegmentation_cmd.override(
        resources=flyte.Resources(cpu=2, memory="4Gi", disk="50Gi"),
        cache=flyte.Cache(behavior="auto", salt="cpu=2,mem=4Gi"),
    )
    results = await task(
        bundle=bundle,
        transcript_assignment=None,
        viz_polygons=None,
        nuclei=await File.from_local(nuclei),
        cells=None,
        coordinate_transform=None,
        expansion_distance="5",
        prefix=prefix,
    )
    await check_outs(await xeniumranger_outs(results, prefix, "import-segmentation"), prefix, [])


tests = [test_xeniumranger_relabel, test_xeniumranger_resegment, test_xeniumranger_importsegmentation_nuclei]
