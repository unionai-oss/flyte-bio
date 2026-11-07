"""Tests for flyte_bio.modules.cellranger.

The cases and expected md5s are upstream's cellranger/mkvdjref and
cellranger/vdj tests (modules commit 9a1cc9f, Cell Ranger 10.1.0):

- mkvdjref builds a reference (``homo_sapiens_chr2_reference``, upstream's
  name) from Ensembl 103's chromosome 22 and its whole-genome GTF, both
  gunzipped first, and checks ``fasta/regions.fa`` and ``reference.json``;
- vdj runs a subsampled 10x T cell sample against the GRCh38 5.0.0 V(D)J
  reference's ``reference.json``, ``fasta/regions.fa`` and
  ``fasta/supp_regions.fa``, copied into a ``vdj_reference`` folder as
  upstream does (here archived as a ``.tar.gz``), and checks the outputs
  upstream snapshots: ``vdj_contig_info.pb``, ``metrics_summary.csv``,
  ``vdj_reference/reference.json`` and ``web_summary.html``.

They run with the module's default resources (4 CPUs; 8 GiB for mkvdjref,
16 GiB for vdj), a little above upstream's tests' 6 GB.
"""

import asyncio
import shutil
import tarfile
import tempfile
from pathlib import Path

from flyte.io import File

from flyte_bio.modules.cellranger import cellranger_mkvdjref, cellranger_vdj
from flyte_bio.modules.gunzip import gunzip
from tests.framework import assert_md5, env, fixture

ENSEMBL = "http://ftp.ensembl.org/pub/release-103/"
CHR22 = ENSEMBL + "fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.22.fa.gz"
GTF = ENSEMBL + "gtf/homo_sapiens/Homo_sapiens.GRCh38.103.gtf.gz"
CELLRANGER = "genomics/homo_sapiens/10xgenomics/cellranger/"
VDJ_REFERENCE = CELLRANGER + "references/vdj/refdata-cellranger-vdj-GRCh38-alts-ensembl-5.0.0/"
VDJ_REFERENCE_FILES = ("reference.json", "fasta/regions.fa", "fasta/supp_regions.fa")
FASTQS = CELLRANGER + "5k_cmvpos_tcells/fastqs/vdj/subsampled_5k_human_antiCMV_T_TBNK_connect_VDJ_S1_L001_"
SAMPLE = "subsampled_5k_human_antiCMV_T_TBNK_connect_VDJ"


@env.task(cache="auto")
async def cellranger_vdj_test_reference() -> File:
    """Upstream's test V(D)J reference files as ``vdj_reference.tar.gz``, a ``vdj_reference`` folder."""
    root = Path(tempfile.mkdtemp())
    for name in VDJ_REFERENCE_FILES:
        local = await (await fixture(VDJ_REFERENCE + name)).download()
        dest = root / "vdj_reference" / name
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.move(local, dest)
    archive = root / "vdj_reference.tar.gz"
    with tarfile.open(archive, "w:gz") as tar:
        tar.add(root / "vdj_reference", arcname="vdj_reference")
    return await File.from_local(str(archive))


@env.task
async def test_cellranger_mkvdjref() -> None:
    fasta, gtf = await asyncio.gather(gunzip(await fixture(CHR22)), gunzip(await fixture(GTF)))
    reference = await cellranger_mkvdjref("homo_sapiens_chr2_reference", fasta=fasta, gtf=gtf)
    await assert_md5(reference.regions, "1e22fd84ab26c89d3c285140df2daf5f", label="regions.fa")
    await assert_md5(reference.reference_json, "ed1c323eb2f4d2acbc595c0aaee955ae", label="reference.json")


@env.task
async def test_cellranger_vdj() -> None:
    reads = [await fixture(FASTQS + "R1_001.fastq.gz"), await fixture(FASTQS + "R2_001.fastq.gz")]
    result = await cellranger_vdj(reads, await cellranger_vdj_test_reference(), prefix=SAMPLE)
    for path, md5 in (
        ("vdj_contig_info.pb", "e36273c501265c8581df0195bd70eb36"),
        ("metrics_summary.csv", "7994a492bd179334628e65b322c9df15"),
        ("vdj_reference/reference.json", "384f6efabad59f6da5c89b862aee71a8"),
        ("web_summary.html", "6fe0002bd4babae46e5c6560685113ff"),
    ):
        file = File.from_existing_remote(result.outs.path.rstrip("/") + "/" + path)
        await assert_md5(file, md5, label=path)
    assert result.airr_rearrangement is not None, "cellranger vdj wrote no airr_rearrangement.tsv"


tests = [test_cellranger_mkvdjref, test_cellranger_vdj]
