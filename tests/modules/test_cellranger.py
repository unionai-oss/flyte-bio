"""Tests for flyte_bio.modules.cellranger.

The cases and expected md5s are upstream's cellranger tests (modules commit
9a1cc9f; Cell Ranger 10.1.0, 9.0.1 for mkfastq). nf-test hashes ``.gz``
files decompressed, so these do too.

- mkgtf filters the test genome's GTF with upstream's test biotypes.
- mkref builds ``homo_sapiens_reference`` from the test genome and its
  (unfiltered) GTF and checks the files upstream snapshots inside it.
- count runs a subsampled 10x T cell gene-expression sample (as
  ``test_10x``, ``--chemistry SC3Pv3 --create-bam true``) against a reference
  built from the mkgtf-filtered GTF, and checks every ``outs`` file upstream
  snapshots. For the BAM upstream snapshots nft-bam's ``getSamLinesMD5``,
  which hashes htsjdk's SAM text and ``samtools view`` doesn't reproduce,
  so it gets content checks instead (:func:`check_count_bam`).
- multi runs upstream's cases (10k PBMC 5' GEX with antibody, V(D)J-T and
  V(D)J-B; the 5k CMV+ T cells with V(D)J and antibody; the A549 CRISPR
  screen on its own chromosome 1 slice; the GEM-X Flex kidney singleplex)
  and checks the outputs and config lines upstream asserts. Several
  ``outs`` files share a name (one per sample or library), so they are
  compared as multisets of name and md5. In the Flex case the raw matrices
  (``raw_feature_bc_matrix.h5``, ``sample_raw_feature_bc_matrix.h5``) differ
  from upstream's md5s while every other output matches, so they are checked
  to hold every filtered cell with the filtered counts. Upstream's two
  CellPlex (CMO multiplexing) cases can't be reproduced here (see
  ``CMO_NOT_REPRODUCED``): their config is checked and they report
  NotTested. Upstream's stub cases check only the generated config; they
  are ported as checks of :func:`cellranger_multi_config`, which builds it.
- mkfastq demultiplexes 10x's tiny BCL run (tile 1101) with its simple and
  Illumina sample sheets and checks the FASTQs and ``IndexMetricsOut.bin``.
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
16 GiB for the others), a little above upstream's tests' 6 GB.
"""

import asyncio
import gzip
import hashlib
import shutil
import tarfile
import tempfile
from collections import Counter
from pathlib import Path

from flyte.io import File

from flyte_bio.modules.cellranger import (
    MultiLibrary,
    cellranger_count,
    cellranger_mkfastq,
    cellranger_mkgtf,
    cellranger_mkref,
    cellranger_mkvdjref,
    cellranger_multi,
    cellranger_multi_config,
    cellranger_vdj,
    multi_reference_paths,
)
from flyte_bio.modules.gunzip import gunzip
from flyte_bio.modules.samtools import samtools_view
from tests.framework import NotTested, assert_md5, env, fixture

ENSEMBL = "http://ftp.ensembl.org/pub/release-103/"
CHR22 = ENSEMBL + "fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.22.fa.gz"
GTF = ENSEMBL + "gtf/homo_sapiens/Homo_sapiens.GRCh38.103.gtf.gz"
CELLRANGER = "genomics/homo_sapiens/10xgenomics/cellranger/"
VDJ_REFERENCE = CELLRANGER + "references/vdj/refdata-cellranger-vdj-GRCh38-alts-ensembl-5.0.0/"
VDJ_REFERENCE_FILES = ("reference.json", "fasta/regions.fa", "fasta/supp_regions.fa")
FASTQS = CELLRANGER + "5k_cmvpos_tcells/fastqs/vdj/subsampled_5k_human_antiCMV_T_TBNK_connect_VDJ_S1_L001_"
SAMPLE = "subsampled_5k_human_antiCMV_T_TBNK_connect_VDJ"
GENOME = "genomics/homo_sapiens/genome/"
PROBESET = CELLRANGER + "references/flex/Chromium_Human_Transcriptome_Probe_Set_v1.1.0_GRCh38-2024-A.chr22.csv"
TINY_BCL = "https://cf.10xgenomics.com/supp/cell-exp/"

# Upstream's test ext.args for mkgtf.
BIOTYPES = " ".join(
    f"--attribute=gene_biotype:{biotype}"
    for biotype in (
        "protein_coding",
        "lincRNA",
        "antisense",
        "IG_LV_gene",
        "IG_V_gene",
        "IG_V_pseudogene",
        "IG_D_gene",
        "IG_J_gene",
        "IG_J_pseudogene",
        "IG_C_gene",
        "IG_C_pseudogene",
        "TR_V_gene",
        "TR_V_pseudogene",
        "TR_D_gene",
        "TR_J_gene",
        "TR_J_pseudogene",
        "TR_C_gene",
    )
)


def file_md5(path: Path) -> str:
    """The md5 nf-test snapshots for ``path``: of its decompressed content when it ends in ``.gz``."""
    h = hashlib.md5()
    opener = gzip.open if path.name.endswith(".gz") else open
    with opener(path, "rb") as fh:
        while chunk := fh.read(1 << 20):
            h.update(chunk)
    return h.hexdigest()


def files_under(root: Path) -> dict[str, Path]:
    return {str(p.relative_to(root)): p for p in sorted(root.rglob("*")) if p.is_file()}


async def extracted(archive: File) -> dict[str, Path]:
    """``archive`` (a ``.tar.gz``) extracted locally: its files by path inside it."""
    root = Path(tempfile.mkdtemp())
    local = await archive.download()
    with tarfile.open(local) as tar:
        tar.extractall(root, filter="data")
    return files_under(root)


def named(files: dict[str, Path], match) -> list[str]:
    """``name:md5`` of the files whose path matches, sorted (nf-test's snapshot entries)."""
    return sorted(f"{Path(path).name}:{file_md5(p)}" for path, p in files.items() if match(path))


def assert_entries(actual: list[str], expected: list[str], label: str) -> None:
    if Counter(actual) != Counter(expected):
        raise AssertionError(f"{label}: expected {sorted(expected)}, got {sorted(actual)}")


def single_entry_mismatch(files: dict[str, Path], name: str, md5: str) -> list[str]:
    """[] when exactly one file is named ``name`` and has ``md5``, else a description of what was found."""
    actual = named(files, lambda p: Path(p).name == name)
    return [] if actual == [f"{name}:{md5}"] else [f"{name}: expected {md5}, got {actual}"]


def read_h5_matrix(path: Path) -> tuple[list[str], list[str], dict[str, dict[int, int]]]:
    """A Cell Ranger HDF5 matrix: its feature ids, barcodes, and each barcode's nonzero counts by feature row."""
    import h5py

    with h5py.File(path, "r") as h5:
        matrix = h5["matrix"]
        barcodes = [b.decode() for b in matrix["barcodes"][:]]
        features = [f.decode() for f in matrix["features"]["id"][:]]
        data, indices, indptr = matrix["data"][:], matrix["indices"][:], matrix["indptr"][:]
    columns = {
        barcode: {int(indices[i]): int(data[i]) for i in range(indptr[j], indptr[j + 1])}
        for j, barcode in enumerate(barcodes)
    }
    return features, barcodes, columns


def check_raw_holds_filtered(raw: Path, filtered: Path) -> None:
    raw_features, raw_barcodes, raw_columns = read_h5_matrix(raw)
    features, barcodes, columns = read_h5_matrix(filtered)
    assert raw_features == features, f"{raw.name}: features differ from {filtered.name}'s"
    assert barcodes, f"{filtered.name}: no cells"
    assert len(raw_barcodes) >= len(barcodes), f"{raw.name}: fewer barcodes than {filtered.name}"
    differing = [barcode for barcode in barcodes if raw_columns.get(barcode) != columns[barcode]]
    assert not differing, f"{raw.name}: {len(differing)} filtered cells differ, e.g. {differing[:3]}"


async def text_of(file: File) -> str:
    async with file.open("rb") as fh:
        return bytes(await fh.read()).decode()


# --- references -------------------------------------------------------------


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


async def filtered_reference(gtf_id: str, gtf: str, fasta: str, reference_name: str) -> File:
    """As upstream's test setup: mkgtf (``<gtf_id>.filtered``), then mkref with the filtered GTF."""
    filtered = await cellranger_mkgtf(await fixture(gtf), prefix=f"{gtf_id}.filtered", args=BIOTYPES)
    reference = await cellranger_mkref(await fixture(fasta), filtered, reference_name)
    return reference.reference


@env.task
async def test_cellranger_mkgtf() -> None:
    gtf = await cellranger_mkgtf(await fixture(GENOME + "genome.gtf"), prefix="homo_sapiens.filtered", args=BIOTYPES)
    assert gtf.name == "homo_sapiens.filtered.gtf", gtf.name
    await assert_md5(gtf, "6adb86eec59d00e503ee961261182066", label="homo_sapiens.filtered.gtf")


@env.task
async def test_cellranger_mkref() -> None:
    fasta, gtf = await fixture(GENOME + "genome.fasta"), await fixture(GENOME + "genome.gtf")
    reference = await cellranger_mkref(fasta, gtf, "homo_sapiens_reference")
    files = await extracted(reference.reference)
    ref = "homo_sapiens_reference/"
    assert_entries(
        named(files, lambda p: p.startswith(ref + "fasta/")),
        ["genome.fa:f315020d899597c1b57e5fe9f60f4c3e", "genome.fa.fai:3520cd30e1b100e55f578db9c855f685"],
        "fasta/",
    )
    assert_entries(
        named(files, lambda p: p.startswith(ref + "genes/")),
        ["genes.gtf.gz:50fc877b1c53b36b3b413aff88bda48c"],
        "genes/",
    )
    for path, md5 in (
        ("reference.json", "6f6c61f2b75b26d7b1f1b733f0dd1ebc"),
        ("star/SA", "bcf3e1a855783105150b46c905465333"),
        ("star/SAindex", "b93fb07d342e6c32a00ebc4311c0ad38"),
        ("star/Genome", "22102926fadf5890e905ca71b2da3f35"),
        ("star/geneInfo.tab", "8b608537307443ffaee4927d2b428805"),
        ("star/exonGeTrInfo.tab", "72d4dd88d25e3c5b0bb72b12b4ac99a8"),
        ("star/exonInfo.tab", "0d560290fab688b7268d88d5494bf9fe"),
        ("star/sjdbList.out.tab", "9e4f991abbbfeb3935a2bb21b9e258f1"),
    ):
        assert ref + path in files, f"mkref wrote no {path}"
        actual = file_md5(files[ref + path])
        assert actual == md5, f"{path}: md5 mismatch (expected {md5}, got {actual})"


async def check_count_bam(bam: File) -> None:
    """Content checks of the BAM, for upstream's ``getSamLinesMD5``.

    nft-bam hashes htsjdk's SAM text of every record, which ``samtools view``
    doesn't reproduce byte for byte (it gives b704262..., not upstream's
    9eff679f..., while every other output matches upstream). So: the header
    says coordinate-sorted, it holds the 10,000 records ``samtools view``
    counted, and mapped records carry Cell Ranger's raw UMI tag (``UR``).
    """
    sam = await samtools_view(bam=bam, args="-h")
    async with sam.open("rb") as fh:
        lines = bytes(await fh.read()).decode().splitlines()
    header = [line for line in lines if line.startswith("@")]
    records = [line.split("\t") for line in lines if not line.startswith("@")]
    assert any(line.startswith("@HD") and "SO:coordinate" in line for line in header), header[:3]
    assert len(records) == 10000, f"{len(records)} records, expected 10000"
    mapped = [r for r in records if not int(r[1]) & 4]
    assert mapped, "no mapped records"
    untagged = [r[0] for r in mapped if not any(t.startswith("UR:Z:") for t in r[11:])]
    assert not untagged, f"{len(untagged)} mapped records lack a UR tag, e.g. {untagged[:3]}"


@env.task
async def test_cellranger_count() -> None:
    reference = await filtered_reference(
        "homo_sapiens", GENOME + "genome.gtf", GENOME + "genome.fasta", "homo_sapiens_reference"
    )
    gex = CELLRANGER + "5k_cmvpos_tcells/fastqs/gex_1/subsampled_5k_human_antiCMV_T_TBNK_connect_GEX_1_S1_L001_"
    reads = [await fixture(gex + "R1_001.fastq.gz"), await fixture(gex + "R2_001.fastq.gz")]
    result = await cellranger_count(reads, reference, prefix="test_10x", args="--chemistry SC3Pv3 --create-bam true")
    files = await extracted(result.archive)
    unhashed = ("web_summary.html", "possorted_genome_bam.bam.bai", "possorted_genome_bam.bam", "barcodes.tsv.gz")
    assert_entries(
        named(files, lambda p: Path(p).name not in unhashed),
        [
            "filtered_feature_bc_matrix.h5:2f21402adbbbbe0793b6f70dd5502d7d",
            "features.tsv.gz:9f93621be0bede2b75596ad255607633",
            "matrix.mtx.gz:9ef2a9fab653234876c89236f710cb98",
            "metrics_summary.csv:14a544ef1204f4f9d128bc15febb94cf",
            "molecule_info.h5:25db52735645139ad09a83e5a0a17911",
            "raw_feature_bc_matrix.h5:dc30fe8b6d50a268b1d6a965998e016f",
            "features.tsv.gz:9f93621be0bede2b75596ad255607633",
            "matrix.mtx.gz:5846de1b0df533e2826ab175f4f610f4",
        ],
        "outs",
    )
    for name in ("web_summary.html", "barcodes.tsv.gz", "possorted_genome_bam.bam.bai"):
        assert any(Path(p).name == name for p in files), f"cellranger count wrote no {name}"
    assert result.bam is not None, "cellranger count wrote no possorted_genome_bam.bam"
    await check_count_bam(result.bam)


# --- multi ------------------------------------------------------------------


def check_config(text: str, present: tuple[str, ...], absent: tuple[str, ...]) -> None:
    for part in present:
        assert part in text, f"config lacks {part!r}:\n{text}"
    for part in absent:
        assert part not in text, f"config has {part!r}:\n{text}"


async def multi_files(archive: File) -> dict[str, Path]:
    """The outs archive's files, by path from the run folder's parent (``/<prefix>/outs/...``)."""
    return {"/" + path: p for path, p in (await extracted(archive)).items()}


async def library(path: str, lanes: tuple[str, ...] = ("L001",), **options) -> MultiLibrary:
    reads = [await fixture(f"{path}_{lane}_{read}_001.fastq.gz") for lane in lanes for read in ("R1", "R2")]
    return MultiLibrary(reads=reads, options=options)


async def chr22_reference() -> File:
    return await filtered_reference(
        "homo_sapiens_chr22", GENOME + "genome.gtf", GENOME + "genome.fasta", "homo_sapiens_chr22_reference"
    )


PBMC = CELLRANGER + "10k_pbmc/"
PBMC_CMO = CELLRANGER + "10k_pbmc_cmo/"
CMV_TCELLS = CELLRANGER + "5k_cmvpos_tcells/"
A549 = CELLRANGER + "sc3_v3_5k_a549_gex_crispr/"
NO_VDJ = ("[vdj]", "VDJ")
NO_FEATURE = ("[feature]", "Antibody Capture", "CRISPR Guide Capture", "fb_reference_copy.csv")
NO_CMO = ("cmo-set,", "[samples]", "Multiplexing Capture")


async def gex_5p(path: str) -> MultiLibrary:
    return await library(
        path, **{"expect-cells": "1000", "chemistry": "SC5P-R2", "create-bam": False, "no-secondary": True}
    )


@env.task
async def test_cellranger_multi_pbmc() -> None:
    prefix = "subsampled_sc5p_v2_hs_PBMC_10k"
    result = await cellranger_multi(
        prefix,
        gex=await gex_5p(PBMC + "fastqs/5gex/5gex/subsampled_sc5p_v2_hs_PBMC_10k_5gex_S1"),
        ab=await library(PBMC + "fastqs/5gex/5fb/subsampled_sc5p_v2_hs_PBMC_10k_5fb_S1"),
        vdj_t=await library(PBMC + "fastqs/tcell/subsampled_sc5p_v2_hs_PBMC_10k_t_S1"),
        vdj_b=await library(PBMC + "fastqs/bcell/subsampled_sc5p_v2_hs_PBMC_10k_b_S1"),
        gex_reference=await chr22_reference(),
        vdj_reference=await cellranger_vdj_test_reference(),
        fb_reference=await fixture(PBMC + "sc5p_v2_hs_PBMC_10k_multi_5gex_5fb_b_t_feature_ref.csv"),
    )
    files = await multi_files(result.archive)

    def name_is(name: str, nonempty: bool = False):
        return lambda p: Path(p).name == name and (not nonempty or files[p].stat().st_size > 0)

    # The subsampled T cell reads produce no contigs, so only non-empty contig files are hashed.
    assert_entries(
        named(files, name_is("all_contig.fasta", True)), ["all_contig.fasta:494adf38d81f749fc87409c982a1dd2d"], "fasta"
    )
    assert_entries(
        named(files, name_is("all_contig.fasta.fai")), ["all_contig.fasta.fai:95e967ffde85602135d22764664d9e3e"], "fai"
    )
    assert_entries(
        named(files, name_is("all_contig.fastq", True)), ["all_contig.fastq:e5b0cec096f4d57981916d2247a9f9b7"], "fastq"
    )
    assert_entries(
        named(files, lambda p: "consensus.bam.bai" in Path(p).name),
        ["consensus.bam.bai:125978e04496fe193a5c162603a7e6f6"] * 2,
        "consensus.bam.bai",
    )
    assert_entries(named(files, name_is("regions.fa")), ["regions.fa:1953406d913fb0fd753c6b56e847da83"], "regions.fa")
    assert_entries(
        named(files, lambda p: ".json" in Path(p).name),
        [
            "cell_barcodes.json:d751713988987e9331980363e24189ce",
            "cell_barcodes.json:d751713988987e9331980363e24189ce",
            "all_contig_annotations.json:cf33174fc4dbbadc5e74051bc13a5bd2",
            "reference.json:384f6efabad59f6da5c89b862aee71a8",
            "all_contig_annotations.json:d751713988987e9331980363e24189ce",
        ],
        "json",
    )
    vdj_t = sorted(Path(p).name for p in files if "/outs/vdj_t/" in p)
    expected = [
        "all_contig.bam",
        "all_contig.bam.bai",
        "all_contig.fasta",
        "all_contig.fastq",
        "all_contig_annotations.bed",
        "all_contig_annotations.csv",
        "all_contig_annotations.json",
    ]
    assert vdj_t == expected, f"outs/vdj_t/: expected {expected}, got {vdj_t} (all: {sorted(files)})"
    check_config(
        await text_of(result.config),
        present=(
            "[gene-expression]",
            "fastq_all/gex",
            "[vdj]",
            "fastq_all/vdj_t,,VDJ-T\n",
            "fastq_all/vdj_b,,VDJ-B",
            "[feature]",
            "Antibody Capture",
            "fb_reference_copy.csv",
        ),
        absent=("probe-set,", ",,VDJ\n", "CRISPR Guide Capture", *NO_CMO),
    )


CMO_NOT_REPRODUCED = (
    "cellranger multi stops on upstream's CellPlex test sample here (Deplex Error: no cell multiplexing tag "
    "sequences detected): on this cluster the gene-expression library calls 19 cells, and the marginal tag "
    "caller assigns no tag to any of them, though the CMO reads carry CMO301/CMO302. Upstream's run on its CI "
    "runners called tags; 4 CPUs and 16 or 6 GiB give the same failure here. The generated config is checked."
)


async def check_multi_cmo(probeset: bool) -> None:
    """Upstream's CellPlex cases: the config cellranger_multi would write, then NotTested (see CMO_NOT_REPRODUCED).

    Upstream snapshots ``assignment_confidence_table.csv`` (1bccc6e7...), ``tag_calls_summary.csv``
    (55d2c04c...), ``cells_per_tag.json`` (3f28cb95...), two ``metrics_summary.csv`` (292174f1..., 319f8b46...)
    and two ``sample_filtered_feature_bc_matrix.h5`` (929c9c72..., ca293555...).
    """
    prefix = "subsampled_SC3_v3_NextGem_DI_CellPlex_Human_PBMC_10K"
    gex_options = {"expect-cells": "1000", "chemistry": "SC3Pv3", "create-bam": False, "no-secondary": True}
    files = {"cmo_reference": "10k_pbmc_cmo_count_feature_reference.csv", "cmo_barcodes": "cmo_barcodes.csv"}
    if probeset:
        files["gex_frna_probeset"] = PROBESET.rsplit("/", 1)[-1]
    config = cellranger_multi_config(
        prefix,
        {"gex": MultiLibrary([], gex_options), "cmo": MultiLibrary([])},
        multi_reference_paths("homo_sapiens_chr22_reference.tar.gz", None, files),
    )
    check_config(
        config + "\n",
        present=("[gene-expression]", "fastq_all/gex", "cmo-set,", "[samples]", "Multiplexing Capture"),
        absent=("probe-set,", *NO_VDJ, *NO_FEATURE),
    )
    assert "$(cat references/cmo/barcodes/cmo_barcodes.csv)" in config, config
    raise NotTested(CMO_NOT_REPRODUCED)


@env.task
async def test_cellranger_multi_cmo() -> None:
    await check_multi_cmo(probeset=False)


@env.task
async def test_cellranger_multi_cmo_with_probeset() -> None:
    """Upstream's "pipeline scenario where different references might be available for different samples"."""
    await check_multi_cmo(probeset=True)


@env.task
async def test_cellranger_multi_vdj() -> None:
    prefix = "subsampled_5k_human_antiCMV_T_TBNK_connect"
    result = await cellranger_multi(
        prefix,
        gex=await gex_5p(CMV_TCELLS + "fastqs/gex_1/subsampled_5k_human_antiCMV_T_TBNK_connect_GEX_1_S1"),
        vdj=await library(CMV_TCELLS + "fastqs/vdj/subsampled_5k_human_antiCMV_T_TBNK_connect_VDJ_S1"),
        ab=await library(CMV_TCELLS + "fastqs/ab/subsampled_5k_human_antiCMV_T_TBNK_connect_AB_S2", lanes=("L004",)),
        gex_reference=await chr22_reference(),
        vdj_reference=await cellranger_vdj_test_reference(),
        gex_frna_probeset=await fixture(PROBESET),
        fb_reference=await fixture(
            CMV_TCELLS + "5k_human_antiCMV_T_TBNK_connect_Multiplex_count_feature_reference.csv"
        ),
    )
    files = await multi_files(result.archive)
    assert_entries(named(files, lambda p: Path(p).name == "all_contig.fasta.fai"), [], "fai")
    assert_entries(
        named(files, lambda p: "consensus.bam.bai" in Path(p).name),
        ["consensus.bam.bai:125978e04496fe193a5c162603a7e6f6"],
        "consensus.bam.bai",
    )
    assert_entries(
        named(files, lambda p: Path(p).name == "regions.fa"), ["regions.fa:1953406d913fb0fd753c6b56e847da83"], "regions"
    )
    assert_entries(
        named(files, lambda p: ".json" in Path(p).name),
        [
            "cell_barcodes.json:d751713988987e9331980363e24189ce",
            "reference.json:384f6efabad59f6da5c89b862aee71a8",
            "all_contig_annotations.json:d751713988987e9331980363e24189ce",
        ],
        "json",
    )
    check_config(
        await text_of(result.config),
        present=(
            "[gene-expression]",
            "fastq_all/gex",
            "[vdj]",
            ",,VDJ\n",
            "[feature]",
            "Antibody Capture",
            "fb_reference_copy.csv",
        ),
        absent=("probe-set,", "VDJ-", "CRISPR Guide Capture", *NO_CMO),
    )


@env.task
async def test_cellranger_multi_crispr() -> None:
    prefix = "subsampled_SC3_v3_NextGem_DI_CRISPR_A549_5K"
    reference = await filtered_reference(
        "homo_sapiens_chr1",
        A549 + "reference/genes_chr1_32292083_32686211.gtf",
        A549 + "reference/genome_chr1_32292083_32686211.fa",
        "homo_sapiens_chr1_32292083_32686211_reference",
    )
    lanes = ("L001", "L002")
    result = await cellranger_multi(
        prefix,
        gex=await library(
            A549 + "fastqs/gex/subsampled_SC3_v3_NextGem_DI_CRISPR_A549_5K_gex_S5",
            lanes,
            **{"expect-cells": "5000", "chemistry": "SC3Pv3", "create-bam": False},
        ),
        crispr=await library(A549 + "fastqs/crispr/subsampled_SC3_v3_NextGem_DI_CRISPR_A549_5K_crispr_S4", lanes),
        gex_reference=reference,
        fb_reference=await fixture(
            A549 + "SC3_v3_NextGem_DI_CRISPR_A549_5K_Multiplex_count_feature_reference_chr1_32292083_32686211.csv"
        ),
    )
    files = await multi_files(result.archive)
    mismatches = []
    for name, md5 in (
        ("metrics_summary.csv", "99b49df8fbb21832ff02451f970bc393"),
        ("raw_feature_bc_matrix.h5", "39fbcc1701143097d6899186ba660e67"),
        ("sample_filtered_feature_bc_matrix.h5", "005ca583fe9fe93004a40892bff63122"),
        ("protospacer_calls_per_cell.csv", "a94184bddcd8107421f62986b9e0984a"),
    ):
        mismatches += single_entry_mismatch(files, name, md5)
    assert not mismatches, "; ".join(mismatches)
    check_config(
        await text_of(result.config),
        present=("[gene-expression]", "fastq_all/gex", "[feature]", "CRISPR Guide Capture", "fb_reference_copy.csv"),
        absent=("probe-set,", *NO_VDJ, "Antibody Capture", *NO_CMO),
    )


@env.task
async def test_cellranger_multi_flex() -> None:
    prefix = "subsampled_10k_human_kidney_gemx_flex"
    flex = CELLRANGER + "singleplex_flex/Human_Kidney_GEM-X_Flex_S1_L001_"
    reads = [await fixture(flex + "R1_001.subsampled.fastq.gz"), await fixture(flex + "R2_001.subsampled.fastq.gz")]
    result = await cellranger_multi(
        prefix,
        gex=MultiLibrary(reads=reads, options={"create-bam": False}),
        gex_reference=await chr22_reference(),
        gex_frna_probeset=await fixture(PROBESET),
    )
    files = await multi_files(result.archive)
    mismatches = []
    for name, md5 in (
        ("metrics_summary.csv", "0baa073812fdd8df5c2e057a35b4746c"),
        ("filtered_feature_bc_matrix.h5", "ca5b896c260a17b2e88920c5f4564adb"),
        ("sample_filtered_feature_bc_matrix.h5", "bb7ac3fc77d4662b7315f229b8d9d4fd"),
    ):
        mismatches += single_entry_mismatch(files, name, md5)
    assert not mismatches, "; ".join(mismatches)
    # Upstream's raw matrices (2303d5a3..., 95e2f320...) differ here while the filtered ones match:
    # the raw ones must hold every filtered cell with the same counts.
    by_name = {Path(path).name: p for path, p in files.items()}
    check_raw_holds_filtered(by_name["raw_feature_bc_matrix.h5"], by_name["filtered_feature_bc_matrix.h5"])
    check_raw_holds_filtered(
        by_name["sample_raw_feature_bc_matrix.h5"], by_name["sample_filtered_feature_bc_matrix.h5"]
    )
    check_config(
        await text_of(result.config),
        present=("[gene-expression]", "fastq_all/gex", "probe-set,"),
        absent=(*NO_VDJ, *NO_FEATURE, *NO_CMO),
    )


# Upstream's stub cases, which check only the config: here, cellranger_multi_config itself.
PBMC_ID = "subsampled_sc5p_v2_hs_PBMC_10k"


def stub_config(libraries: dict[str, MultiLibrary], gex_reference: bool, **files: str) -> list[str]:
    references = multi_reference_paths(
        "homo_sapiens_chr22_reference.tar.gz" if gex_reference else None, "vdj_reference_stub", files
    )
    return cellranger_multi_config(PBMC_ID, libraries, references).split("\n")


@env.task
async def test_cellranger_multi_config_vdj_chains() -> None:
    lines = stub_config(
        {
            "gex": MultiLibrary(
                [], {"expect-cells": "1000", "chemistry": "SC5P-PE", "create-bam": False, "no-secondary": True}
            ),
            "vdj_t": MultiLibrary([], {"r1-length": 26}),
            "vdj_t_gd": MultiLibrary([]),
            "vdj_b": MultiLibrary([]),
        },
        gex_reference=True,
        gex_frna_probeset="probeset.csv",
        vdj_primer_index="gd_primers.txt",
    )
    expected = [
        "[gene-expression]",
        "reference,$PWD/references/gex/homo_sapiens_chr22_reference",
        "chemistry,SC5P-PE",
        "expect-cells,1000",
        "no-secondary,true",
        "create-bam,false",
        "[vdj]",
        "reference,$PWD/references/vdj/vdj_reference_stub",
        "inner-enrichment-primers,$PWD/references/vdj/primers/gd_primers.txt",
        "r1-length,26",
        "[libraries]",
        "fastq_id,fastqs,lanes,feature_types",
        f"{PBMC_ID},$PWD/fastq_all/gex,,Gene Expression",
        f"{PBMC_ID},$PWD/fastq_all/vdj_t,,VDJ-T",
        f"{PBMC_ID},$PWD/fastq_all/vdj_t_gd,,VDJ-T-GD",
        f"{PBMC_ID},$PWD/fastq_all/vdj_b,,VDJ-B",
    ]
    assert lines == expected, lines


@env.task
async def test_cellranger_multi_config_antibody_without_gex() -> None:
    lines = stub_config(
        {"ab": MultiLibrary([]), "vdj_t_gd": MultiLibrary([]), "vdj_b": MultiLibrary([])},
        gex_reference=True,
        vdj_primer_index="gd_primers.txt",
        fb_reference="feature_ref.csv",
    )
    expected = [
        "[gene-expression]",
        "reference,$PWD/references/gex/homo_sapiens_chr22_reference",
        "create-bam,true",
        "[feature]",
        "reference,$PWD/fb_reference_copy.csv",
        "[vdj]",
        "reference,$PWD/references/vdj/vdj_reference_stub",
        "inner-enrichment-primers,$PWD/references/vdj/primers/gd_primers.txt",
        "[libraries]",
        "fastq_id,fastqs,lanes,feature_types",
        f"{PBMC_ID},$PWD/fastq_all/vdj_t_gd,,VDJ-T-GD",
        f"{PBMC_ID},$PWD/fastq_all/vdj_b,,VDJ-B",
        f"{PBMC_ID},$PWD/fastq_all/ab,,Antibody Capture",
    ]
    assert lines == expected, lines


@env.task
async def test_cellranger_multi_config_gamma_delta_only() -> None:
    lines = stub_config({"vdj_t_gd": MultiLibrary([])}, gex_reference=True, vdj_primer_index="gd_primers.txt")
    expected = [
        "[vdj]",
        "reference,$PWD/references/vdj/vdj_reference_stub",
        "inner-enrichment-primers,$PWD/references/vdj/primers/gd_primers.txt",
        "[libraries]",
        "fastq_id,fastqs,lanes,feature_types",
        f"{PBMC_ID},$PWD/fastq_all/vdj_t_gd,,VDJ-T-GD",
    ]
    assert lines == expected, lines


@env.task
async def test_cellranger_multi_config_gamma_delta_without_primers() -> None:
    try:
        stub_config({"vdj_t_gd": MultiLibrary([])}, gex_reference=False)
    except ValueError as err:
        assert "VDJ-T-GD libraries require inner enrichment primers" in str(err), str(err)
    else:
        raise AssertionError("a VDJ-T-GD library without inner enrichment primers was accepted")


# --- mkfastq ----------------------------------------------------------------


async def check_mkfastq(samplesheet: str, index_metrics: str) -> None:
    run_dir = await fixture(TINY_BCL + "cellranger-tiny-bcl-1.2.0.tar.gz")
    result = await cellranger_mkfastq(
        await fixture(TINY_BCL + samplesheet), run_dir, prefix="test", args="--tiles 1101"
    )
    files = await extracted(result.fastq)
    assert_entries(
        named(files, lambda p: "_S1_" in Path(p).name and not Path(p).name.startswith("Undetermined")),
        [
            "test_sample_S1_L001_R1_001.fastq.gz:df3ddf85d83bd722219c646c7815dbdc",
            "test_sample_S1_L001_R2_001.fastq.gz:c8a387b517ff5b7e56b85ccbf882efd2",
            "test_sample_S1_L001_I1_001.fastq.gz:018b0dbd8daf76fb9f6dfa4729ea0453",
        ],
        "sample FASTQs",
    )
    assert_entries(
        named(files, lambda p: Path(p).name.startswith("Undetermined")),
        [
            "Undetermined_S0_L001_I1_001.fastq.gz:c93946d67bdb1d5b2ab955fd1dc9cb31",
            "Undetermined_S0_L001_R1_001.fastq.gz:faa79ce215ae72fdc3375cda066c6757",
            "Undetermined_S0_L001_R2_001.fastq.gz:d8955db4dfbf7aee89b6bdc3b7629f8f",
        ],
        "undetermined FASTQs",
    )
    interop = await result.interop.list_files()
    names = sorted(f.path.rstrip("/").rsplit("/", 1)[-1] for f in interop)
    assert names == ["IndexMetricsOut.bin"], names
    await assert_md5(interop[0], index_metrics, label="IndexMetricsOut.bin")


@env.task
async def test_cellranger_mkfastq_simple() -> None:
    await check_mkfastq("cellranger-tiny-bcl-simple-1.2.0.csv", "437d55bac453ddf603ef084ba2dab82e")


@env.task
async def test_cellranger_mkfastq_illumina() -> None:
    await check_mkfastq("cellranger-tiny-bcl-samplesheet-1.2.0.csv", "e2ac28695537336b53a509f48e93ed1f")


# --- V(D)J ------------------------------------------------------------------


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


tests = [
    test_cellranger_mkgtf,
    test_cellranger_mkref,
    test_cellranger_count,
    test_cellranger_multi_pbmc,
    test_cellranger_multi_cmo,
    test_cellranger_multi_cmo_with_probeset,
    test_cellranger_multi_vdj,
    test_cellranger_multi_crispr,
    test_cellranger_multi_flex,
    test_cellranger_multi_config_vdj_chains,
    test_cellranger_multi_config_antibody_without_gex,
    test_cellranger_multi_config_gamma_delta_only,
    test_cellranger_multi_config_gamma_delta_without_primers,
    test_cellranger_mkfastq_simple,
    test_cellranger_mkfastq_illumina,
    test_cellranger_mkvdjref,
    test_cellranger_vdj,
]
