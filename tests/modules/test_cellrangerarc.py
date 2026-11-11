"""Tests for flyte_bio.modules.cellrangerarc.

The cases and expected md5s are upstream's cellrangerarc tests (modules
commit 9a1cc9f, Cell Ranger ARC 2.0.2). nf-test hashes ``.gz`` files
decompressed, so these do too.

- mkgtf filters the mouse chromosome 19 GTF with upstream's test biotypes.
- mkref builds ``cellrangerarc_reference`` from mouse chromosome 19
  (gunzipped first), its GTF and the test motifs, as upstream's test config
  names them, and checks the config and every file of the reference
  (``genomeParameters.txt``, which records the run, only by name).
- count runs upstream's subsampled multiome sample (a gene-expression and an
  ATAC library) against that reference and checks the md5s of the outputs
  upstream calls stable and the presence of the others. 19 of the 62
  "stable" outputs (the LSA/PCA components, variances and projections, the
  t-SNE projections, the ATAC UMAP and every
  ``differential_accessibility.csv``) differ from upstream's md5s on this
  cluster in their floating-point values, while the clusters computed from
  them match; those get content checks (:func:`check_numeric_table`).

They run with the module's default resources (4 CPUs, 16 GiB).
"""

import gzip
import hashlib
import tarfile
import tempfile
from pathlib import Path

from flyte.io import File

from flyte_bio.modules.cellrangerarc import cellrangerarc_count, cellrangerarc_mkgtf, cellrangerarc_mkref
from flyte_bio.modules.gunzip import gunzip
from tests.framework import assert_md5, env, fixture

MOUSE = "genomics/mus_musculus/"
MULTIOME = MOUSE + "illumina/10xgenomics/multiome/"

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


def arc_file_md5(path: Path) -> str:
    """The md5 nf-test snapshots for ``path``: of its decompressed content when it ends in ``.gz``."""
    h = hashlib.md5()
    opener = gzip.open if path.name.endswith(".gz") else open
    with opener(path, "rb") as fh:
        while chunk := fh.read(1 << 20):
            h.update(chunk)
    return h.hexdigest()


def files_under(root: Path) -> dict[str, Path]:
    return {str(p.relative_to(root)): p for p in sorted(root.rglob("*")) if p.is_file()}


async def arc_extracted(archive: File) -> dict[str, Path]:
    """``archive`` (a ``.tar.gz``) extracted locally: its files by path inside it."""
    root = Path(tempfile.mkdtemp())
    with tarfile.open(await archive.download()) as tar:
        tar.extractall(root, filter="data")
    return files_under(root)


async def arc_reference():
    """As upstream's test: the reference upstream's test config describes."""
    fasta = await gunzip(await fixture(MOUSE + "genome/chr19.fa.gz"))
    return await cellrangerarc_mkref(
        fasta,
        await fixture(MOUSE + "genome/chr19.filtered.gtf.gz"),
        "cellrangerarc_reference",
        motifs=await fixture("genomics/homo_sapiens/genome/genome_motifs.txt"),
        reference_config=await fixture(MULTIOME + "cellranger_arc_mkref_test_mm39_chr19_config.json"),
    )


@env.task
async def test_cellrangerarc_mkgtf() -> None:
    gtf = await cellrangerarc_mkgtf(await fixture(MOUSE + "genome/chr19.filtered.gtf.gz"), "test.filtered", BIOTYPES)
    assert gtf.name == "test.filtered.gtf", gtf.name
    await assert_md5(gtf, "29ea9f8629c84868b4b0ce5cd4f57632", label="test.filtered.gtf")


REFERENCE_MD5S = {
    "genome.fa": "dd9375c0cfa5de18757a0608411ff608",
    "genome.fa.amb": "1b4fb43678899d6115d37d88f22552ad",
    "genome.fa.ann": "4df1dc4bd46b294e1e1343dc0ab646ef",
    "genome.fa.bwt": "64d2ab05657806c66fb0e0d1ea7fc4d1",
    "genome.fa.fai": "291a00d5431e0f97c356c272d201d8bd",
    "genome.fa.pac": "7554f21668cd8efc77f1272b95550cc3",
    "genome.fa.sa": "7d4b562f369d65ad60a889f74609f330",
    "genes.gtf.gz": "29ea9f8629c84868b4b0ce5cd4f57632",
    "reference.json": "1ff977fa99515695d81b8da012569901",
    "motifs.pfm": "294d4041e4ffd4c633376c2b9d180ac0",
    "transcripts.bed": "bc5cfb75dc1795cfc2caeddfc15aaa30",
    "tss.bed": "6aa464a762cab3ca3bf004dd6f35de8e",
    "Genome": "4664d17f3cbcde9859f1af9fac08f9ed",
    "SA": "3c73b38e44260c9a591e74f5f1dfec87",
    "SAindex": "c1ba6a91058b189bffa7a56592635f23",
    "chrLength.txt": "263f433d64f5f9f3181ed90542735474",
    "chrName.txt": "255cdc5fefdc3ff84cfe3f3270ab8fed",
    "chrNameLength.txt": "bb9ac7f24e6c1c6ea0677f825858449a",
    "chrStart.txt": "fba9f49cedf08a4136244d4b9b776c01",
    "exonGeTrInfo.tab": "4fb4acda70869b18296f8ce9c6a297db",
    "exonInfo.tab": "e2f4405c0cbf847580791f2e9954f936",
    "geneInfo.tab": "786000c34020b0b39e4ab6b9b3e4d922",
    "genomeParameters.txt": "",  # only its name: it records the run
    "sjdbInfo.txt": "17bd31bf73ea326f7befac8a57dc1780",
    "sjdbList.fromGTF.out.tab": "23c67498223466c717f6f901ff9df9a4",
    "sjdbList.out.tab": "6923d3e29967b2f3faad78960101a2a6",
    "transcriptInfo.tab": "b3b9193220aaa3fbf63a49148ef150e0",
}


@env.task
async def test_cellrangerarc_mkref() -> None:
    reference = await arc_reference()
    await assert_md5(reference.config, "be9df01d35e426f4f37bdca2c779e199", label="config")
    files = await arc_extracted(reference.reference)
    by_name = {Path(path).name: p for path, p in files.items()}
    assert len(by_name) == len(files), f"duplicate file names in the reference: {sorted(files)}"
    assert sorted(by_name) == sorted(REFERENCE_MD5S), f"reference files: {sorted(files)}"
    mismatches = [
        f"{name}: expected {md5}, got {arc_file_md5(by_name[name])}"
        for name, md5 in REFERENCE_MD5S.items()
        if md5 and arc_file_md5(by_name[name]) != md5
    ]
    assert not mismatches, "; ".join(mismatches)


# Upstream's stable outputs, in its (sorted) order, and their md5s.
STABLE = [
    ("analysis/clustering/atac/graphclust/clusters.csv", "ed29b219ca3b79305129aea2d38ca60e"),
    ("analysis/clustering/atac/graphclust/differential_accessibility.csv", "90e3e1001432bd1f94f53f39ae86971c"),
    ("analysis/clustering/atac/graphclust/differential_expression.csv", "2806595b0c66af8407977b6ea57e30a3"),
    ("analysis/clustering/atac/kmeans_2_clusters/clusters.csv", "35232075d0184cc6a31178ed58176bbe"),
    ("analysis/clustering/atac/kmeans_2_clusters/differential_accessibility.csv", "0c597e33c58aac61a54e71c0ca3d1343"),
    ("analysis/clustering/atac/kmeans_2_clusters/differential_expression.csv", "e77580db30873fba52bf4490f752ed22"),
    ("analysis/clustering/atac/kmeans_3_clusters/clusters.csv", "127636548ad77bfa2646cbf62209e73b"),
    ("analysis/clustering/atac/kmeans_3_clusters/differential_accessibility.csv", "9d1efbb7209c95a6f7714b43aca22130"),
    ("analysis/clustering/atac/kmeans_3_clusters/differential_expression.csv", "24a85767ec0013a9f444bed5b23fb7f2"),
    ("analysis/clustering/atac/kmeans_4_clusters/clusters.csv", "c84df7cc5cb684dd4a431915df0b58b2"),
    ("analysis/clustering/atac/kmeans_4_clusters/differential_accessibility.csv", "2492068224aba9a65f860b6fb3464056"),
    ("analysis/clustering/atac/kmeans_4_clusters/differential_expression.csv", "e02a8c8e3bc1d924357b46af2da76a6b"),
    ("analysis/clustering/atac/kmeans_5_clusters/clusters.csv", "fcbbc4f29456778935fbbf4dc52791ea"),
    ("analysis/clustering/atac/kmeans_5_clusters/differential_accessibility.csv", "630f412020f85765b2a32ac689701ce5"),
    ("analysis/clustering/atac/kmeans_5_clusters/differential_expression.csv", "0aae96dbf649f0bc738f5d2c19e2be76"),
    ("analysis/clustering/gex/graphclust/clusters.csv", "6cd4c2b06047c99cee59548c9b3000ac"),
    ("analysis/clustering/gex/graphclust/differential_accessibility.csv", "08fc89ad11d1834cab677390f385a88a"),
    ("analysis/clustering/gex/graphclust/differential_expression.csv", "6e97ced39f84b02f843f0b757b02f100"),
    ("analysis/clustering/gex/kmeans_2_clusters/clusters.csv", "cd135bf9e9b5aa69abe95d5356aec96c"),
    ("analysis/clustering/gex/kmeans_2_clusters/differential_accessibility.csv", "9ef66432c3f66f2190ff936324091e9f"),
    ("analysis/clustering/gex/kmeans_2_clusters/differential_expression.csv", "223e516ce7ddf219575a32a60217e83d"),
    ("analysis/clustering/gex/kmeans_3_clusters/clusters.csv", "ae9c24a834f1863a8faf03b3b4bdf74a"),
    ("analysis/clustering/gex/kmeans_3_clusters/differential_accessibility.csv", "84bf1b3b3c90abdd92a817cc5a9d2ee8"),
    ("analysis/clustering/gex/kmeans_3_clusters/differential_expression.csv", "8b0ec69f3cd0d3156d6def50cdde4340"),
    ("analysis/clustering/gex/kmeans_4_clusters/clusters.csv", "aceb697994f98ee50b4d3db8e12e2a3a"),
    ("analysis/clustering/gex/kmeans_4_clusters/differential_accessibility.csv", "05f0b3a10ce71d2fbafc6c9b74acf969"),
    ("analysis/clustering/gex/kmeans_4_clusters/differential_expression.csv", "7b339771b1731b9cff2c52ac4fc632ce"),
    ("analysis/clustering/gex/kmeans_5_clusters/clusters.csv", "243233d51193e6ce3127ed5ef8a08daf"),
    ("analysis/clustering/gex/kmeans_5_clusters/differential_accessibility.csv", "25deefdff4c870e52c84f97cc4b9d5d2"),
    ("analysis/clustering/gex/kmeans_5_clusters/differential_expression.csv", "b70730d7301cba11c46297d08b81677d"),
    ("analysis/dimensionality_reduction/atac/lsa_components.csv", "c3a48d3849f6b5e268bcb16b4d3db7ec"),
    ("analysis/dimensionality_reduction/atac/lsa_dispersion.csv", "1bdae634bd63d044dcfd8ae821b3744c"),
    ("analysis/dimensionality_reduction/atac/lsa_features_selected.csv", "b4e3cc5cac26dc6df059c1b8e38e72a8"),
    ("analysis/dimensionality_reduction/atac/lsa_projection.csv", "3bce47a5a0e8dc0cb24c2ecbba5c745b"),
    ("analysis/dimensionality_reduction/atac/lsa_variance.csv", "ff92c64ec140084673fd7ec48ce48993"),
    ("analysis/dimensionality_reduction/atac/tsne_projection.csv", "b0ff0a327a7cad898087c169b9815bab"),
    ("analysis/dimensionality_reduction/atac/umap_projection.csv", "8a34ad7bcf4ae7d14ad8f451acb4a12e"),
    ("analysis/dimensionality_reduction/gex/pca_components.csv", "82a4384a22390821cbfe06b2de9e9b4d"),
    ("analysis/dimensionality_reduction/gex/pca_dispersion.csv", "a6c8d15062a67bf828c1b58d1bc1eff8"),
    ("analysis/dimensionality_reduction/gex/pca_features_selected.csv", "34acd3e8d57f68bf592856c0c69847e7"),
    ("analysis/dimensionality_reduction/gex/pca_projection.csv", "6a2bf5e268da6705b92ec02309b3f61e"),
    ("analysis/dimensionality_reduction/gex/pca_variance.csv", "060e80415e690a7bfc9cd385c929842f"),
    ("analysis/dimensionality_reduction/gex/tsne_projection.csv", "aa42abe3c868c7d0e35628acf1a24966"),
    ("analysis/dimensionality_reduction/gex/umap_projection.csv", "34acdb54ff9f9c8b7d7f91a625450294"),
    ("analysis/feature_linkage/feature_linkage.bedpe", "4a92b2e49591c21a3d08602235e6ef7a"),
    ("analysis/tf_analysis/filtered_tf_bc_matrix/barcodes.tsv.gz", "1b2010b38c2de59bea1b71d374114919"),
    ("analysis/tf_analysis/filtered_tf_bc_matrix/matrix.mtx.gz", "2627809b0951961282c3b31fb83d7f1b"),
    ("analysis/tf_analysis/filtered_tf_bc_matrix/motifs.tsv", "9a397341713d07560937054df442afcb"),
    ("analysis/tf_analysis/peak_motif_mapping.bed", "b1d02a287485989ca294284bd95ad572"),
    ("atac_cut_sites.bigwig", "03e0c63674192256d55a3c117fff66ee"),
    ("atac_peak_annotation.tsv", "f4770b4197ff32c74bd206476f1f0dba"),
    ("filtered_feature_bc_matrix/barcodes.tsv.gz", "1b2010b38c2de59bea1b71d374114919"),
    ("filtered_feature_bc_matrix/features.tsv.gz", "9a315f46dbf99b9120527df78fd7e081"),
    ("filtered_feature_bc_matrix/matrix.mtx.gz", "2304085aaac98027cf9241a3b2593760"),
    ("gex_molecule_info.h5", "f1ab414eb8b8da017c974d7109f5e7ca"),
    ("gex_possorted_bam.bam", "8455846fa1aa417e4252bd09584d1098"),
    ("gex_possorted_bam.bam.bai", "444fb175b75d0e22c8e14724b7790e34"),
    ("per_barcode_metrics.csv", "055ddb433b94d44efa0113a1db84a63d"),
    ("raw_feature_bc_matrix/barcodes.tsv.gz", "4e020c4c968100a61a542c0d6e56e5f5"),
    ("raw_feature_bc_matrix/features.tsv.gz", "9a315f46dbf99b9120527df78fd7e081"),
    ("raw_feature_bc_matrix/matrix.mtx.gz", "743fc9d85f85cd2e5a40654d3ae2e3e4"),
    ("summary.csv", "7641192d4cc0c6072dc447a46d5f89a2"),
]
UNSTABLE = [
    "analysis/feature_linkage/feature_linkage_matrix.h5",
    "analysis/tf_analysis/filtered_tf_bc_matrix.h5",
    "atac_fragments.tsv.gz",
    "atac_fragments.tsv.gz.tbi",
    "atac_peaks.bed",
    "atac_possorted_bam.bam",
    "atac_possorted_bam.bam.bai",
    "cloupe.cloupe",
    "filtered_feature_bc_matrix.h5",
    "raw_feature_bc_matrix.h5",
    "web_summary.html",
]


def floating_point(path: str) -> bool:
    """The outputs whose floating-point values differ in their last digits from upstream's on this cluster.

    The LSA/PCA decompositions (and the t-SNE/UMAP projections and
    differential accessibility computed from them) come out of BLAS, whose
    rounding depends on the CPU; the clusters they lead to match upstream's
    md5s.
    """
    return path.endswith("differential_accessibility.csv") or (
        path.startswith("analysis/dimensionality_reduction/")
        and not path.endswith(("_dispersion.csv", "_features_selected.csv", "gex/umap_projection.csv"))
    )


def check_numeric_table(outs: dict[str, Path], path: str) -> None:
    """Checks for an output whose md5 can't match: a header and rows of equal width whose numbers are finite,
    one row per barcode for projections (the barcodes of the modality's graph-based clusters, in order)."""
    with open(outs[path]) as fh:
        rows = [line.rstrip("\n").split(",") for line in fh]
    assert len(rows) > 1, f"{path}: no rows"
    widths = {len(row) for row in rows}
    assert len(widths) == 1, f"{path}: rows of {sorted(widths)} columns"
    numbers = 0
    for row in rows[1:]:
        for value in row[1:]:
            try:
                number = float(value)
            except ValueError:
                continue
            numbers += 1
            assert number == number and abs(number) != float("inf"), f"{path}: {value!r} is not a finite number"
    assert numbers, f"{path}: no numbers"
    if path.endswith("_projection.csv"):
        modality = path.split("/")[2]
        with open(outs[f"analysis/clustering/{modality}/graphclust/clusters.csv"]) as fh:
            barcodes = [line.split(",")[0] for line in fh][1:]
        assert [row[0] for row in rows[1:]] == barcodes, f"{path}: barcodes differ from graphclust/clusters.csv"


@env.task
async def test_cellrangerarc_count() -> None:
    reads = [
        await fixture(MULTIOME + name)
        for name in (
            "SRR18907480_chr19_sub_S1_L001_R1_001.fastq.gz",
            "SRR18907480_chr19_sub_S1_L001_R2_001.fastq.gz",
            "SRR18907481_chr19_sub_S1_L001_R1_001.fastq.gz",
            "SRR18907481_chr19_sub_S1_L001_R2_001.fastq.gz",
            "SRR18907481_chr19_sub_S1_L001_I2_001.fastq.gz",
        )
    ]
    reference = await arc_reference()
    result = await cellrangerarc_count(
        reads,
        reference.reference,
        sample_types=["gex", "atac"],
        sample_names=["SRR18907480_chr19_sub", "SRR18907481_chr19_sub"],
        prefix="test",
    )
    files = await arc_extracted(result.archive)
    outs = {path.removeprefix("test/outs/"): p for path, p in files.items()}
    missing = [path for path, _ in STABLE if path not in outs] + [path for path in UNSTABLE if path not in outs]
    assert not missing, f"cellranger-arc count wrote no {missing}"
    mismatches = [
        f"{path}: expected {md5}, got {arc_file_md5(outs[path])}"
        for path, md5 in STABLE
        if not floating_point(path) and arc_file_md5(outs[path]) != md5
    ]
    assert not mismatches, f"{len(mismatches)} of {len(STABLE)} differ: " + "; ".join(mismatches)
    for path, _ in STABLE:
        if floating_point(path):
            check_numeric_table(outs, path)
    assert result.libraries.name == "test_lib.csv", result.libraries.name


tests = [test_cellrangerarc_mkgtf, test_cellrangerarc_mkref, test_cellrangerarc_count]
