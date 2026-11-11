"""Tests for flyte_bio.modules.plasmidid.

The case is upstream's plasmidid test (modules commit 9a1cc9f, plasmidID
1.6.5) with its config's ``-k 0.8``. As upstream: the final HTML and table
and the matched reference's FASTA by md5; the files of ``images/``,
``logs/`` (without mash's ``mashclust_*.log``), ``data/``, ``database/``
and ``kmer/`` by name.
"""

from flyte.io import Dir

from flyte_bio.modules.plasmidid import plasmidid
from tests.framework import assert_md5, env, fixture

EXPECTED_NAMES = {
    "data": [
        "database_reconstruct_test.length",
        "pID_highlights.conf",
        "pID_text_annotation.coordinates",
        "test.bedgraph_term",
        "test.coverage_adapted_clustered_ac",
        "test.err",
        "test.fna",
        "test.gbk",
        "test.gff",
        "test.gff.bed",
        "test.gff.coordinates",
        "test.gff.forward.bed",
        "test.gff.forward.coordinates",
        "test.gff.reverse.bed",
        "test.gff.reverse.coordinates",
        "test.karyotype_individual.txt",
        "test.karyotype_summary.txt",
        "test.log",
        "test.plasmids.bed",
        "test.plasmids.blast",
        "test.plasmids.blast.links",
        "test.plasmids.complete",
        "test.plasmids.links",
    ],
    "database": [
        "test.err",
        "test.fna",
        "test.gbk",
        "test.gff",
        "test.gff.bed",
        "test.gff.forward.bed",
        "test.gff.renamed",
        "test.gff.reverse.bed",
        "test.log",
    ],
    "fasta_files": ["MT192765.1_term.fasta"],
    "images": [
        "test_MT192765.1.png",
        "test_MT192765.1_individual.circos.conf",
        "test_individual.circos.conf",
        "test_summary.circos.conf",
        "test_summary.png",
    ],
    "kmer": [
        "database.filtered_0.8",
        "database.filtered_0.8_term.0.5.clusters.tab",
        "database.filtered_0.8_term.0.5.representative.fasta",
        "database.filtered_0.8_term.0.5.representative.fasta.blast.tmp.ndb",
        "database.filtered_0.8_term.0.5.representative.fasta.blast.tmp.nhr",
        "database.filtered_0.8_term.0.5.representative.fasta.blast.tmp.nin",
        "database.filtered_0.8_term.0.5.representative.fasta.blast.tmp.not",
        "database.filtered_0.8_term.0.5.representative.fasta.blast.tmp.nsq",
        "database.filtered_0.8_term.0.5.representative.fasta.blast.tmp.ntf",
        "database.filtered_0.8_term.0.5.representative.fasta.blast.tmp.nto",
        "database.filtered_0.8_term.fasta",
        "database.filtered_0.8_term.mash.distances.tab",
        "database.msh",
        "database.screen.tab",
    ],
    "logs": ["draw_circos_images.log", "plasmidID.log"],
}


async def plasmidid_subdir_files(report: Dir) -> dict[str, dict]:
    """The report's files, by subdirectory (top-level files under "")."""
    root = report.path.rstrip("/") + "/"
    found: dict[str, dict] = {}
    async for file in report.walk(recursive=True):
        relative = file.path[len(root) :] if file.path.startswith(root) else file.name
        subdir, _, name = relative.rpartition("/")
        found.setdefault(subdir, {})[name] = file
    return found


@env.task
async def test_plasmidid_contigs() -> None:
    # upstream case: plasmidid "test-plasmidid"
    result = await plasmidid(
        await fixture("genomics/sarscov2/illumina/fasta/contigs.fasta"),
        await fixture("genomics/sarscov2/genome/genome.fasta"),
        prefix="test",
        args="-k 0.8",
    )
    await assert_md5(result.html, "415774b2cff3a50e4304c875cf72529d", label=result.html.name)
    await assert_md5(result.tab, "e7e6be8e76cac944979bb3dfc21da8bd", label=result.tab.name)
    files = await plasmidid_subdir_files(result.report)
    for subdir, expected in EXPECTED_NAMES.items():
        names = sorted(n for n in files.get(subdir, {}) if not (n.startswith("mashclust_") and n.endswith(".log")))
        assert names == sorted(expected), f"{subdir}/: {names}"
    await assert_md5(
        files["fasta_files"]["MT192765.1_term.fasta"], "8a8537dd3b21e6905f9367b51c3b3074", label="MT192765.1_term.fasta"
    )


tests = [test_plasmidid_contigs]
