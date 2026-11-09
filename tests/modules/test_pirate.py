"""Tests for flyte_bio.modules.pirate.

The case is upstream's pirate test (modules commit 9a1cc9f, PIRATE 1.0.5).
Upstream treats every output as unstable and snapshots only the names of the
entries in ``<prefix>_results`` (files and directories), and that there's no
core alignment; this test asserts the same names.
"""

from flyte_bio.modules.pirate import pirate
from tests.framework import env, fixture

GFF = "genomics/prokaryotes/candidatus_portiera_aleyrodidarum/genome/gff/"
ENTRIES = [
    "PIRATE.gene_families.ordered.tsv",
    "PIRATE.gene_families.tsv",
    "PIRATE.genomes_per_allele.tsv",
    "PIRATE.log",
    "PIRATE.pangenome_summary.txt",
    "PIRATE.unique_alleles.tsv",
    "binary_presence_absence.fasta",
    "binary_presence_absence.nwk",
    "cluster_alleles.tab",
    "co-ords",
    "genome2loci.tab",
    "genome_list.txt",
    "link_clusters.log",
    "loci_list.tab",
    "loci_paralog_categories.tab",
    "modified_gffs",
    "pan_sequences.fasta",
    "pangenome.connected_blocks.tsv",
    "pangenome.edges",
    "pangenome.gfa",
    "pangenome.order.tsv",
    "pangenome.reversed.tsv",
    "pangenome.syntenic_blocks.tsv",
    "pangenome.temp",
    "pangenome_iterations",
    "pangenome_log.txt",
    "paralog_clusters.tab",
    "representative_sequences.faa",
    "representative_sequences.ffn",
    "split_groups.log",
]


@env.task
async def test_pirate_three_genomes() -> None:
    # upstream case: pirate "test-pirate"
    gff = [await fixture(f"{GFF}test{n}.gff") for n in (1, 2, 3)]
    result = await pirate(gff, prefix="test")
    base = result.results.path.rstrip("/")
    entries = sorted({f.path[len(base) :].lstrip("/").split("/")[0] async for f in result.results.walk()})
    assert entries == ENTRIES, f"test_results holds {entries}, expected {ENTRIES}"
    assert result.aln is None, "pirate wrote a core alignment without -a"


tests = [test_pirate_three_genomes]
