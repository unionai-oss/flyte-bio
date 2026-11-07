"""Tests for flyte_bio.modules.multiqcsav.

Mirrors upstream's multiqcsav "NovaSeq6000" case (modules commit 9a1cc9f,
MultiQC 1.35, MultiQC-SAV 0.2.0): with ``-p``, as upstream's test config,
over a NovaSeq run folder plus a FastQC zip, with seqinspector's MultiQC
config. As upstream, the report's name and the names in ``multiqc_data/``
and ``multiqc_plots/{pdf,png,svg}/`` are checked; their content embeds the
run date.
"""

from flyte.io import Dir

from flyte_bio.modules.multiqcsav import multiqcsav
from tests.framework import env, fixture

RUN_FOLDER = "genomics/sarscov2/illumina/bcl/200624_A00834_0183_BHMTFYDRXX.tar.gz"
FASTQC_ZIP = "genomics/sarscov2/illumina/fastqc/test_fastqc.zip"
CONFIG = "https://raw.githubusercontent.com/nf-core/seqinspector/1.0.0/assets/multiqc_config.yml"

DATA = [
    "fastqc-status-check-heatmap.txt",
    "fastqc_overrepresented_sequences_plot.txt",
    "fastqc_per_base_n_content_plot.txt",
    "fastqc_per_base_sequence_quality_plot.txt",
    "fastqc_per_sequence_gc_content_plot_Counts.txt",
    "fastqc_per_sequence_gc_content_plot_Percentages.txt",
    "fastqc_per_sequence_quality_scores_plot.txt",
    "fastqc_sequence_counts_plot.txt",
    "fastqc_sequence_duplication_levels_plot.txt",
    "fastqc_sequence_length_distribution_plot.txt",
    "fastqc_top_overrepresented_sequences_table.txt",
    "llms-full.txt",
    "multiqc.log",
    "multiqc.parquet",
    "multiqc_citations.txt",
    "multiqc_data.json",
    "multiqc_fastqc.txt",
    "multiqc_general_stats.txt",
    "multiqc_sav.txt",
    "multiqc_sources.txt",
    "sav-intensity-vs-cycle-plot.txt",
    "sav-lane-metrics-summary-table.txt",
    "sav-pf-vs-occ-plot.txt",
    "sav-qscore-heatmap-plot.txt",
    "sav-qscore-histogram-plot.txt",
    "sav-read-metrics-summary-table.txt",
    "sav-summary-clusters-reads-lane-plot_Clusters.txt",
    "sav-summary-clusters-reads-lane-plot_Reads.txt",
]

PLOTS = [
    "fastqc-status-check-heatmap",
    "fastqc_overrepresented_sequences_plot",
    "fastqc_per_base_n_content_plot",
    "fastqc_per_base_sequence_quality_plot",
    "fastqc_per_sequence_gc_content_plot_Counts",
    "fastqc_per_sequence_gc_content_plot_Percentages",
    "fastqc_per_sequence_quality_scores_plot",
    "fastqc_sequence_counts_plot-cnt",
    "fastqc_sequence_counts_plot-pct",
    "fastqc_sequence_duplication_levels_plot",
    "fastqc_sequence_length_distribution_plot",
    "fastqc_top_overrepresented_sequences_table",
    "sav-intensity-vs-cycle-plot",
    "sav-lane-metrics-summary-table",
    "sav-pf-vs-occ-plot",
    "sav-qscore-heatmap-plot",
    "sav-qscore-histogram-plot",
    "sav-read-metrics-summary-table",
    "sav-summary-clusters-reads-lane-plot_Clusters-cnt",
    "sav-summary-clusters-reads-lane-plot_Clusters-pct",
    "sav-summary-clusters-reads-lane-plot_Reads-cnt",
    "sav-summary-clusters-reads-lane-plot_Reads-pct",
]


async def children(results: Dir, folder: str) -> set[str]:
    """Names of the entries directly under ``folder/`` in ``results``."""
    marker = f"/{folder}/"
    names = set()
    async for f in results.walk():
        if marker in f.path:
            names.add(f.path.split(marker, 1)[1].split("/", 1)[0])
    return names


@env.task
async def test_multiqcsav_novaseq6000() -> None:
    # upstream case: multiqcsav "NovaSeq6000"
    run_folder, zipped, config = await fixture(RUN_FOLDER), await fixture(FASTQC_ZIP), await fixture(CONFIG)
    result = await multiqcsav(run_folder, extra=[zipped], configs=[config], args="-p")
    assert result.report.name == "multiqc_report.html", result.report.name
    data = await children(result.results, "multiqc_data")
    assert data == set(DATA), f"multiqc_data: missing {set(DATA) - data}, extra {data - set(DATA)}"
    for kind in ("pdf", "png", "svg"):
        plots = await children(result.results, f"multiqc_plots/{kind}")
        expected = {f"{p}.{kind}" for p in PLOTS}
        assert plots == expected, f"multiqc_plots/{kind}: missing {expected - plots}, extra {plots - expected}"


tests = [test_multiqcsav_novaseq6000]
