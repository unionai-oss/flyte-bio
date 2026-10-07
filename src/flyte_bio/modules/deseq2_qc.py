"""deseq2_qc — sample-level QC of a gene count matrix with DESeq2.

Exposes :func:`deseq2_qc`, which runs upstream's ``deseq2_qc.r`` on a gene
by-sample count table: a variance-stabilised (or rlog) transform, PCA and
sample-distance clustering, with plots, the R objects and size factors. As
upstream, the PCA and distance tables are then prefixed with MultiQC
headers relabelled for ``label`` (the aligner name, e.g. ``star_salmon``),
giving ``<label>.pca.vals_mqc.tsv`` and ``<label>.sample.dists_mqc.tsv``.

The R script is vendored verbatim (it's an optparse CLI) and staged, with
upstream's two MultiQC headers, into the stock R/DESeq2 container.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.scripts import path

DESEQ2_QC_IMAGE = "community.wave.seqera.io/library/r-base_r-optparse_r-ggplot2_r-rcolorbrewer_pruned:9e75394d0bc21987"

DESEQ2_QC_CPUS = 2
DEFAULT_RESOURCES = flyte.Resources(cpu=DESEQ2_QC_CPUS, memory="8Gi")

# The script's own options are typed inputs (not a free `args` string) so an
# empty --sample_suffix stays an empty argument, as upstream's `''` does.
deseq2_qc_cmd = shell.create(
    name="deseq2_qc",
    image=DESEQ2_QC_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "script": File,
        "counts": File,
        "pca_header": File,
        "clustering_header": File,
        "label": str,
        "prefix": str,
        "id_col": int,
        "count_col": int,
        "sample_suffix": str,
        "vst": bool,
    },
    outputs={"results": Dir},
    script=rf"""
        SCRIPT=({{inputs.script}}); COUNTS=({{inputs.counts}})
        PCA_HEADER=({{inputs.pca_header}}); CLUSTERING_HEADER=({{inputs.clustering_header}})
        LABEL={{inputs.label}}
        label_lower=${{LABEL,,}}; label_upper=${{LABEL^^}}
        VST=()
        [ {{inputs.vst}} = true ] && VST=(--vst TRUE)
        cd {{outputs.results}}
        Rscript "${{SCRIPT[0]}}" \
            --count_file "${{COUNTS[0]}}" \
            --outdir ./ \
            --cores {DESEQ2_QC_CPUS} \
            --outprefix {{inputs.prefix}} \
            --id_col {{inputs.id_col}} \
            --sample_suffix {{inputs.sample_suffix}} \
            --count_col {{inputs.count_col}} \
            "${{VST[@]}}"
        if [ -f "R_sessionInfo.log" ]; then
            # Handle PCA files
            sed "s/deseq2_pca/${{label_lower}}_deseq2_pca/g" <"${{PCA_HEADER[0]}}" > pca_header.tmp
            sed -i -e "s/DESeq2 PCA/${{label_upper}} DESeq2 PCA/g" pca_header.tmp
            cat pca_header.tmp *.pca.vals.txt > ${{label_lower}}.pca.vals_mqc.tsv
            rm pca_header.tmp
            # Handle clustering files
            sed "s/deseq2_clustering/${{label_lower}}_deseq2_clustering/g" <"${{CLUSTERING_HEADER[0]}}" > clustering_header.tmp
            sed -i -e "s/DESeq2 sample/${{label_upper}} DESeq2 sample/g" clustering_header.tmp
            cat clustering_header.tmp *.sample.dists.txt > ${{label_lower}}.sample.dists_mqc.tsv
            rm clustering_header.tmp
        fi
    """,
)


env = flyte.TaskEnvironment.from_task(
    "deseq2_qc",
    deseq2_qc_cmd.as_task(),
)


@dataclass
class Deseq2QCResult:
    pca_multiqc: File  # <label>.pca.vals_mqc.tsv
    dists_multiqc: File  # <label>.sample.dists_mqc.tsv
    results: Dir  # everything: plots, RData, pca/dists tables, size_factors/, R log


async def deseq2_qc(
    counts: File,
    label: str = "star_salmon",
    prefix: str = "deseq2",
    id_col: int = 1,
    count_col: int = 3,
    sample_suffix: str = "",
    vst: bool = True,
    pca_header: File | None = None,
    clustering_header: File | None = None,
) -> Deseq2QCResult:
    """DESeq2 PCA / sample-distance QC of ``counts``; defaults are upstream's for star_salmon."""
    script = await File.from_local(str(path("deseq2_qc.r")))
    if pca_header is None:
        pca_header = await File.from_local(str(path("deseq2_pca_header.txt")))
    if clustering_header is None:
        clustering_header = await File.from_local(str(path("deseq2_clustering_header.txt")))
    results = await deseq2_qc_cmd(
        script=script,
        counts=counts,
        pca_header=pca_header,
        clustering_header=clustering_header,
        label=label,
        prefix=prefix,
        id_col=id_col,
        count_col=count_col,
        sample_suffix=sample_suffix,
        vst=vst,
    )
    lower = label.lower()
    pca = await results.get_file(f"{lower}.pca.vals_mqc.tsv")
    dists = await results.get_file(f"{lower}.sample.dists_mqc.tsv")
    if pca is None or dists is None:
        raise FileNotFoundError(f"deseq2_qc wrote no {lower}.pca.vals_mqc.tsv / .sample.dists_mqc.tsv")
    return Deseq2QCResult(pca_multiqc=pca, dists_multiqc=dists, results=results)
