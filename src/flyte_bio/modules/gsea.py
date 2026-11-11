"""gsea — gene set enrichment analysis with the Broad's GSEA.

Exposes :func:`gsea`, which runs ``gsea-cli GSEA`` on an expression matrix
(``.gct``) with phenotype labels (``.cls``), comparing ``target`` with
``reference``, over one or more gene set files (``.gmt``, pooled into one
run), optionally collapsing probes to symbols through a ``.chip`` file.
``args`` go to gsea-cli (upstream's tests use ``-rnd_seed 10``).

As upstream, the timestamps GSEA puts in its file names are removed, every
file name is prefixed with ``prefix`` (as is: upstream's names read
``<prefix>gsea_report_for_<level>.tsv``), and per-gene-set reports are
renamed ``gene_sets_<prefix><set>.{tsv,html,png}``. The run's report files
are returned one by one: ``<prefix>.Gsea.rpt`` (GSEA's ``rpt_label`` is the
prefix without a trailing dot), ``<prefix>index.html``,
``<prefix>heat_map_corr_plot.html``,
``<prefix>gsea_report_for_<level>.{tsv,html}`` for each level,
``<prefix>ranked_gene_list*.tsv``, ``<prefix>gene_set_sizes.tsv``, its
plots (``global_es_histogram``, ``heat_map_1``, ``pvalues_vs_nes_plot``,
``ranked_list_corr_2``, ``butterfly_plot``), the ``*_snapshot*.html``
pages and any ``.zip`` report. The per-gene-set files (``gene_sets_*``,
``*enplot*.png``, ``*gset_rnd_es_dist*.png``; several per gene set, so
hundreds per run) come in one ``<prefix>gene_sets.tar.gz``.

GSEA is distributed by the Broad Institute under a BSD-style licence; gene
sets from MSigDB carry their own terms.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

GSEA_IMAGE = "community.wave.seqera.io/library/gsea:4.3.2--a7421d7504fd7c81"
# Upstream: process_single.
DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory="8Gi")

# Upstream's commands, after its inputs are set as shell variables.
GSEA_SH = r"""
# Run GSEA

gsea-cli GSEA \
    -res "$gct" \
    -cls "${cls}#${target}_versus_${reference}" \
    -gmx "$gmx" \
    "${chip_command[@]}" \
    -out . \
    --rpt_label "$rpt_label" \
    $args

# Un-timestamp the outputs for path consistency
mv "$rpt_label".Gsea.*/* .
timestamp=$(cat *.rpt | grep producer_timestamp | awk '{print $2}')

for pattern in _${timestamp} .${timestamp}; do
    find . -name "*${pattern}*" | sed "s|^\./||" | while read -r f; do
        mv "$f" "${f//$pattern/}"
    done
done
sed -i.bak "s/[_\.]$timestamp//g" *.rpt *.html && rm *.bak

# Prefix files that currently lack it
ls -p | grep -v / | grep -v "$prefix" | while read -r f; do
    mv $f ${prefix}${f}
    sed -i.bak "s/$f/${prefix}${f}/g" *.rpt *.html && rm *.bak
done

# Rename files so that they can be properly referenced by the output channels
# Function to rename files based on the given pattern
rename_files() {
    local pattern=$1
    local exclude_patterns=$2
    local extension=$3

    # Find files matching the pattern but not matching the exclusion patterns
    find . -type f -name "$pattern" | while read -r file; do
        # Exclude files based on the provided exclusion patterns
        if ! echo "$file" | grep -qE "$exclude_patterns"; then
            # Rename the file by adding the prefix "gene_sets_"
            mv "$file" "$(dirname "$file")/gene_sets_$(basename "$file")"
        fi
    done
}

# Pattern and exclusion for .tsv files
tsv_pattern="*.tsv"
tsv_exclude="gene_set_size|gsea_report|ranked_gene_list"

# Pattern and exclusion for .html files
html_pattern="*.html"
html_exclude="gsea_report|heat_map_corr_plot|index|pos_snapshot|neg_snapshot"

# Pattern and exclusion for .png files
png_pattern="*.png"
png_exclude="butterfly|enplot|global_es_histogram|gset_rnd_es_dist|heat_map|pvalues_vs_nes_plot|ranked_list_corr"

# Rename .tsv files
rename_files "$tsv_pattern" "$tsv_exclude" ".tsv"

# Rename .html files
rename_files "$html_pattern" "$html_exclude" ".html"

# Rename .png files
rename_files "$png_pattern" "$png_exclude" ".png"
""".strip()

gsea_cmd = shell.create(
    name="gsea",
    image=GSEA_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "gct": File,
        "cls": File,
        "gene_sets": list[File],
        "chip": File | None,
        "reference": str,
        "target": str,
        "prefix": str,
        "args": str,
    },
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        GCT=({inputs.gct})
        CLS=({inputs.cls})
        GENE_SETS=({inputs.gene_sets})
        shopt -s nullglob; CHIP=({inputs.chip}); shopt -u nullglob
        reference={inputs.reference}
        target={inputs.target}
        prefix={inputs.prefix}
        args={inputs.args}
        OUT={outputs.results}
        export HOME=$(mktemp -d)
        W=$(mktemp -d); cd "$W"
        ln -s "${GCT[0]}" .; gct=$(basename "${GCT[0]}")
        ln -s "${CLS[0]}" .; cls=$(basename "${CLS[0]}")
        gmx=""
        for f in "${GENE_SETS[@]}"; do ln -s "$f" .; gmx="${gmx:+$gmx,}$(basename "$f")"; done
        chip_command=()
        if [ ${#CHIP[@]} -gt 0 ]; then
            ln -s "${CHIP[0]}" .; chip_command=(-chip "$(basename "${CHIP[0]}")" -collapse true)
        fi
        rpt_label="${prefix%.}"

GSEA_SH

        shopt -s nullglob
        PER_SET=(gene_sets_* *enplot*.png *gset_rnd_es_dist*.png)
        tar -czf "$OUT/${prefix}gene_sets.tar.gz" "${PER_SET[@]}" || exit 1
        rm -f "${PER_SET[@]}"
        for f in *.rpt *index.html *heat_map_corr_plot.html *gsea_report_for_*.tsv *gsea_report_for_*.html \
            *ranked_gene_list*.tsv *gene_set_sizes.tsv *global_es_histogram.png *heat_map_1.png \
            *pvalues_vs_nes_plot.png *ranked_list_corr_2.png *butterfly_plot.png *_snapshot*.html *.zip; do
            if [ -f "$f" ] && [ ! -L "$f" ]; then mv "$f" "$OUT"/; fi
        done
        true
    """.replace("GSEA_SH", GSEA_SH),
)

env = flyte.TaskEnvironment.from_task("gsea", gsea_cmd.as_task())


@dataclass
class GseaResult:
    rpt: File  # <rpt_label>.Gsea.rpt
    index_html: File  # <prefix>index.html
    heat_map_corr_plot: File  # <prefix>heat_map_corr_plot.html
    report_tsv_ref: File  # <prefix>gsea_report_for_<reference>.tsv
    report_html_ref: File  # <prefix>gsea_report_for_<reference>.html
    report_tsv_target: File  # <prefix>gsea_report_for_<target>.tsv
    report_html_target: File  # <prefix>gsea_report_for_<target>.html
    ranked_gene_list: File  # <prefix>ranked_gene_list*.tsv
    gene_set_sizes: File  # <prefix>gene_set_sizes.tsv
    histogram: File  # <prefix>global_es_histogram.png
    heatmap: File  # <prefix>heat_map_1.png
    pvalues_vs_nes_plot: File  # <prefix>pvalues_vs_nes_plot.png
    ranked_list_corr: File  # <prefix>ranked_list_corr_2.png
    butterfly_plot: File | None  # <prefix>butterfly_plot.png
    snapshots: list[File]  # <prefix>*_snapshot*.html
    archive: File | None  # a .zip report, when args ask for one
    gene_sets: File  # <prefix>gene_sets.tar.gz: gene_sets_*, *enplot*.png, *gset_rnd_es_dist*.png
    results: Dir  # all of the above


async def gsea(
    gct: File,
    cls: File,
    gene_sets: list[File],
    reference: str,
    target: str,
    prefix: str,
    chip: File | None = None,
    args: str = "",
) -> GseaResult:
    """GSEA of ``target`` versus ``reference`` (two of ``cls``'s classes) over ``gene_sets``."""
    results = await gsea_cmd(
        gct=gct, cls=cls, gene_sets=gene_sets, chip=chip, reference=reference, target=target, prefix=prefix, args=args
    )
    files = sorted(await results.list_files(), key=lambda f: f.name)

    def matching(suffix: str) -> list[File]:
        return [f for f in files if f.name.endswith(suffix)]

    def one(suffix: str) -> File:
        found = matching(suffix)
        if len(found) != 1:
            raise FileNotFoundError(f"GSEA wrote {len(found)} *{suffix} files, expected one")
        return found[0]

    ranked = [f for f in files if "ranked_gene_list" in f.name and f.name.endswith(".tsv")]
    if len(ranked) != 1:
        raise FileNotFoundError(f"GSEA wrote {len(ranked)} ranked gene lists, expected one")
    butterfly = matching("butterfly_plot.png")
    archives = matching(".zip")
    return GseaResult(
        rpt=one(".rpt"),
        index_html=one("index.html"),
        heat_map_corr_plot=one("heat_map_corr_plot.html"),
        report_tsv_ref=one(f"gsea_report_for_{reference}.tsv"),
        report_html_ref=one(f"gsea_report_for_{reference}.html"),
        report_tsv_target=one(f"gsea_report_for_{target}.tsv"),
        report_html_target=one(f"gsea_report_for_{target}.html"),
        ranked_gene_list=ranked[0],
        gene_set_sizes=one("gene_set_sizes.tsv"),
        histogram=one("global_es_histogram.png"),
        heatmap=one("heat_map_1.png"),
        pvalues_vs_nes_plot=one("pvalues_vs_nes_plot.png"),
        ranked_list_corr=one("ranked_list_corr_2.png"),
        butterfly_plot=butterfly[0] if butterfly else None,
        snapshots=[f for f in files if "_snapshot" in f.name and f.name.endswith(".html")],
        archive=archives[0] if archives else None,
        gene_sets=one("gene_sets.tar.gz"),
        results=results,
    )
