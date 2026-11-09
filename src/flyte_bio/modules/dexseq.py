"""dexseq — differential exon usage (DEU) and transcript usage (DTU) with DEXSeq.

Exposes:

- :func:`dexseq_prepareannotation` — flatten a GTF into DEXSeq's exonic-part
  GFF with DEXSeq's ``dexseq_prepare_annotation.py``: ``<prefix>.gff``.
  Overlapping genes are merged into aggregate genes (IDs joined with ``+``)
  unless ``args`` has ``-r no``. The script iterates over Python sets, so the
  order of the aggregate gene IDs depends on Python's hash seed: pass
  ``random_seed`` (as upstream, ``PYTHONHASHSEED``) for a reproducible GFF.
- :func:`dexseq_count` — count a BAM's reads per exonic part with DEXSeq's
  ``dexseq_count.py`` (``--format bam --order pos``, paired unless
  ``single_end``): ``<prefix>.txt``. As upstream, the quotes HTSeq 2 puts
  around the feature IDs are removed, so ``DEXSeqDataSetFromHTSeq`` can read
  the table.
- :func:`dexseq_exon` — test for DEU between two levels of a sample-sheet
  column, from per-sample exon counts and the flattened GFF.
- :func:`dexseq_dtu` — test for DTU between two levels, from one table of
  transcript counts (feature and gene ID columns plus one column per
  sample; estimated counts are rounded).

:func:`dexseq_exon` and :func:`dexseq_dtu` write ``<prefix>.dexseq.results.tsv``,
``<prefix>.dexseq.gene_qvalues.tsv``, the ``.dataset.rds``, ``.results.rds``
and ``.gene_qvalues.rds`` R objects, ``<prefix>.R_sessionInfo.log`` and a
``versions.yml``; :func:`dexseq_exon` also plots the ``n_plot`` (default 10)
genes with the smallest p-values to ``<prefix>.dexseq.plots.pdf``. Their R
code is upstream's templates, unchanged but for taking the templated values
(prefix, inputs, contrast and CPUs) as command-line arguments. ``args`` takes
the templates' ``--option value`` overrides (``--sample_id_col``,
``--count_suffix``, ``--n_plot``, ``--feature_id_col``, ``--gene_id_col``, ...).
As upstream, :func:`dexseq_exon` reads the counts as ``counts/<sample><count
suffix>``.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

PREPAREANNOTATION_IMAGE = "community.wave.seqera.io/library/bioconductor-dexseq_htseq_python:dfa2c806cf6767c6"
COUNT_IMAGE = "community.wave.seqera.io/library/bioconductor-dexseq_htseq:4266bf085c9883ac"
DEXSEQ_IMAGE = "community.wave.seqera.io/library/bioconductor-dexseq:1.56.0--548bf2f0349db29a"

# Upstream: process_single for the Python scripts, process_medium for the R analyses.
SINGLE_RESOURCES = flyte.Resources(cpu=1, memory="6Gi")
ANALYSIS_RESOURCES = flyte.Resources(cpu=4, memory="12Gi")

# Shared by both R templates: upstream's helper functions.
R_HELPERS = r"""
#' Parse out options from a string without recourse to optparse
#'
#' @param x Long-form argument list like --opt1 val1 --opt2 val2
#'
#' @return named list of options and values similar to optparse

parse_args <- function(x) {
    args_list <- unlist(strsplit(x, ' ?--')[[1]])[-1]
    args_vals <- lapply(args_list, function(x) scan(text = x, what = 'character', quiet = TRUE))

    # Ensure the option vectors are length 2 (key/ value) to catch empty ones
    args_vals <- lapply(args_vals, function(z) { length(z) <- 2; z })

    parsed_args <- structure(lapply(args_vals, function(x) x[2]), names = lapply(args_vals, function(x) x[1]))
    parsed_args[! is.na(parsed_args)]
}

#' Flexibly read CSV or TSV files
#'
#' @param file Input file
#'
#' @return Data frame

read_delim_flexible <- function(file) {
    ext <- tolower(tail(strsplit(basename(file), split = "\\.")[[1]], 1))

    if (ext == "tsv" || ext == "txt") {
        separator <- "\t"
    } else if (ext == "csv") {
        separator <- ","
    } else {
        stop(paste("Unknown separator for", ext))
    }

    read.delim(file, sep = separator, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
}

#' Keep the result columns that can be written to a flat table
#'
#' @param x DEXSeqResults object
#'
#' @return DEXSeqResults object without the list columns

keep_valid_columns <- function(x) {
    keep <- vapply(x, function(col) is.numeric(col) || is.character(col), FUN.VALUE = logical(1))
    x[, keep, drop = FALSE]
}
""".strip()

# Upstream's sample-sheet reading, common to both templates.
R_READ_SAMPLES = r"""
sample_sheet <- read_delim_flexible(opt[['sample_file']])

for (column in c(opt[['sample_id_col']], opt[['contrast_variable']])) {
    if (! column %in% colnames(sample_sheet)) {
        stop(paste0("Column '", column, "' is not in the sample sheet"))
    }
}

# The sample sheet can have several rows per sample, e.g. for sequencing runs

sample_sheet <- sample_sheet[! duplicated(sample_sheet[[opt[['sample_id_col']]]]), , drop = FALSE]

levels <- c(opt[['reference_level']], opt[['target_level']])
missing_levels <- setdiff(levels, sample_sheet[[opt[['contrast_variable']]]])

if (length(missing_levels) > 0) {
    stop(paste0(
        "Levels missing from the '", opt[['contrast_variable']], "' column of the sample sheet: ",
        paste(missing_levels, collapse = ", ")
    ))
}

# Keep the samples of the two levels being compared, with the reference first so
# that the fold changes are those of the target over the reference

sample_sheet <- sample_sheet[sample_sheet[[opt[['contrast_variable']]]] %in% levels, , drop = FALSE]

sample_data <- data.frame(
    sample = as.character(sample_sheet[[opt[['sample_id_col']]]]),
    condition = factor(sample_sheet[[opt[['contrast_variable']]]], levels = levels)
)
""".strip()

# Upstream's test fitting and outputs, common to both templates.
R_FIT_AND_WRITE = r"""
dxd <- estimateSizeFactors(dxd)
dxd <- estimateDispersions(dxd, BPPARAM = SerialParam())
dxd <- testForDEU(dxd, BPPARAM = BPPARAM)
dxd <- estimateExonFoldChanges(dxd, BPPARAM = BPPARAM)

results <- DEXSeqResults(dxd)
gene_qvalues <- perGeneQValue(results)

prefix <- opt[['output_prefix']]

saveRDS(dxd, file = paste0(prefix, '.dexseq.dataset.rds'))
saveRDS(results, file = paste0(prefix, '.dexseq.results.rds'))
saveRDS(gene_qvalues, file = paste0(prefix, '.dexseq.gene_qvalues.rds'))

write.table(
    as.data.frame(keep_valid_columns(results)),
    file = paste0(prefix, '.dexseq.results.tsv'),
    sep = '\t',
    quote = FALSE,
    row.names = FALSE
)

write.table(
    data.frame(groupID = names(gene_qvalues), padj = unname(gene_qvalues)),
    file = paste0(prefix, '.dexseq.gene_qvalues.tsv'),
    sep = '\t',
    quote = FALSE,
    row.names = FALSE
)
""".strip()


def r_versions(process: str) -> str:
    """Upstream's versions.yml writer, under upstream's process name."""
    return rf"""
sink(paste0(prefix, '.R_sessionInfo.log'))
print(sessionInfo())
sink()

writeLines(
    c(
        '"{process}":',
        paste('    r-base:', paste(R.version[['major']], R.version[['minor']], sep = '.')),
        paste('    bioconductor-dexseq:', as.character(packageVersion('DEXSeq')))
    ),
    'versions.yml'
)
""".strip()


# Arguments: prefix, sample sheet, gff, contrast variable, reference, target, cores, overrides.
DEXSEQ_EXON_R = "\n\n".join(
    [
        "cli <- commandArgs(trailingOnly = TRUE)",
        R_HELPERS,
        r"""
#' Plot the exon usage of the genes with the smallest p-values
#'
#' @param x DEXSeqResults object
#' @param file Path of the PDF file to write
#' @param n Number of genes to plot

save_plot_dexseq <- function(x, file, n) {
    x_valid <- x[is.finite(x[['pvalue']]) & is.finite(x[['padj']]), , drop = FALSE]

    if (nrow(x_valid) == 0) {
        warning("No valid results to plot.")
        return(invisible(NULL))
    }

    n <- min(n, nrow(x_valid))
    ids <- x_valid[['groupID']][order(x_valid[['pvalue']])[seq_len(n)]]

    pdf(file = file, width = 11.7, height = 8.3)
    on.exit(dev.off())

    for (i in ids) {
        tryCatch(
            plotDEXSeq(x, geneID = i, displayTranscripts = TRUE),
            error = function(e) warning(paste("Error plotting gene", i, ":", conditionMessage(e)))
        )
    }

    invisible(NULL)
}

opt <- list(
    output_prefix = cli[1],
    sample_file = cli[2],
    count_dir = 'counts',
    annotation_file = cli[3],
    contrast_variable = cli[4],
    reference_level = cli[5],
    target_level = cli[6],
    sample_id_col = 'sample',
    count_suffix = '.txt',
    n_plot = 10,
    cores = as.integer(cli[7])
)
opt_types <- lapply(opt, class)

# Apply parameter overrides

args_opt <- parse_args(cli[8])
for (ao in names(args_opt)) {
    if (! ao %in% names(opt)) {
        stop(paste("Invalid option:", ao))
    }
    opt[[ao]] <- as(args_opt[[ao]], opt_types[[ao]])
}

if (is.na(opt[['n_plot']]) || opt[['n_plot']] < 0) {
    stop("Option --n_plot must be a non-negative number")
}

if (opt[['reference_level']] == opt[['target_level']]) {
    stop("The reference and target levels must differ")
}

library(BiocParallel)
library(DEXSeq)

if (! is.na(opt[['cores']]) && opt[['cores']] > 1) {
    BPPARAM <- MulticoreParam(workers = opt[['cores']])
} else {
    BPPARAM <- SerialParam()
}
""".strip(),
        R_READ_SAMPLES,
        r"""
count_files <- file.path(opt[['count_dir']], paste0(sample_data[['sample']], opt[['count_suffix']]))
missing_files <- count_files[! file.exists(count_files)]

if (length(missing_files) > 0) {
    stop(paste("Missing count files:", paste(basename(missing_files), collapse = ", ")))
}

dxd <- DEXSeqDataSetFromHTSeq(
    countfiles = count_files,
    sampleData = sample_data,
    design = ~ sample + exon + condition:exon,
    flattenedfile = opt[['annotation_file']]
)
""".strip(),
        R_FIT_AND_WRITE,
        r"""
if (opt[['n_plot']] > 0) {
    save_plot_dexseq(results, file = paste0(prefix, '.dexseq.plots.pdf'), n = opt[['n_plot']])
}
""".strip(),
        r_versions("DEXSEQ_EXON"),
    ]
)

# Arguments: prefix, sample sheet, count table, contrast variable, reference, target, cores, overrides.
DEXSEQ_DTU_R = "\n\n".join(
    [
        "cli <- commandArgs(trailingOnly = TRUE)",
        R_HELPERS,
        r"""
opt <- list(
    output_prefix = cli[1],
    sample_file = cli[2],
    count_file = cli[3],
    contrast_variable = cli[4],
    reference_level = cli[5],
    target_level = cli[6],
    sample_id_col = 'sample',
    feature_id_col = 'feature_id',
    gene_id_col = 'gene_id',
    cores = as.integer(cli[7])
)
opt_types <- lapply(opt, class)

# Apply parameter overrides

args_opt <- parse_args(cli[8])
for (ao in names(args_opt)) {
    if (! ao %in% names(opt)) {
        stop(paste("Invalid option:", ao))
    }
    opt[[ao]] <- as(args_opt[[ao]], opt_types[[ao]])
}

if (opt[['reference_level']] == opt[['target_level']]) {
    stop("The reference and target levels must differ")
}

library(BiocParallel)
library(DEXSeq)

if (! is.na(opt[['cores']]) && opt[['cores']] > 1) {
    BPPARAM <- MulticoreParam(workers = opt[['cores']])
} else {
    BPPARAM <- SerialParam()
}
""".strip(),
        R_READ_SAMPLES,
        r"""
count_table <- read_delim_flexible(opt[['count_file']])

for (column in c(opt[['feature_id_col']], opt[['gene_id_col']])) {
    if (! column %in% colnames(count_table)) {
        stop(paste0("Column '", column, "' is not in the count table"))
    }
}

missing_samples <- setdiff(sample_data[['sample']], colnames(count_table))

if (length(missing_samples) > 0) {
    stop(paste("Samples missing from the count table:", paste(missing_samples, collapse = ", ")))
}

# DEXSeq models integer counts, so estimated counts such as those of tximport are rounded

counts <- round(as.matrix(count_table[, sample_data[['sample']], drop = FALSE]))
rownames(counts) <- count_table[[opt[['feature_id_col']]]]

dxd <- DEXSeqDataSet(
    countData = counts,
    sampleData = sample_data,
    design = ~ sample + exon + condition:exon,
    featureID = count_table[[opt[['feature_id_col']]]],
    groupID = count_table[[opt[['gene_id_col']]]]
)
""".strip(),
        R_FIT_AND_WRITE,
        r_versions("DEXSEQ_DTU"),
    ]
)


def run_r(script: str, name: str, arguments: str) -> str:
    """Shell lines that write ``script`` to ``<name>.R`` and run it, putting R's error in the termination log."""
    return rf"""
        cat > {name}.R <<'END_OF_SCRIPT'
{script}
END_OF_SCRIPT
        if ! Rscript {name}.R {arguments} 2> {name}.err; then
            cat {name}.err >&2
            {{ tail -c 2048 {name}.err > /dev/termination-log; }} 2> /dev/null
            exit 1
        fi
        cat {name}.err >&2"""


dexseq_prepareannotation_cmd = shell.create(
    name="dexseq_prepareannotation",
    image=PREPAREANNOTATION_IMAGE,
    resources=SINGLE_RESOURCES,
    inputs={"gtf": File, "prefix": str, "random_seed": str, "args": str},
    defaults={"random_seed": "", "args": ""},
    outputs={"results": Dir},
    script=r"""
        GTF=({inputs.gtf})
        P={inputs.prefix}
        SEED={inputs.random_seed}
        ARGS={inputs.args}
        cd {outputs.results}
        DEXSEQ_SCRIPTS=$(Rscript -e 'cat(system.file("python_scripts", package = "DEXSeq"))')
        if [ -n "$SEED" ]; then export PYTHONHASHSEED="$SEED"; fi
        python "$DEXSEQ_SCRIPTS/dexseq_prepare_annotation.py" $ARGS "${GTF[0]}" "$P.gff" || exit 1
    """,
)

dexseq_count_cmd = shell.create(
    name="dexseq_count",
    image=COUNT_IMAGE,
    resources=SINGLE_RESOURCES,
    inputs={"bam": File, "gff": File, "prefix": str, "single_end": bool, "args": str},
    defaults={"single_end": False, "args": ""},
    outputs={"results": Dir},
    script=r"""
        BAM=({inputs.bam}); GFF=({inputs.gff})
        P={inputs.prefix}
        ARGS={inputs.args}
        PAIRED=yes
        if [ {inputs.single_end} = true ]; then PAIRED=no; fi
        cd {outputs.results}
        DEXSEQ_SCRIPTS=$(Rscript -e 'cat(system.file("python_scripts", package = "DEXSeq"))')
        python "$DEXSEQ_SCRIPTS/dexseq_count.py" --paired $PAIRED --format bam --order pos $ARGS \
            "${GFF[0]}" "${BAM[0]}" "$P.txt" || exit 1
        sed -i 's/"//g' "$P.txt"
    """,
)

dexseq_exon_cmd = shell.create(
    name="dexseq_exon",
    image=DEXSEQ_IMAGE,
    resources=ANALYSIS_RESOURCES,
    inputs={
        "samplesheet": File,
        "counts": list[File],
        "gff": File,
        "contrast_variable": str,
        "reference": str,
        "target": str,
        "prefix": str,
        "args": str,
    },
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        SHEET=({{inputs.samplesheet}}); COUNTS=({{inputs.counts}}); GFF=({{inputs.gff}})
        CONTRAST={{inputs.contrast_variable}}; REFERENCE={{inputs.reference}}; TARGET={{inputs.target}}
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        mkdir counts
        for F in "${{COUNTS[@]}}"; do ln -s "$F" counts/; done
        ln -s "${{SHEET[0]}}" .; ln -s "${{GFF[0]}}" .
        SHEET_NAME=$(basename "${{SHEET[0]}}"); GFF_NAME=$(basename "${{GFF[0]}}")
        {
        run_r(
            DEXSEQ_EXON_R,
            "dexseq_exon",
            '"$P" "$SHEET_NAME" "$GFF_NAME" "$CONTRAST" "$REFERENCE" "$TARGET" "$TASK_CPUS" "$ARGS"',
        )
    }
        shopt -s nullglob
        mv ./*.dexseq.* ./*.R_sessionInfo.log versions.yml "$OUT"/
    """,
)

dexseq_dtu_cmd = shell.create(
    name="dexseq_dtu",
    image=DEXSEQ_IMAGE,
    resources=ANALYSIS_RESOURCES,
    inputs={
        "samplesheet": File,
        "counts": File,
        "contrast_variable": str,
        "reference": str,
        "target": str,
        "prefix": str,
        "args": str,
    },
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        SHEET=({{inputs.samplesheet}}); COUNTS=({{inputs.counts}})
        CONTRAST={{inputs.contrast_variable}}; REFERENCE={{inputs.reference}}; TARGET={{inputs.target}}
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        ln -s "${{SHEET[0]}}" .; ln -s "${{COUNTS[0]}}" .
        SHEET_NAME=$(basename "${{SHEET[0]}}"); COUNTS_NAME=$(basename "${{COUNTS[0]}}")
        {
        run_r(
            DEXSEQ_DTU_R,
            "dexseq_dtu",
            '"$P" "$SHEET_NAME" "$COUNTS_NAME" "$CONTRAST" "$REFERENCE" "$TARGET" "$TASK_CPUS" "$ARGS"',
        )
    }
        mv ./*.dexseq.* ./*.R_sessionInfo.log versions.yml "$OUT"/
    """,
)

# One environment per image.
dexseq_prepareannotation_env = flyte.TaskEnvironment.from_task(
    "dexseq_prepareannotation", dexseq_prepareannotation_cmd.as_task()
)
dexseq_count_env = flyte.TaskEnvironment.from_task("dexseq_count", dexseq_count_cmd.as_task())
dexseq_analysis_env = flyte.TaskEnvironment.from_task(
    "dexseq_analysis", dexseq_exon_cmd.as_task(), dexseq_dtu_cmd.as_task()
)
env = flyte.TaskEnvironment(
    name="dexseq",
    depends_on=[dexseq_prepareannotation_env, dexseq_count_env, dexseq_analysis_env],
)


async def pick(results: Dir, name: str, tool: str) -> File:
    found = await results.get_file(name)
    if found is None:
        raise FileNotFoundError(f"{tool} wrote no {name}")
    return found


async def dexseq_prepareannotation(
    gtf: File, prefix: str = "annotation", random_seed: int | None = None, args: str = ""
) -> File:
    """Flatten ``gtf`` into DEXSeq's exonic-part GFF, ``<prefix>.gff``.

    ``random_seed`` sets Python's hash seed, which fixes the order of the
    aggregate gene IDs; ``args`` takes the script's options (e.g. ``-r no``).
    """
    seed = "" if random_seed is None else str(random_seed)
    results = await dexseq_prepareannotation_cmd(gtf=gtf, prefix=prefix, random_seed=seed, args=args)
    return await pick(results, f"{prefix}.gff", "dexseq_prepare_annotation")


async def dexseq_count(bam: File, gff: File, prefix: str, single_end: bool = False, args: str = "") -> File:
    """Count ``bam``'s reads per exonic part of ``gff``: ``<prefix>.txt``.

    ``args`` takes the script's options, e.g. ``--stranded no`` (its default
    is ``yes``) or ``--minaqual 10``.
    """
    results = await dexseq_count_cmd(bam=bam, gff=gff, prefix=prefix, single_end=single_end, args=args)
    return await pick(results, f"{prefix}.txt", "dexseq_count")


@dataclass
class DexseqResult:
    results: File  # <prefix>.dexseq.results.tsv
    gene_qvalues: File  # <prefix>.dexseq.gene_qvalues.tsv
    dataset_rds: File
    results_rds: File
    gene_qvalues_rds: File
    plots: File | None  # <prefix>.dexseq.plots.pdf (dexseq_exon, unless --n_plot 0 or nothing to plot)
    session_info: File  # <prefix>.R_sessionInfo.log
    versions: File
    outputs: Dir  # everything


async def read_dexseq_results(outputs: Dir, prefix: str, tool: str) -> DexseqResult:
    return DexseqResult(
        results=await pick(outputs, f"{prefix}.dexseq.results.tsv", tool),
        gene_qvalues=await pick(outputs, f"{prefix}.dexseq.gene_qvalues.tsv", tool),
        dataset_rds=await pick(outputs, f"{prefix}.dexseq.dataset.rds", tool),
        results_rds=await pick(outputs, f"{prefix}.dexseq.results.rds", tool),
        gene_qvalues_rds=await pick(outputs, f"{prefix}.dexseq.gene_qvalues.rds", tool),
        plots=await outputs.get_file(f"{prefix}.dexseq.plots.pdf"),
        session_info=await pick(outputs, f"{prefix}.R_sessionInfo.log", tool),
        versions=await pick(outputs, "versions.yml", tool),
        outputs=outputs,
    )


async def dexseq_exon(
    samplesheet: File,
    counts: list[File],
    gff: File,
    contrast_variable: str,
    reference: str,
    target: str,
    prefix: str,
    args: str = "",
) -> DexseqResult:
    """Test for differential exon usage of ``target`` over ``reference`` (levels of ``contrast_variable``).

    ``samplesheet`` (CSV or TSV) has a ``sample`` column; ``counts`` are the
    samples' exon counts (:func:`dexseq_count`), named ``<sample>.txt``;
    ``gff`` is the flattened annotation they were counted against.
    """
    outputs = await dexseq_exon_cmd(
        samplesheet=samplesheet,
        counts=counts,
        gff=gff,
        contrast_variable=contrast_variable,
        reference=reference,
        target=target,
        prefix=prefix,
        args=args,
    )
    return await read_dexseq_results(outputs, prefix, "dexseq_exon")


async def dexseq_dtu(
    samplesheet: File,
    counts: File,
    contrast_variable: str,
    reference: str,
    target: str,
    prefix: str,
    args: str = "",
) -> DexseqResult:
    """Test for differential transcript usage of ``target`` over ``reference`` (levels of ``contrast_variable``).

    ``counts`` (CSV or TSV) has a feature ID column (``feature_id``, or
    ``--feature_id_col`` in ``args``), a ``gene_id`` column and one column
    per sample of ``samplesheet``.
    """
    outputs = await dexseq_dtu_cmd(
        samplesheet=samplesheet,
        counts=counts,
        contrast_variable=contrast_variable,
        reference=reference,
        target=target,
        prefix=prefix,
        args=args,
    )
    return await read_dexseq_results(outputs, prefix, "dexseq_dtu")
