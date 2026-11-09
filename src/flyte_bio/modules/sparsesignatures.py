"""sparsesignatures — SparseSignatures, de novo mutational signatures with an SBS5 background.

Exposes :func:`sparsesignatures`, which pools the joint mutation tables in
``tsv_join`` (``Indiv``, ``chr``, ``from``, ``to``, ``ref``, ``alt``, ``NV``
…; mutations with ``NV`` 0 are dropped), counts trinucleotide contexts against
``genome`` (``GRCh37``: BSgenome hs37d5; ``GRCh38``: BSgenome UCSC hg38),
picks the number of signatures and the sparsity by cross-validation, and fits
the best configuration with NMF Lasso. It writes, with a ``versions.yml``:

- ``<prefix>_mut_counts.rds``: the patient-by-context count matrix;
- ``<prefix>_cv_means_mse.rds`` and ``<prefix>_best_params_config.rds``: the
  cross-validation's mean squared errors and its best K and lambda;
- ``<prefix>_nmf_Lasso_out.rds``: the fit;
- ``<prefix>_plot_all.{rds,pdf}``: exposures and signatures.

``args`` is upstream's ``ext.args`` map as Groovy prints it,
``[K:2:3, nmf_runs:3, seed:42, …]``; the template's defaults cover every key
(seed 42 unless set). Upstream's default ``num_processes:all`` would start a
worker per core of the node, so when ``args`` doesn't set ``num_processes``
the task passes its own CPU count instead.

The R code is upstream's template, written out by an unquoted bash heredoc:
Groovy's ``\\$`` and ``\\\\`` escapes mean the same to bash, so only the
dotted placeholders are renamed (``$task.ext.args`` → ``${task_ext_args}``).
The tables are staged as ``1.tsv``, ``2.tsv`` … (upstream stages them as ``*.tsv``).
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

SPARSESIGNATURES_IMAGE = (
    "quay.io/biocontainers/mulled-v2-d410175e2fbd9c47aa685bb5dfb87cfad76d408b"
    ":a995bb98b7122825523ffed7ae131cb006e56cbe-0"
)
# The versions export loads both BSgenome packages.
DEFAULT_RESOURCES = flyte.Resources(cpu=2, memory="12Gi")

# Upstream's templates/main_script.R, for an unquoted heredoc (see the module docstring).
SPARSESIGNATURES_TEMPLATE = r"""
#!/usr/bin/env Rscript

pkgs <- c("SparseSignatures", "dplyr", "tidyr", "tibble", "purrr", "stringr", "ggplot2", "patchwork")
sapply(pkgs, require, character.only = TRUE)

parse_args <- function(x) {
  # Remove brackets
  x = gsub("\\\\[","",x)
  x = gsub("\\\\]","",x)

  # Split into key:value pairs
  args_list = unlist(strsplit(x, ", (?=[^)]*(?:\\\\(|\$))", perl=TRUE))

  # Convert to tibble and split key:value
  tibble(arg = args_list) %>%
    separate(arg, into = c("key", "value"), sep = ":", extra = "merge", fill = "right") %>%
    dplyr::mutate(across(everything(), ~str_trim(.))) %>%
    dplyr::filter(!is.na(key) & !is.na(value)) %>%
    tibble::deframe()
}

opt = list(
    prefix = ifelse('${task_ext_prefix}' == 'null', '${meta_id}', '${task_ext_prefix}'),
    genome = "${genome}",
    K = "2:10",
    nmf_runs = "10",
    iterations = "30",
    max_iterations_lasso = "10000",
    num_processes = "all",
    cross_validation_entries = "0.01",
    cross_validation_repetitions = "50",
    cross_validation_iterations = "5",
    lambda_values_alpha = "c(0.00, 0.01, 0.05, 0.10)",
    lambda_values_beta = "c(0.01, 0.05, 0.1, 0.2)",
    lambda_rate_alpha = "0",
    seed = "NULL"
)

args_opt = parse_args('${task_ext_args}')
for ( ao in names(args_opt)) opt[[ao]] = args_opt[[ao]]


# Script #

num_procs_string <- opt[["num_processes"]]

if (is.null(num_procs_string)) {
    stop("Missing required option: num_processes")
} else if (num_procs_string == "all") {
    n_procs <- Inf
} else {
    n_procs <- as.integer(opt[["num_processes"]])
}


patients_tsv = strsplit("$tsv_join", " ")[[1]]
tables = lapply(patients_tsv, FUN = function(p_table){
    read.delim(p_table, sep = "\\t", header=T) %>%
        mutate(across(everything(), as.character))
}
)
# filter to keep only mutations with NV different from 0
multisample_table = dplyr::bind_rows(tables) %>% dplyr::filter(NV != 0)

#Extract input data information
input_data = multisample_table[,c("Indiv","chr","from","to","ref","alt")]
input_data = setNames(input_data, c("sample","chrom","start","end","ref","alt"))
input_data[["end"]] = input_data[["start"]]
input_data = input_data %>% mutate(start = as.integer(start), end = as.integer(end))

#Generate the patient vs mutation count matrix from mutation data
#Load reference human-genome specification.
#The user must select, among the available choices, the reference genome consistent with the mutation dataset.

load_genome = function(genome, input_data) {
    if (genome == "GRCh37") {
        library(BSgenome.Hsapiens.1000genomes.hs37d5)
        bsg = BSgenome.Hsapiens.1000genomes.hs37d5::hs37d5

	if (any(grepl("^chr", input_data[["chrom"]]))) {
            input_data <- input_data %>% mutate(chrom = str_remove(chrom, "^chr"))
        }

    } else if (genome == "GRCh38") {
        library(BSgenome.Hsapiens.UCSC.hg38)
        bsg = BSgenome.Hsapiens.UCSC.hg38

	if (all(!grepl("^chr", input_data[["chrom"]]))) {
            input_data <- input_data %>% mutate(chrom = paste0("chr", chrom))
        }

    }
    return(list(bsg = bsg, input_data = input_data))
}

data_list = load_genome(opt[["genome"]], input_data)
bsg = data_list[["bsg"]]
input_data = data_list[["input_data"]]

mut_counts = SparseSignatures::import.trinucleotides.counts(data=input_data, reference=bsg)

saveRDS(object = mut_counts, file = paste0(opt[["prefix"]], "_mut_counts.rds"))

# Load a reference SBS5 background signature from COSMIC
data(background)


seed_val <- if (!is.null(opt[["seed"]]) && opt[["seed"]] != "NULL") {
    as.integer(opt[["seed"]])
} else {
    42  # fallback default
}
set.seed(seed_val)

# Estimate the initial values of beta
starting_betas = SparseSignatures::startingBetaEstimation(x = mut_counts,
                                                        K = eval(parse(text=opt[["K"]])),
                                                        background_signature = background,
                                                        nmf_runs = as.integer(opt[["nmf_runs"]]),
                                                        seed = seed_val)

# Find the optimal number of signatures and sparsity level: rely on cross-validation
# higher number of CV repetitions corresponds to more accurate parameter estimates

cv_out = SparseSignatures::nmfLassoCV(
  x = mut_counts,
  K = eval(parse(text=opt[["K"]])),
  starting_beta = starting_betas,
  background_signature = background,
  nmf_runs = as.integer(opt[["nmf_runs"]]),
  lambda_values_alpha = eval(parse(text=opt[["lambda_values_alpha"]])),
  lambda_values_beta = eval(parse(text=opt[["lambda_values_beta"]])),
  cross_validation_entries = as.numeric(opt[["cross_validation_entries"]]),
  cross_validation_iterations = as.integer(opt[["cross_validation_iterations"]]),
  cross_validation_repetitions = as.integer(opt[["cross_validation_repetitions"]]),
  iterations = as.integer(opt[["iterations"]]),
  max_iterations_lasso = as.integer(opt[["max_iterations_lasso"]]),
  num_processes = n_procs,
  seed = seed_val
)

str(cv_out)


# Analyze the mean squared error results averaging over cross-validation repetitions
cv_mses = cv_out[["grid_search_mse"]][1,,]
cv_means_mse = matrix(sapply(cv_mses, FUN = mean),
                    nrow = dim(cv_mses)[1]
)

dimnames(cv_means_mse) = dimnames(cv_mses)

# Find the combination of parameters that yields the lowest MSE

print("CV MEAN MSE")
print(cv_means_mse)

print("MIN MEANS MSE")
print(min(cv_means_mse, na.rm = TRUE))

min_ii = which(cv_means_mse == min(cv_means_mse, na.rm = TRUE), arr.ind = TRUE)
min_Lambda_beta = rownames(cv_means_mse)[min_ii[1]]
min_Lambda_beta = as.numeric(gsub("_Lambda_Beta", "", min_Lambda_beta))
min_K = colnames(cv_means_mse)[min_ii[2]]
min_K = as.numeric(gsub("_Signatures", "", min_K))
best_params_config = data.frame(min_K, min_Lambda_beta)

saveRDS(object = cv_means_mse, file = paste0(opt[["prefix"]], "_cv_means_mse.rds"))
saveRDS(object = best_params_config, file = paste0(opt[["prefix"]], "_best_params_config.rds"))

# Discovering the signatures within the dataset: NMF Lasso
# Compute the signatures for the best configuration.

print(paste("MIN K =", min_K))

nmf_Lasso_out = SparseSignatures::nmfLasso(
  x = mut_counts,
  K = min_K,
  background_signature = background,
  lambda_rate_alpha = eval(parse(text=opt[["lambda_rate_alpha"]])),
  lambda_rate_beta = min_Lambda_beta,
  iterations = as.integer(opt[["iterations"]]),
  max_iterations_lasso = as.integer(opt[["max_iterations_lasso"]]),
  seed = seed_val
)

saveRDS(object = nmf_Lasso_out, file =  paste0(opt[["prefix"]], "_nmf_Lasso_out.rds"))

# Signature visualization
signatures = nmf_Lasso_out[["beta"]]
plot_signatures = SparseSignatures::signatures.plot(beta=signatures, xlabels=FALSE)

plot_exposure = nmf_Lasso_out[["alpha"]] %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var="PatientID") %>%
    tidyr::pivot_longer(cols=!"PatientID", names_to="Signatures", values_to="Exposures") %>%

    ggplot() +
    geom_bar(aes(x=PatientID, y=Exposures, fill=Signatures),
            position="stack", stat="identity") +
    theme(axis.text.x=element_text(angle=90,hjust=1),
            panel.background=element_blank(),
            axis.line=element_line(colour="black"))

plt_all = patchwork::wrap_plots(plot_exposure, plot_signatures, ncol=2) + patchwork::plot_annotation(title = "${meta_id}")
ggplot2::ggsave(plot = plt_all, filename = paste0(opt[["prefix"]], "_plot_all.pdf"), width = 210, height = 297, units="mm", dpi = 200, device = "pdf")
saveRDS(object = plt_all, file = paste0(opt[["prefix"]], "_plot_all.rds"))

# version export
f = file("versions.yml","w")
SparseSignatures_version = sessionInfo()\$otherPkgs\$SparseSignatures\$Version
dplyr_version = sessionInfo()\$otherPkgs\$dplyr\$Version
ggplot2_version = sessionInfo()\$otherPkgs\$ggplot2\$Version
patchwork_version = sessionInfo()\$otherPkgs\$patchwork\$Version
writeLines(paste0('"', "${task_process}", '"', ":"), f)
writeLines(paste("    SparseSignatures:", SparseSignatures_version), f)
writeLines(paste("    dplyr:", dplyr_version), f)
writeLines(paste("    ggplot2:", ggplot2_version), f)
writeLines(paste("    patchwork:", patchwork_version), f)
library(BSgenome.Hsapiens.1000genomes.hs37d5)
bsg_hs37_version = sessionInfo()\$otherPkgs\$BSgenome.Hsapiens.1000genomes.hs37d5\$Version
writeLines(paste("    BSgenome.Hsapiens.1000genomes.hs37d5:", bsg_hs37_version), f)
library(BSgenome.Hsapiens.UCSC.hg38)
bsg_hg38_version = sessionInfo()\$otherPkgs\$BSgenome.Hsapiens.UCSC.hg38\$Version
writeLines(paste("    BSgenome.Hsapiens.UCSC.hg38:", bsg_hg38_version), f)
close(f)
"""  # noqa: E501

sparsesignatures_cmd = shell.create(
    name="sparsesignatures",
    image=SPARSESIGNATURES_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"tsv_join": list[File], "genome": str, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        TSV=({{inputs.tsv_join}})
        genome={{inputs.genome}}
        meta_id={{inputs.prefix}}
        task_ext_prefix=$meta_id
        task_ext_args={{inputs.args}}
        if [[ "$task_ext_args" != *num_processes:* ]]; then
            if [ -z "$task_ext_args" ]; then
                task_ext_args="[num_processes:$TASK_CPUS]"
            else
                task_ext_args="${{task_ext_args%]}}, num_processes:$TASK_CPUS]"
            fi
        fi
        task_process=SPARSE_SIGNATURES
        W=$(mktemp -d); cd "$W"
        tsv_join=""
        for I in "${{!TSV[@]}}"; do
            ln -s "${{TSV[$I]}}" "$((I + 1)).tsv"
            tsv_join="$tsv_join${{tsv_join:+ }}$((I + 1)).tsv"
        done
        cat > main_script.R <<END_OF_TEMPLATE"""
    + SPARSESIGNATURES_TEMPLATE
    + r"""END_OF_TEMPLATE
        Rscript main_script.R || exit 1
        mv ./"$meta_id"_* versions.yml {outputs.results}/
    """,
)

env = flyte.TaskEnvironment.from_task("sparsesignatures", sparsesignatures_cmd.as_task())


@dataclass
class SparsesignaturesResult:
    mut_counts_rds: File  # <prefix>_mut_counts.rds
    cv_means_mse_rds: File  # <prefix>_cv_means_mse.rds
    best_params_config_rds: File  # <prefix>_best_params_config.rds
    nmf_lasso_out_rds: File  # <prefix>_nmf_Lasso_out.rds
    plot_rds: File  # <prefix>_plot_all.rds
    plot_pdf: File  # <prefix>_plot_all.pdf
    versions: File
    results: Dir


async def sparsesignatures(
    tsv_join: list[File], genome: str, prefix: str = "sample", args: str = ""
) -> SparsesignaturesResult:
    """Extract mutational signatures from the joint tables ``tsv_join`` (``genome``: GRCh37 or GRCh38)."""
    results = await sparsesignatures_cmd(tsv_join=tsv_join, genome=genome, prefix=prefix, args=args)

    async def pick(name: str) -> File:
        found = await results.get_file(name)
        if found is None:
            raise FileNotFoundError(f"sparsesignatures wrote no {name}")
        return found

    return SparsesignaturesResult(
        mut_counts_rds=await pick(f"{prefix}_mut_counts.rds"),
        cv_means_mse_rds=await pick(f"{prefix}_cv_means_mse.rds"),
        best_params_config_rds=await pick(f"{prefix}_best_params_config.rds"),
        nmf_lasso_out_rds=await pick(f"{prefix}_nmf_Lasso_out.rds"),
        plot_rds=await pick(f"{prefix}_plot_all.rds"),
        plot_pdf=await pick(f"{prefix}_plot_all.pdf"),
        versions=await pick("versions.yml"),
        results=results,
    )
