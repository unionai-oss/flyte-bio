"""viber — VIBER, variational Binomial mixture clustering of read counts across tumour samples.

Exposes :func:`viber`, which clusters the diploid (``1:1``), non-blacklisted
mutations shared by ``tumour_samples`` on their read counts (``NV`` of
``DP``). The input is either an ``.rds`` multi-sample CNAqc object
(``m_cnaqc``) or a CSV joint table. It fits VIBER, then applies VIBER's
cluster-filtering heuristic (keeping the plain fit when the heuristic removes
every cluster), and writes, with a ``versions.yml``:

- ``<prefix>_viber_best_fit.rds`` and ``<prefix>_viber_best_heuristic_fit.rds``;
- ``<prefix>_viber_best_fit_plots.rds`` and ``<prefix>_viber_best_heuristic_fit_plots.rds``;
- ``<prefix>_viber_report.{rds,pdf,png}``.

``args`` is upstream's ``ext.args`` map as Groovy prints it; the script needs
``K``, ``binomial_cutoff``, ``dimensions_cutoff``, ``pi_cutoff``,
``re_assign`` and ``qc_chr``, and the default is upstream's test settings.
VIBER's initialisation is unseeded, as upstream.

The R code is upstream's template, written out by an unquoted bash heredoc:
Groovy's ``\\$`` and ``\\\\`` escapes mean the same to bash, so only the
dotted placeholders are renamed (``$task.ext.args`` → ``${task_ext_args}``,
``$meta.patient`` → ``${meta_patient}``). ``tumour_samples`` renders as
Groovy prints a list, ``[A, B]``.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

VIBER_IMAGE = "community.wave.seqera.io/library/r-cnaqc_r-viber_r-cli_r-dplyr_pruned:f00a3b14d68a7bac"
# Upstream's process_single; the template turns easypar's parallelism off.
DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory="8Gi")

DEFAULT_ARGS = "[K:10, binomial_cutoff:0.05, dimensions_cutoff:1, pi_cutoff:0.02, re_assign:FALSE, qc_chr:FALSE]"

# Upstream's templates/viber_main_script.R, for an unquoted heredoc (see the module docstring).
VIBER_TEMPLATE = r"""
#!/usr/bin/env Rscript

options(easypar.parallel=FALSE)

parse_args = function(x) {
    x = gsub("\\\\[","",x)
    x = gsub("\\\\]","",x)
    # giving errors when we have lists like c(xxx, xxx) since it will separate it
    # args_list = unlist(strsplit(x, ', ')[[1]])
    args_list = unlist(strsplit(x, ", (?=[^)]*(?:\\\\(|\$))", perl=TRUE))
    # args_vals = lapply(args_list, function(x) strsplit(x, split=":")[[1]])
    args_vals = lapply(args_list, function(x) {
        x_splt = strsplit(x, split=":")[[1]]
        c(x_splt[1],  paste(x_splt[2:length(x_splt)], collapse=":"))
    })

    # Ensure the option vectors are length 2 (key/ value) to catch empty ones
    args_vals = lapply(args_vals, function(z){ length(z) = 2; z})

    parsed_args = structure(lapply(args_vals, function(x) x[2]), names = lapply(args_vals, function(x) x[1]))
    parsed_args[! is.na(parsed_args)]
}

opt = list(
    prefix = ifelse('${task_ext_prefix}' == 'null', '${meta_id}', '${task_ext_prefix}')
)
args_opt = parse_args('${task_ext_args}')
for ( ao in names(args_opt)) opt[[ao]] = args_opt[[ao]]
print(opt)

# Script #####

library(VIBER)
library(dplyr)
library(tidyr)
library(ggplot2)
library(CNAqc)
library(ggpubr)

samples = substr("$tumour_samples", 2, nchar("$tumour_samples")-1)
samples = strsplit(samples, ", ")[[1]]

cli::cli_text("Patient ${meta_patient} with samples $tumour_samples.")
cli::cli_text("Input file: $rds_join.")

if ( grepl(".rds\$", tolower("$rds_join")) ) {
  input_obj = readRDS("$rds_join")
  if (as.logical(opt[["qc_chr"]])){
    chroms <- intersect(names(input_obj), paste0("chr", 1:22))
    joint_table <- lapply(chroms,function(c){
      shared = input_obj[[c]] %>% get_sample(sample=samples, which_obj="shared")
      joint_table_chr = lapply(names(shared),
                               function(sample_name) {
                                 table_s = shared[[sample_name]] %>%
                                   # CNAqc::subset_by_segment_karyotype("1:1") %>%
                                   CNAqc::Mutations() %>%
                                   dplyr::filter(karyotype=="1:1") %>%
                                   dplyr::mutate(sample_id=sample_name) %>%
                                   dplyr::filter(blacklisted==FALSE)
                                 if (nrow(table_s) == 0) {
                                   cli::cli_alert_warning("Sample {sample_name} has no diploid mutations!")
                                 }
                                 return(table_s)
                               }) %>% dplyr::bind_rows()
    })  %>% dplyr::bind_rows()
  } else {
    if (class(input_obj) == "m_cnaqc") {
      shared = input_obj %>% get_sample(sample=samples, which_obj="shared")
      joint_table = lapply(names(shared),
                           function(sample_name) {
                             table_s = shared[[sample_name]] %>%
                               CNAqc::subset_by_segment_karyotype("1:1") %>%
                               CNAqc::Mutations() %>%
                               dplyr::mutate(sample_id=sample_name) %>%
                               dplyr::filter(chr %in% c(paste0('chr', 1:22), as.character(1:22))) %>%
                               dplyr::filter(blacklisted==FALSE)

                             if (nrow(table_s) == 0) {
                               cli::cli_alert_warning("Sample {sample_name} has no diploid mutations!")
                             }
                             return(table_s)
                           }) %>% dplyr::bind_rows()
    } else {
      cli::cli_alert_warning("Object of class {class(input_obj)} not supported.")
      return()
    }
  }
} else {
  joint_table = read.csv("$rds_join")
  joint_table = joint_table %>%
    dplyr::filter(karyotype=="1:1") %>%
    dplyr::filter(chr %in% c(paste0('chr', 1:22), as.character(1:22))) %>%
    dplyr::filter(blacklisted==FALSE)
}

if (nrow(joint_table) == 0) {
  cli::cli_alert_warning("No samples contain diploid mutations!")
}

cli::cli_text("Joint table created.")

## Read input joint table
input_tab = joint_table %>%
  dplyr::mutate(VAF=replace(VAF, VAF==0, 1e-7))

## Convert the input table into longer format
reads_data = input_tab %>%
  dplyr::select(chr, from, ref, alt, NV, DP, VAF, sample_id,driver_label,is_driver) %>%
  dplyr::rename(gene=driver_label) %>%
  dplyr::rename(driver=is_driver) %>%
  tidyr::pivot_wider(names_from="sample_id",
                     values_from=c("NV","DP","VAF"), names_sep=".",values_fill=0)

## Extract DP (depth)
dp = reads_data %>%
  dplyr::select(dplyr::starts_with("DP")) %>%
  dplyr::rename_all(function(x) stringr::str_remove_all(x,"DP."))

## Extract NV (alt_counts)
nv = reads_data %>%
  dplyr::select(dplyr::starts_with("NV")) %>%
  dplyr::rename_all(function(x) stringr::str_remove_all(x,"NV."))

# Standard fit
viber_K = as.integer(opt[["K"]])

cli::cli_text("Starting standard fit")
st_fit = VIBER::variational_fit(nv, dp,
                                K=viber_K,
                                data=reads_data)
st_fit[["description"]]="${meta_patient}"
cli::cli_text("End standard fit")
best_fit = best_fit_heuristic = st_fit

# If all clusters are removed -> keep the origianl best fit
tryCatch(expr = {
    # Apply the heuristic
    best_fit_heuristic = VIBER::choose_clusters(st_fit,
        binomial_cutoff=as.numeric(opt[["binomial_cutoff"]]),
        dimensions_cutoff=as.integer(opt[["dimensions_cutoff"]]),
        pi_cutoff=as.numeric(opt[["pi_cutoff"]]),
        re_assign=as.logical(opt[["re_assign"]])
        )
    }, error = function(e) {
        print(e)
        best_fit_heuristic <<- st_fit
    } )

# Save fits
saveRDS(best_fit, file=paste0(opt[["prefix"]], "_viber_best_fit.rds"))
saveRDS(best_fit_heuristic, file = paste0(opt[["prefix"]], "_viber_best_heuristic_fit.rds"))

# Save plots
n_samples = ncol(best_fit[["x"]]) - 1
if (n_samples >1) {
  print("multisample mode on")
  plot_fit = plot(best_fit)
  plot_fit_heuristic = plot(best_fit_heuristic)
} else {
  plot_fit = plot_mixing_proportions(best_fit)
  plot_fit_heuristic = plot_mixing_proportions(best_fit_heuristic)
}

saveRDS(plot_fit, file=paste0(opt[["prefix"]], "_viber_best_fit_plots.rds"))
saveRDS(plot_fit_heuristic, file=paste0(opt[["prefix"]], "_viber_best_heuristic_fit_plots.rds"))

# Save report plot
n_samples = ncol(best_fit[["x"]]) - 1
marginals = multivariate = ggplot()

try(expr = {marginals <<- VIBER::plot_1D(best_fit)} )
try(expr = {multivariate = plot(best_fit)})
try(expr = {multivariate = ggpubr::ggarrange(plotlist = multivariate)})
top_p = ggpubr::ggarrange(plotlist = list(marginals, multivariate), widths=ifelse(n_samples>2, c(1,2), c(2,1)))

mix_p = VIBER::plot_mixing_proportions(best_fit)
binom = VIBER::plot_peaks(best_fit)
elbo = VIBER::plot_ELBO(best_fit)
bottom_p = ggpubr::ggarrange(plotlist = list(mix_p, binom, elbo), widths = c(1,2,1), nrow = 1)

report_fig = ggpubr::ggarrange(top_p, bottom_p, nrow = 2, heights = ifelse(n_samples>2, c(3,1), c(2,1)))
saveRDS(report_fig, file=paste0(opt[["prefix"]], "_viber_report.rds"))
ggplot2::ggsave(plot=report_fig, filename=paste0(opt[["prefix"]], "_viber_report.pdf"), height=210, width=210, units="mm", dpi = 200)
ggplot2::ggsave(plot=report_fig, filename=paste0(opt[["prefix"]], "_viber_report.png"), height=210, width=210, units="mm", dpi = 200)

# version export
f = file("versions.yml","w")
cnaqc_version = sessionInfo()\$otherPkgs\$CNAqc\$Version
viber_version = sessionInfo()\$otherPkgs\$VIBER\$Version
cli_version = sessionInfo()\$otherPkgs\$cli\$Version
dplyr_version = sessionInfo()\$otherPkgs\$dplyr\$Version
tidyr_version = sessionInfo()\$otherPkgs\$tidyr\$Version
ggplot2_version = sessionInfo()\$otherPkgs\$ggplot2\$Version
ggpubr_version = sessionInfo()\$otherPkgs\$ggpubr\$Version
writeLines(paste0('"', "${task_process}", '"', ":"), f)
writeLines(paste("    CNAqc:", cnaqc_version), f)
writeLines(paste("    VIBER:", viber_version), f)
writeLines(paste("    cli:", cli_version), f)
writeLines(paste("    dplyr:", dplyr_version), f)
writeLines(paste("    tidyr:", tidyr_version), f)
writeLines(paste("    ggplot2:", ggplot2_version), f)
writeLines(paste("    ggpubr:", ggpubr_version), f)
close(f)
"""  # noqa: E501

viber_cmd = shell.create(
    name="viber",
    image=VIBER_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"rds_join": File, "patient": str, "tumour_samples": str, "prefix": str, "args": str},
    defaults={"args": DEFAULT_ARGS},
    outputs={"results": Dir},
    script=r"""
        IN=({inputs.rds_join})
        meta_patient={inputs.patient}
        tumour_samples={inputs.tumour_samples}
        meta_id={inputs.prefix}
        task_ext_prefix=$meta_id
        task_ext_args={inputs.args}
        task_ext_args=${task_ext_args:-null}
        task_process=VIBER
        W=$(mktemp -d); cd "$W"
        ln -s "${IN[0]}" .
        rds_join=$(basename "${IN[0]}")
        cat > viber_main_script.R <<END_OF_TEMPLATE"""
    + VIBER_TEMPLATE
    + r"""END_OF_TEMPLATE
        Rscript viber_main_script.R || exit 1
        mv ./"$meta_id"_viber_* versions.yml {outputs.results}/
    """,
)

env = flyte.TaskEnvironment.from_task("viber", viber_cmd.as_task())


@dataclass
class ViberResult:
    best_fit_rds: File  # <prefix>_viber_best_fit.rds
    best_heuristic_fit_rds: File  # <prefix>_viber_best_heuristic_fit.rds
    best_fit_plots_rds: File  # <prefix>_viber_best_fit_plots.rds
    best_heuristic_fit_plots_rds: File  # <prefix>_viber_best_heuristic_fit_plots.rds
    report_rds: File  # <prefix>_viber_report.rds
    report_pdf: File
    report_png: File
    versions: File
    results: Dir


async def viber(
    rds_join: File, patient: str, tumour_samples: list[str], prefix: str = "sample", args: str = DEFAULT_ARGS
) -> ViberResult:
    """Cluster ``tumour_samples``' diploid mutations in ``rds_join`` (an m_cnaqc RDS or a CSV) with VIBER."""
    results = await viber_cmd(
        rds_join=rds_join,
        patient=patient,
        tumour_samples="[" + ", ".join(tumour_samples) + "]",
        prefix=prefix,
        args=args,
    )

    async def pick(name: str) -> File:
        found = await results.get_file(f"{prefix}_viber_{name}")
        if found is None:
            raise FileNotFoundError(f"viber wrote no {prefix}_viber_{name}")
        return found

    versions = await results.get_file("versions.yml")
    if versions is None:
        raise FileNotFoundError("viber wrote no versions.yml")
    return ViberResult(
        best_fit_rds=await pick("best_fit.rds"),
        best_heuristic_fit_rds=await pick("best_heuristic_fit.rds"),
        best_fit_plots_rds=await pick("best_fit_plots.rds"),
        best_heuristic_fit_plots_rds=await pick("best_heuristic_fit_plots.rds"),
        report_rds=await pick("report.rds"),
        report_pdf=await pick("report.pdf"),
        report_png=await pick("report.png"),
        versions=versions,
        results=results,
    )
