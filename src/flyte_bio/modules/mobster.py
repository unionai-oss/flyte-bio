"""mobster — MOBSTER, subclonal deconvolution of one tumour sample's VAFs with a neutral tail.

Exposes :func:`mobster`, which fits MOBSTER's mixture of Beta clusters and a
Pareto tail to the diploid (``1:1``), non-blacklisted mutations of each tumour
sample in ``tumour_sample`` (comma-separated). The input is either an
``.rds`` multi-sample CNAqc object (``m_cnaqc``) or a CSV joint table with
``sample_id``, ``karyotype``, ``VAF``, ``purity`` … columns. It writes, with a
``versions.yml``:

- ``<prefix>_mobster_fit.rds``: all fits; ``<prefix>_mobster_best_fit.rds``:
  the best one, with its evolutionary parameters when it has a tail;
- ``<prefix>_mobster_best_plots.rds``: the best fit's plot;
- ``<prefix>_mobster_report.{rds,pdf,png}``: the model-selection report.

``args`` is upstream's ``ext.args`` map as Groovy prints it,
``[K:2:3, tail:TRUE, pi_cutoff:0.02, n_cutoff:10, min_VAF:0.05, qc_chr:FALSE]``
(the default, upstream's test settings); the script needs all six. The fit
takes MOBSTER's random initial conditions: the template sets no seed.

The R code is upstream's template, written out by an unquoted bash heredoc:
Groovy's ``\\$`` and ``\\\\`` escapes mean the same to bash, so only the
dotted placeholders are renamed (``$task.ext.args`` → ``${task_ext_args}``,
``$meta.patient`` → ``${meta_patient}``).
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

MOBSTER_IMAGE = "community.wave.seqera.io/library/r-cnaqc_r-mobster_r-cli_r-dplyr_r-ggplot2:cc406f90216d5a33"
# The fits run serially (parallel = FALSE).
DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory="8Gi")

DEFAULT_ARGS = "[K:2:3, tail:TRUE, pi_cutoff:0.02, n_cutoff:10, min_VAF:0.05, qc_chr:FALSE]"

# Upstream's templates/main_script.R, for an unquoted heredoc (see the module docstring).
MOBSTER_TEMPLATE = r"""
#!/usr/bin/env Rscript

parse_args = function(x) {
    x = gsub("\\\\[","",x)
    x = gsub("\\\\]","",x)
    # giving errors when we have lists like c(xxx, xxx) since it will separate it
    # args_list = unlist(strsplit(x, ', ')[[1]])
    args_list = unlist(strsplit(x, ", (?=[^)]*(?:\\\\(|\$))", perl=TRUE))
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


# Script #####

library(CNAqc)
library(mobster)
library(dplyr)
library(cli)
library(ggplot2)

description = "${meta_patient}"
samples = strsplit(x="${meta_tumour_sample}", ",") %>% unlist()  # list of samples

## read mCNAqc object
if ( grepl(".rds\$", tolower("$rds_join")) ) {
  obj = readRDS("$rds_join")
  if (as.logical(opt[["qc_chr"]])){
    chroms <- intersect(names(obj), paste0("chr", 1:22))
    input_table <- lapply(chroms,function(c){
      original = obj[[c]] %>% get_sample(sample=samples, which_obj="original")
      joint_table_chr = lapply(names(original),
                               function(sample_name) {
                                 purity = original[[sample_name]][["purity"]]
                                 table_s = original[[sample_name]] %>%
                                   # CNAqc::subset_by_segment_karyotype("1:1") %>%
                                   CNAqc::Mutations() %>%
                                   dplyr::filter(karyotype=="1:1") %>%
                                   dplyr::mutate(sample_id=sample_name, purity=purity) %>%
                                   dplyr::filter(blacklisted==FALSE)
                                 if (nrow(table_s) == 0) {
                                   cli::cli_alert_warning("Sample {sample_name} has no diploid mutations!")
                                 }
                                 return(table_s)
                               }) %>% dplyr::bind_rows()
    })  %>% dplyr::bind_rows()
  } else {
    if (class(obj) == "m_cnaqc") {
      original = obj %>% get_sample(sample=samples, which_obj="original")
      input_table = lapply(names(original),
                           function(sample_name) {
                             purity = original[[sample_name]][["purity"]]
                             table_s = original[[sample_name]] %>%
                               # keep only mutations on the diploid karyotype
                               CNAqc::subset_by_segment_karyotype("1:1") %>%
                               CNAqc::Mutations() %>%
                               dplyr::mutate(sample_id=sample_name, purity=purity) %>%
                               dplyr::filter(blacklisted==FALSE)
                             if (nrow(table_s) == 0) {
                               cli::cli_alert_warning("Sample {sample_name} has no diploid mutations!")
                             }
                             return(table_s)
                           }) %>% dplyr::bind_rows()
    } else {
      cli::cli_abort("Object of class {class($rds_join)} not supported.")
    }
  }
} else {
  input_table = read.csv("$rds_join")
  input_table = input_table %>%
    dplyr::filter(karyotype=="1:1") %>%
    dplyr::filter(chr %in% c(paste0('chr', 1:22), as.character(1:22))) %>%
    dplyr::filter(blacklisted==FALSE)
}

## Function to run a single MOBSTER fit
run_mobster_fit = function(joint_table, descr) {
    purity = unique(joint_table[["purity"]]) %>% as.numeric()

    if(length(purity) > 1)
        cli::cli_abort("Multiple values of purity are found in one sample (sample ID: {descr}).")

    # remove mutations with VAF < min_VAF
    inp_tb = joint_table %>%
        dplyr::filter(VAF < 1) %>%
        dplyr::filter(VAF>=as.numeric(opt[["min_VAF"]])) %>%
        dplyr::filter(karyotype=="1:1")

    mobster_fit(x = inp_tb,
                K = eval(parse(text=opt[["K"]])),
                tail = eval(parse(text=opt[["tail"]])),
                pi_cutoff = as.numeric(opt[["pi_cutoff"]]),
                N_cutoff = as.integer(opt[["n_cutoff"]]),
                parallel = FALSE,
                description = descr)
}

## Run MOBSTER on each sample
lapply(samples, function(sample_name) {

    fit = run_mobster_fit(joint_table=input_table %>% dplyr::filter(sample_id == !!sample_name),
                          descr=paste(description, sample_name, sep=":"))

    if (fit[["best"]][["fit.tail"]]) {
        evoparams = evolutionary_parameters(fit)
        fit[["evolutionary_parameters"]] = evoparams
        fit[["best"]][["evolutionary_parameters"]] = evoparams
    }

    best_fit = fit[["best"]]
    plot_fit = plot(best_fit)

    saveRDS(object=fit, file=paste0(opt[["prefix"]], "_mobster_fit.rds"))
    saveRDS(object=best_fit, file=paste0(opt[["prefix"]], "_mobster_best_fit.rds"))
    saveRDS(object=plot_fit, file=paste0(opt[["prefix"]], "_mobster_best_plots.rds"))

    # save report plots
    report_fig = mobster::plot_model_selection(fit)
    saveRDS(report_fig, file=paste0(opt[["prefix"]], "_mobster_report.rds"))
    ggplot2::ggsave(plot=report_fig, filename=paste0(opt[["prefix"]], "_mobster_report.pdf"), height=210, width=210, units="mm", dpi = 200)
    ggplot2::ggsave(plot=report_fig, filename=paste0(opt[["prefix"]], "_mobster_report.png"), height=210, width=210, units="mm", dpi = 200)
})


# Version export #####
f = file("versions.yml","w")
cnaqc_version = sessionInfo()\$otherPkgs\$CNAqc\$Version
mobster_version = sessionInfo()\$otherPkgs\$mobster\$Version
cli_version = sessionInfo()\$otherPkgs\$cli\$Version
dplyr_version = sessionInfo()\$otherPkgs\$dplyr\$Version
ggplot2_version = sessionInfo()\$otherPkgs\$ggplot2\$Version
writeLines(paste0('"', "${task_process}", '"', ":"), f)
writeLines(paste("    CNAqc:", cnaqc_version), f)
writeLines(paste("    mobster:", mobster_version), f)
writeLines(paste("    cli:", cli_version), f)
writeLines(paste("    dplyr:", dplyr_version), f)
writeLines(paste("    ggplot2:", ggplot2_version), f)
close(f)
"""  # noqa: E501

mobster_cmd = shell.create(
    name="mobster",
    image=MOBSTER_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"rds_join": File, "patient": str, "tumour_sample": str, "prefix": str, "args": str},
    defaults={"args": DEFAULT_ARGS},
    outputs={"results": Dir},
    script=r"""
        IN=({inputs.rds_join})
        meta_patient={inputs.patient}
        meta_tumour_sample={inputs.tumour_sample}
        meta_id={inputs.prefix}
        task_ext_prefix=$meta_id
        task_ext_args={inputs.args}
        task_ext_args=${task_ext_args:-null}
        task_process=MOBSTER
        W=$(mktemp -d); cd "$W"
        ln -s "${IN[0]}" .
        rds_join=$(basename "${IN[0]}")
        cat > main_script.R <<END_OF_TEMPLATE"""
    + MOBSTER_TEMPLATE
    + r"""END_OF_TEMPLATE
        Rscript main_script.R || exit 1
        mv ./"$meta_id"_mobster_* versions.yml {outputs.results}/
    """,
)

env = flyte.TaskEnvironment.from_task("mobster", mobster_cmd.as_task())


@dataclass
class MobsterResult:
    fit_rds: File  # <prefix>_mobster_fit.rds
    best_fit_rds: File  # <prefix>_mobster_best_fit.rds
    best_plots_rds: File  # <prefix>_mobster_best_plots.rds
    report_rds: File  # <prefix>_mobster_report.rds
    report_pdf: File
    report_png: File
    versions: File
    results: Dir


async def mobster(
    rds_join: File, patient: str, tumour_sample: str, prefix: str = "sample", args: str = DEFAULT_ARGS
) -> MobsterResult:
    """Fit MOBSTER to ``tumour_sample``'s diploid mutations in ``rds_join`` (an m_cnaqc RDS or a CSV)."""
    results = await mobster_cmd(
        rds_join=rds_join, patient=patient, tumour_sample=tumour_sample, prefix=prefix, args=args
    )

    async def pick(name: str) -> File:
        found = await results.get_file(f"{prefix}_mobster_{name}")
        if found is None:
            raise FileNotFoundError(f"mobster wrote no {prefix}_mobster_{name}")
        return found

    versions = await results.get_file("versions.yml")
    if versions is None:
        raise FileNotFoundError("mobster wrote no versions.yml")
    return MobsterResult(
        fit_rds=await pick("fit.rds"),
        best_fit_rds=await pick("best_fit.rds"),
        best_plots_rds=await pick("best_plots.rds"),
        report_rds=await pick("report.rds"),
        report_pdf=await pick("report.pdf"),
        report_png=await pick("report.png"),
        versions=versions,
        results=results,
    )
