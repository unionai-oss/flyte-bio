"""tinc — TINC, tumour-in-normal contamination estimate from matched tumour and normal mutations.

Exposes :func:`tinc`, which joins ``tumour_sample``'s and ``normal_sample``'s
mutations (an RDS of per-sample mutation tables), fits TINC against the
copy-number segments of ``cna_rds`` (an RDS with ``segments``) and grades the
normal's contamination. It writes, with a ``versions.yml``:

- ``<prefix>_fit.rds``: TINC's fit;
- ``<prefix>_plot.rds`` / ``<prefix>_plot.pdf``: its plot;
- ``<prefix>_qc.csv``: ``sample,normal_contamination,normal_contamination_flag``,
  flagged 1 when the contamination level reaches ``normal_contamination_level``.

``args`` is upstream's ``ext.args`` map as Groovy prints it; the script needs
``normal_contamination_level``, and the default is upstream's test setting,
``[normal_contamination_level:3]``.

The R code is upstream's template, written out by an unquoted bash heredoc:
Groovy's ``\\$`` and ``\\\\`` escapes mean the same to bash, so only the
dotted placeholders are renamed (``$task.ext.args`` → ``${task_ext_args}``,
``$meta.tumour_sample`` → ``${meta_tumour_sample}``).
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

TINC_IMAGE = "quay.io/biocontainers/r-tinc:0.1.0--r44hdfd78af_0"
# Upstream's process_low asks for 2 CPUs; the script is single-threaded.
DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory="8Gi")

DEFAULT_ARGS = "[normal_contamination_level:3]"

# Upstream's templates/main_script.R, for an unquoted heredoc (see the module docstring).
TINC_TEMPLATE = r"""
#!/usr/bin/env Rscript

require(CNAqc)
require(tidyverse)
require(dplyr)
require(TINC)

parse_args = function(x) {
    x = gsub("\\\\[","",x)
    x = gsub("\\\\]","",x)
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

all_mutations = readRDS("$vcf_rds")
samples = names(all_mutations)

tumor_sample = "${meta_tumour_sample}"
normal_sample = "${meta_normal_sample}"

tumor_mutations = all_mutations[[tumor_sample]][['mutations']] %>%
    dplyr::select(chr, from, to, ref, alt, NV, DP, NR, VAF) %>%
    dplyr::filter(!is.na(DP)) %>%
    dplyr::rename(t_alt_count = NV, t_ref_count = NR, t_tot_count = DP, t_vaf = VAF)

normal_mutations = all_mutations[[normal_sample]][['mutations']] %>%
    dplyr::select(chr, from, to, ref, alt, NV, DP, NR, VAF) %>%
    dplyr::filter(!is.na(DP)) %>%
    dplyr::rename(n_alt_count = NV, n_ref_count = NR, n_tot_count = DP, n_vaf = VAF)

input_mut = dplyr::full_join(tumor_mutations, normal_mutations, by = c("chr", "from", "to", "ref", "alt")) %>%
    dplyr::mutate(t_vaf = case_when(is.na(t_vaf) ~ 1e-5, .default = t_vaf)) %>%
    dplyr::mutate(n_vaf = case_when(is.na(n_vaf) ~ 1e-5, .default = n_vaf)) %>%
    dplyr::mutate(t_vaf = as.numeric(t_vaf), n_vaf = as.numeric(n_vaf)) %>%
    dplyr::filter(!(is.na(t_alt_count))) %>%
    dplyr::filter(!(is.na(n_alt_count))) %>%
    dplyr::mutate(t_alt_count = as.numeric(t_alt_count), t_ref_count = as.numeric(t_ref_count), n_tot_count = as.numeric(n_tot_count), n_ref_count = as.numeric(n_ref_count)) %>%
    dplyr::filter((n_ref_count + n_alt_count) > 0) %>%
    dplyr::filter(t_vaf > 0)

CNAs = readRDS("$cna_rds")[['segments']]
TINC_fit = TINC::autofit(input = input_mut,
                cna = CNAs)

tinc_plot = plot(TINC_fit)

qc_res = TINC:::classification_normal(TINC_fit[['TIN']])

if (qc_res[["level"]] >= eval(parse(text=opt[["normal_contamination_level"]]))) {
    sample_contamination = dplyr::tibble(
        sample = tumor_sample,
        normal_contamination = qc_res[["level"]],
        normal_contamination_flag = 1
    )
} else {
    sample_contamination = dplyr::tibble(sample = tumor_sample,
    normal_contamination = qc_res[["level"]],
    normal_contamination_flag = 0
    )
}

write.table(file = paste0(opt[["prefix"]], "_qc.csv"), sep = ",", x = sample_contamination, col.names = T, row.names = F, quote = F)

saveRDS(file = paste0(opt[["prefix"]], "_plot.rds"), object = tinc_plot)
ggplot2::ggsave(plot = tinc_plot, filename = paste0(opt[["prefix"]], "_plot.pdf"), width = 210, height = 297, units="mm", dpi = 200)
saveRDS(file = paste0(opt[["prefix"]], "_fit.rds"), object = TINC_fit)

# version export
f <- file("versions.yml","w")
tidyverse_version <- sessionInfo()\$otherPkgs\$tidyverse\$Version
cnaqc_version <- sessionInfo()\$otherPkgs\$CNAqc\$Version
tinc_version <- sessionInfo()\$otherPkgs\$TINC\$Version
writeLines(paste0('"', "${task_process}", '"', ":"), f)
writeLines(paste("    CNAqc:", cnaqc_version), f)
writeLines(paste("    tidyverse:", tidyverse_version), f)
writeLines(paste("    TINC:", tinc_version), f)
close(f)
"""  # noqa: E501

tinc_cmd = shell.create(
    name="tinc",
    image=TINC_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "cna_rds": File,
        "vcf_rds": File,
        "tumour_sample": str,
        "normal_sample": str,
        "prefix": str,
        "args": str,
    },
    defaults={"args": DEFAULT_ARGS},
    outputs={"results": Dir},
    script=r"""
        CNA=({inputs.cna_rds}); SNV=({inputs.vcf_rds})
        meta_tumour_sample={inputs.tumour_sample}
        meta_normal_sample={inputs.normal_sample}
        meta_id={inputs.prefix}
        task_ext_prefix=$meta_id
        task_ext_args={inputs.args}
        task_ext_args=${task_ext_args:-null}
        task_process=TINC
        W=$(mktemp -d); cd "$W"
        ln -s "${CNA[0]}" "${SNV[0]}" .
        cna_rds=$(basename "${CNA[0]}"); vcf_rds=$(basename "${SNV[0]}")
        cat > main_script.R <<END_OF_TEMPLATE"""
    + TINC_TEMPLATE
    + r"""END_OF_TEMPLATE
        Rscript main_script.R || exit 1
        mv ./"$meta_id"_fit.rds ./"$meta_id"_plot.rds ./"$meta_id"_plot.pdf ./"$meta_id"_qc.csv versions.yml \
            {outputs.results}/ || exit 1
    """,
)

env = flyte.TaskEnvironment.from_task("tinc", tinc_cmd.as_task())


@dataclass
class TincResult:
    fit_rds: File  # <prefix>_fit.rds
    plot_rds: File  # <prefix>_plot.rds
    plot_pdf: File  # <prefix>_plot.pdf
    qc_csv: File  # <prefix>_qc.csv
    versions: File
    results: Dir


async def tinc(
    cna_rds: File,
    vcf_rds: File,
    tumour_sample: str,
    normal_sample: str,
    prefix: str = "sample",
    args: str = DEFAULT_ARGS,
) -> TincResult:
    """Estimate ``normal_sample``'s contamination by ``tumour_sample`` with TINC."""
    results = await tinc_cmd(
        cna_rds=cna_rds,
        vcf_rds=vcf_rds,
        tumour_sample=tumour_sample,
        normal_sample=normal_sample,
        prefix=prefix,
        args=args,
    )

    async def pick(name: str) -> File:
        found = await results.get_file(name)
        if found is None:
            raise FileNotFoundError(f"tinc wrote no {name}")
        return found

    return TincResult(
        fit_rds=await pick(f"{prefix}_fit.rds"),
        plot_rds=await pick(f"{prefix}_plot.rds"),
        plot_pdf=await pick(f"{prefix}_plot.pdf"),
        qc_csv=await pick(f"{prefix}_qc.csv"),
        versions=await pick("versions.yml"),
        results=results,
    )
