"""pyclonevi — PyClone-VI, clustering of mutations' cellular prevalence across tumour samples.

Exposes :func:`pyclonevi`, which turns a tab-separated joint mutation table
(``chr``, ``from``, ``alt``, ``karyotype``, ``Indiv``, ``NR``, ``NV``,
``purity``, ``blacklisted`` …) into PyClone-VI's input, keeps the rows of
``tumour_samples``, fits PyClone-VI (``pyclone-vi fit``, then
``write-results-file``) and writes a cluster table for
:func:`flyte_bio.modules.ctree.ctree`. Outputs, with a ``versions.yml``:

- ``<prefix>_pyclone_input_all_samples.tsv`` and ``<prefix>_pyclone_input.tsv``:
  the PyClone-VI input, before and after keeping ``tumour_samples``;
- ``<prefix>_all_fits.h5`` and ``<prefix>_best_fit.txt``: the fits and the best one;
- ``<prefix>_cluster_table.csv``: the ctree input (tab-separated, despite its name).

``args`` is upstream's ``ext.args`` map as Groovy prints it; ``n_cluster``
(PyClone-VI's ``-c``) is required, and the default is upstream's test setting,
``[n_cluster:20]``. PyClone-VI's random restarts are unseeded, as upstream.

The Python code is upstream's template, written out by an unquoted bash
heredoc: Groovy's ``\\$`` and ``\\\\`` escapes mean the same to bash, so only
the dotted placeholders are renamed (``$task.ext.args`` →
``${task_ext_args}``). ``tumour_samples`` renders as Groovy prints a list,
``[A, B]``.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

PYCLONEVI_IMAGE = "quay.io/biocontainers/pyclone-vi:0.2.0--pyhdfd78af_0"
# pyclone-vi fit runs one thread by default (-t 1), as upstream calls it.
DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory="8Gi")

DEFAULT_ARGS = "[n_cluster:20]"

# Upstream's templates/main_script.py, for an unquoted heredoc (see the module docstring).
PYCLONEVI_TEMPLATE = r"""
#!/usr/bin/env python

import subprocess

import numpy as np
import pandas as pd


def parse_args(x):
    x = x.replace("[", "")
    x = x.replace("]", "")

    args_list = x.split(", ")
    parsed_args = {i.split(":")[0]: i.split(":")[1] for i in args_list if i.split(":")[1] is not None}
    return parsed_args


opt = dict()
opt["prefix"] = "${task_ext_prefix}" if "${task_ext_prefix}" == None else "${meta_id}"
args_opt = parse_args("${task_ext_args}")
for ao_k, ao_v in args_opt.items():
    opt[ao_k] = ao_v

print("$tumour_samples")

# Script #####


# input data preprocessing
def create_pyclone_input(input_data, patient_id, output_data):
    df = pd.read_csv(input_data, sep="\t", header=0)
    df = df[~df["chr"].isin(["chrX", "chrY", "X", "Y"])]

    df = df[~df["blacklisted"]]

    df["normal_cn"] = 2
    df["patient_id"] = patient_id
    df["mutation_id"] = df.apply(lambda row: f"{patient_id}:{row['chr'][3:]}:{row['from']}:{row['alt']}", axis=1)
    df[["major_cn", "minor_cn"]] = df["karyotype"].str.split(":", expand=True)
    df = df.drop(columns=["karyotype"])
    df = df.rename(columns={"Indiv": "sample_id", "NR": "ref_counts", "NV": "alt_counts", "purity": "tumour_content"})

    column_names = [
        "mutation_id",
        "patient_id",
        "sample_id",
        "ref_counts",
        "alt_counts",
        "normal_cn",
        "major_cn",
        "minor_cn",
        "tumour_content",
        "driver_label",
        "is_driver",
    ]
    column_names_red = list(set(column_names) & set(df.columns))

    df = df.loc[:, column_names_red]
    df.to_csv(output_data, sep="\t", index=False, header=df.columns)


def pyclone_ctree(joint, best_fit, ctree_input):
    ## Read pyclone best fit table
    best_fit_file = best_fit
    df_output = pd.read_csv(best_fit_file, sep="\t")

    ## Read pyclone input table
    joint_table = joint
    df_input = pd.read_csv(joint_table, sep="\t")
    ## Caluclate number of mutations per cluster and add to the subset orginal df_output dataframe

    df_output_small = df_output[["mutation_id", "sample_id", "cluster_id", "cellular_prevalence"]]
    df_output_small["cluster_id"] = "C" + df_output_small["cluster_id"].astype(str)
    df_output_small["nMuts"] = df_output_small.groupby("cluster_id")["mutation_id"].transform("nunique")

    ## Create cluster table
    df_clusters = (
        df_output_small[["sample_id", "cellular_prevalence", "cluster_id"]]
        .drop_duplicates()
        .pivot(index="sample_id", columns="cluster_id", values="cellular_prevalence")
    )
    ## Find the clonal cluster and add the colum 'is.clonal' to identify it

    samples = df_output_small["sample_id"].unique()
    top_clusters = []
    for s in samples:
        max_values = df_clusters.loc[s].max()
        indexes = np.where(df_clusters.loc[s] == max_values)[0]
        top_clusters.append(list(df_clusters.columns[indexes]))
    top_clusters_all = [item for sublist in top_clusters for item in sublist]

    clonal_cluster = max(set(top_clusters_all), key=top_clusters_all.count)

    i = df_output_small[df_output_small["cluster_id"] == clonal_cluster].index

    df_output_small["is.clonal"] = "F"
    df_output_small.loc[i, "is.clonal"] = "T"

    ## Merge the input information about driver genes to the clusters
    if "driver_label" not in df_input.columns or "is_driver" not in df_input.columns:
        df_input["driver_label"] = "NA"
        df_input["is_driver"] = False
        df_input.loc[1, "driver_label"] = ""
        df_input.loc[1, "is_driver"] = True

    df_merged = pd.merge(
        df_output_small,
        df_input[["mutation_id", "sample_id", "patient_id", "driver_label", "is_driver"]],
        on=["mutation_id", "sample_id"],
        how="inner",
    )
    df_final = df_merged.drop_duplicates(subset=["sample_id", "driver_label", "cluster_id", "is.clonal", "is_driver"])
    df_final = df_final.rename(
        columns={
            "cellular_prevalence": "CCF",
            "cluster_id": "cluster",
            "driver_label": "variantID",
            "patient_id": "patientID",
            "is_driver": "is.driver",
        }
    )
    df_final["is.driver"] = df_final["is.driver"].replace({False: "F", True: "T"})

    ## Final ctree input dataframe
    final_jt_file = ctree_input
    df_final.to_csv(final_jt_file, sep="\t", index=False, header=True)


if __name__ == "__main__":
    create_pyclone_input(
        input_data="$rds_join", patient_id="${meta_patient}", output_data=opt["prefix"] + "_pyclone_input_all_samples.tsv"
    )

    column_name = "sample_id"
    input_all_samples = pd.read_csv(opt["prefix"] + "_pyclone_input_all_samples.tsv", sep="\t")

    samples_list = "$tumour_samples".replace("[", "")
    samples_list = samples_list.replace("]", "")
    samples_list = samples_list.split(sep=", ")
    input_filtered = input_all_samples[input_all_samples[column_name].isin(samples_list)]
    input_filtered.to_csv(opt["prefix"] + "_pyclone_input.tsv", sep="\t")

    # Run pyclone
    pyclonevi_run = (
        "pyclone-vi fit -i "
        + opt["prefix"]
        + "_pyclone_input.tsv -o "
        + opt["prefix"]
        + "_all_fits.h5 -c "
        + opt["n_cluster"]
    )
    subprocess.run(pyclonevi_run, shell=True)

    pyclonevi_write = (
        "pyclone-vi write-results-file -i " + opt["prefix"] + "_all_fits.h5 -o " + opt["prefix"] + "_best_fit.txt"
    )
    subprocess.run(pyclonevi_write, shell=True)

    # # Pyclone ctree
    pyclone_ctree(
        joint=opt["prefix"] + "_pyclone_input.tsv",
        best_fit=opt["prefix"] + "_best_fit.txt",
        ctree_input=opt["prefix"] + "_cluster_table.csv",
    )

    # Version
    version = (
        subprocess.check_output("pyclone-vi --version | sed 's/.*version //g'", shell=True).decode().split("\\n")[0]
    )

    f = open("versions.yml", "a")
    f.write("PYCLONEVI:\\n")
    f.write(f"    pyclonevi: {version}\\n")
    f.close()
"""  # noqa: E501

pyclonevi_cmd = shell.create(
    name="pyclonevi",
    image=PYCLONEVI_IMAGE,
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
        W=$(mktemp -d); cd "$W"
        ln -s "${IN[0]}" .
        rds_join=$(basename "${IN[0]}")
        cat > main_script.py <<END_OF_TEMPLATE"""
    + PYCLONEVI_TEMPLATE
    + r"""END_OF_TEMPLATE
        python main_script.py || exit 1
        mv ./"$meta_id"_*.tsv ./"$meta_id"_cluster_table.csv ./"$meta_id"_all_fits.h5 ./"$meta_id"_best_fit.txt \
            versions.yml {outputs.results}/ || exit 1
    """,
)

env = flyte.TaskEnvironment.from_task("pyclonevi", pyclonevi_cmd.as_task())


@dataclass
class PycloneviResult:
    cluster_table: File  # <prefix>_cluster_table.csv, the ctree input
    pyclone_input: list[File]  # <prefix>_pyclone_input.tsv, <prefix>_pyclone_input_all_samples.tsv
    all_fits: File  # <prefix>_all_fits.h5
    best_fit: File  # <prefix>_best_fit.txt
    versions: File
    results: Dir


async def pyclonevi(
    rds_join: File, patient: str, tumour_samples: list[str], prefix: str = "sample", args: str = DEFAULT_ARGS
) -> PycloneviResult:
    """Cluster ``tumour_samples``' mutations from the joint table ``rds_join`` with PyClone-VI."""
    results = await pyclonevi_cmd(
        rds_join=rds_join,
        patient=patient,
        tumour_samples="[" + ", ".join(tumour_samples) + "]",
        prefix=prefix,
        args=args,
    )

    async def pick(name: str) -> File:
        found = await results.get_file(name)
        if found is None:
            raise FileNotFoundError(f"pyclonevi wrote no {name}")
        return found

    return PycloneviResult(
        cluster_table=await pick(f"{prefix}_cluster_table.csv"),
        pyclone_input=[
            await pick(f"{prefix}_pyclone_input.tsv"),
            await pick(f"{prefix}_pyclone_input_all_samples.tsv"),
        ],
        all_fits=await pick(f"{prefix}_all_fits.h5"),
        best_fit=await pick(f"{prefix}_best_fit.txt"),
        versions=await pick("versions.yml"),
        results=results,
    )
