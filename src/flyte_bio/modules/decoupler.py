"""decoupler — pathway and transcription-factor activity inference.

Exposes :func:`decoupler`, which runs decoupler-py's ``decouple`` over a
features x samples matrix (or, with ``--transpose TRUE`` in ``args``, one
column of a differential results table, by default ``log2FoldChange``) and a
prior-knowledge network (``source``, ``target``, ``weight``), with the
methods ``--methods`` lists (default ``ulm``) plus their consensus. With
``--ensembl_ids TRUE``, the matrix's Ensembl IDs are first mapped to symbols
through the annotation table's ``--features_id_col`` and
``--features_symbol_col`` columns. For each method and the consensus it
writes upstream's outputs:

- ``<prefix>_<method>_estimate_decoupler.tsv`` and
  ``<prefix>_<method>_pvals_decoupler.tsv``;
- ``<prefix>_<method>_estimate_decoupler_plot.png``;

and ``versions.yml``. The Python is upstream's template, unchanged but for
taking the templated values (matrix, network, annotation, ``args``, prefix
and sample name) as command-line arguments, and, so the shell template
doesn't take them for placeholders, ``{pd.__version__!s}`` for
``{pd.__version__}`` in its f-strings.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

DECOUPLER_IMAGE = "community.wave.seqera.io/library/decoupler-py_matplotlib_pandas_python_scanpy:c539733016e15bbc"
# Upstream: process_medium.
DEFAULT_RESOURCES = flyte.Resources(cpu=4, memory="16Gi")

DECOUPLER_PY = r'''
#!/usr/bin/env python3
import os

# This must execute before decoupler is imported, so the env vars are visible when Numba decides whether to cache
os.environ["NUMBA_DISABLE_CACHE"] = "1"
os.environ["NUMBA_CACHE_DIR"] = "./tmp"
os.environ["MPLCONFIGDIR"] = "./tmp"

import argparse
import shlex
import sys

import decoupler as dc  # noqa: E402
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd  # noqa: E402
import scanpy as sc

mat = pd.read_csv(sys.argv[1], sep="\t", index_col=0)
net = pd.read_csv(sys.argv[2], sep="\t")
annot: str = sys.argv[3]


def parse_ext_args(args_string: str):
    """
    Parse external arguments passed from Nextflow's ext.args.
    Expected arguments:
      --min_n <int>
      --transpose <str> (TRUE or FALSE)
      --contrast <str> (optional, e.g., treatment_vs_control)
      --column <str> (Column name to use for transposition; default: log2FoldChange)
      --ensembl_ids <str> (TRUE to convert ENSEMBL IDs to gene symbols, FALSE to skip)
      --methods <str> (Comma-separated list of methods to use (e.g., 'mlm,ulm'))
    """
    if args_string == "null":
        args_string = ""
    args_list = shlex.split(args_string)
    parser = argparse.ArgumentParser()
    parser.add_argument("--min_n", type=int, default=1, help="Minimum n value")
    parser.add_argument("--transpose", type=str, default="FALSE", help="Transpose DESeq2 data if TRUE")
    parser.add_argument("--column", type=str, default="log2FoldChange", help="Column name to use for transposition")
    parser.add_argument("--ensembl_ids", type=str, default="FALSE", help="Convert ENSEMBL IDs to gene symbols if TRUE")
    parser.add_argument(
        "--methods", type=str, default="ulm", help="Comma-separated list of methods to use (e.g., 'mlm,ulm')"
    )
    parser.add_argument("--features_id_col", type=str, default="gene_id", help="Column name for feature IDs")
    parser.add_argument("--features_symbol_col", type=str, default="gene_name", help="Column name for feature symbols")
    parser.add_argument(
        "--round_digits", type=int, default=None, help="Number of digits to round numeric output columns"
    )
    return parser.parse_args(args_list)


# Parse external arguments
raw_args = sys.argv[4]
parsed_args = parse_ext_args(raw_args)
methods = [m.strip() for m in parsed_args.methods.split(",") if m.strip()]

if parsed_args.ensembl_ids.upper() == "TRUE":
    try:
        if not os.path.exists(sys.argv[3]):
            raise FileNotFoundError(f"Annotation file not found: {sys.argv[3]}")

        annot_df = pd.read_csv(sys.argv[3], sep="\t")

        required_cols = {parsed_args.features_id_col, parsed_args.features_symbol_col}

        missing = required_cols - set(annot_df.columns)
        if missing:
            raise ValueError(
                f"Missing required columns in annotation file: {missing}. Available columns: {list(annot_df.columns)}"
            )

        gene_mapping = dict(zip(annot_df[parsed_args.features_id_col], annot_df[parsed_args.features_symbol_col]))
        new_index = [gene_mapping.get(ens, None) for ens in mat.index]
        mat.index = new_index
        mat = mat[mat.index.notnull()]
        mat = mat[~mat.index.duplicated(keep="first")]
    except Exception as e:
        print(f"ERROR: Failed to process annotation file: {e}")
        sys.exit(1)

if parsed_args.transpose.upper() == "TRUE":
    mat = mat[[parsed_args.column]].T.rename(index={parsed_args.column: sys.argv[6]})

parsedargs = {"args": {}}
parsedargs["min_n"] = parsed_args.min_n


results = dc.decouple(mat=mat, net=net, methods=methods, **parsedargs)

for result in results:
    # Save table
    result_df = results[result].copy()
    if parsed_args.round_digits is not None:
        numeric_columns = result_df.select_dtypes(include="number").columns
        result_df[numeric_columns] = result_df[numeric_columns].round(parsed_args.round_digits)

    result_df.to_csv(sys.argv[5] + "_" + result + "_decoupler.tsv", sep="\t")
    contrast_name = result_df.index[0]
    plt.figure(figsize=(8, 6))
    dc.plot_barplot(result_df, contrast_name, top=25, vertical=False)
    plt.savefig(sys.argv[5] + "_" + result + "_decoupler_plot.png", dpi=300, bbox_inches="tight")
    plt.close()

with open("versions.yml", "a") as version_file:
    version_file.write('"DECOUPLER_DECOUPLER":\n')
    version_file.write(f"    decoupler-py: {dc.__version__!s}\n")
    version_file.write(f"    matplotlib: {mpl.__version__!s}\n")
    version_file.write(f"    pandas: {pd.__version__!s}\n")
    version_file.write(f"    python: {sys.version.split()[0]}\n")
    version_file.write(f"    scanpy: {sc.__version__!s}\n")
'''.strip()

decoupler_cmd = shell.create(
    name="decoupler",
    image=DECOUPLER_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"mat": File, "net": File, "annot": File | None, "prefix": str, "sample": str, "args": str},
    defaults={"args": "null"},
    outputs={"results": Dir},
    script=r"""
        MAT=({inputs.mat})
        NET=({inputs.net})
        shopt -s nullglob; ANNOT=({inputs.annot}); shopt -u nullglob
        P={inputs.prefix}
        OUT={outputs.results}
        W=$(mktemp -d); cd "$W"
        ANNOT_PATH=""
        if [ ${#ANNOT[@]} -gt 0 ]; then ANNOT_PATH="${ANNOT[0]}"; fi
        cat > decoupler_run.py <<'END_OF_SCRIPT'
DECOUPLER_PY
END_OF_SCRIPT
        python3 decoupler_run.py "${MAT[0]}" "${NET[0]}" "$ANNOT_PATH" {inputs.args} "$P" {inputs.sample} || exit 1
        mv ./*_decoupler.tsv ./*_decoupler_plot.png versions.yml "$OUT"/
    """.replace("DECOUPLER_PY", DECOUPLER_PY),
)

env = flyte.TaskEnvironment.from_task("decoupler", decoupler_cmd.as_task())


@dataclass
class DecouplerResult:
    estimates: list[File]  # <prefix>_<method>_estimate_decoupler.tsv, by method
    pvals: list[File]  # <prefix>_<method>_pvals_decoupler.tsv
    plots: list[File]  # <prefix>_<method>_estimate_decoupler_plot.png
    versions: File  # versions.yml


async def decoupler(
    mat: File, net: File, prefix: str, sample: str, annot: File | None = None, args: str = ""
) -> DecouplerResult:
    """decoupler's activity estimates for ``mat`` over ``net``; ``sample`` names a transposed column."""
    # As upstream's template renders an unset ext.args.
    results = await decoupler_cmd(mat=mat, net=net, annot=annot, prefix=prefix, sample=sample, args=args or "null")
    files = sorted(await results.list_files(), key=lambda f: f.name)
    versions = [f for f in files if f.name == "versions.yml"]
    if not versions:
        raise FileNotFoundError("decoupler wrote no versions.yml")
    return DecouplerResult(
        estimates=[f for f in files if f.name.endswith("_estimate_decoupler.tsv")],
        pvals=[f for f in files if f.name.endswith("_pvals_decoupler.tsv")],
        plots=[f for f in files if f.name.endswith("_estimate_decoupler_plot.png")],
        versions=versions[0],
    )
