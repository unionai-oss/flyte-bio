"""scanpy — AnnData conversions and cell-hashing demultiplexing of single-cell count matrices.

Exposes:

- :func:`mtx_to_h5ad_simpleaf` — a sample's simpleaf quantification as
  ``<sample>_<input_type>_matrix.h5ad`` (cells and genes sorted, gene IDs
  without version suffixes, the sample recorded in ``obs``).
- :func:`concat_h5ad` — those per-sample files concatenated into
  ``<prefix>_matrix.h5ad``, with each cell annotated from the samplesheet.
- :func:`scanpy_hashsolo` — demultiplex cell-hashing counts with HashSolo
  (``scanpy.external.pp.hashsolo``), from a 10x HTO matrix directory or an
  h5ad: ``<prefix>_assignment_hashsolo.csv`` (the cells' ``obs`` with the
  calls), ``<prefix>_hashsolo.h5ad``, ``<prefix>_params_hashsolo.csv`` and
  ``versions.yml``.

The first two run upstream's scripts (vendored) on upstream's scanpy image.
HashSolo runs upstream's template on its own scanpy image, unchanged but for
taking the templated values (input path, hashing columns, prefix and
``args``) as command-line arguments. Its ``args`` are upstream's
``ext.args`` (``--priors``, ``--pre_existing_clusters``,
``--clustering_data``, ``--number_of_noise_barcodes``, ``--round_digits``).
As upstream, a matrix directory is linked into the working directory under
``matrix_name`` (``hto`` for the test data) and an h5ad keeps its own name;
the parameters table records that name.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.scripts import path

SCANPY_IMAGE = "community.wave.seqera.io/library/scanpy:1.10.2--e83da2205b92a538"
HASHSOLO_IMAGE = "community.wave.seqera.io/library/matplotlib_pandas_python_pyyaml_scanpy:bf7b7ef27120d15c"

DEFAULT_RESOURCES = flyte.Resources(cpu=2, memory="8Gi")

# The scripts point numba's cache at the working directory, so run in a scratch one.
mtx_to_h5ad_simpleaf_cmd = shell.create(
    name="mtx_to_h5ad_simpleaf",
    image=SCANPY_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"script": File, "quant": Dir, "sample": str, "input_type": str},
    outputs={"results": Dir},
    script=r"""
        S=({inputs.script})
        QUANT={inputs.quant}
        W=$(mktemp -d); cd "$W"
        python "${S[0]}" --input "$QUANT" --sample {inputs.sample} \
            --output {outputs.results}/{inputs.sample}_{inputs.input_type}_matrix.h5ad
    """,
)

# Upstream stages the h5ads flat in the working directory and keys each by its
# path there, so link them into a scratch directory and run from it.
concat_h5ad_cmd = shell.create(
    name="concat_h5ad",
    image=SCANPY_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"script": File, "h5ads": list[File], "samplesheet": File, "prefix": str},
    outputs={"results": Dir},
    script=r"""
        S=({inputs.script}); H=({inputs.h5ads}); SHEET=({inputs.samplesheet})
        W=$(mktemp -d); cd "$W"
        for f in "${H[@]}"; do ln -s "$f" "$(basename "$f")"; done
        python "${S[0]}" --samplesheet "${SHEET[0]}" --output {outputs.results}/{inputs.prefix}_matrix.h5ad
    """,
)

HASHSOLO_PY = r'''
#!/usr/bin/env python3

import argparse
import os
import platform
import shlex
import sys

import yaml

os.environ["MPLCONFIGDIR"] = "./tmp/mpl"
os.environ["NUMBA_CACHE_DIR"] = "./tmp/numba"

import matplotlib
import pandas as pd
import scanpy as sc
import scanpy.external as sce


class Arguments:
    # adopted from mygene module (Suzanne Jin)
    """
    Parses the arguments, including the ones coming from the command line.
    """

    def __init__(self) -> None:
        self.data = sys.argv[1]

        self.use_10x = True
        if self.data.endswith(".h5ad"):
            self.use_10x = False

        cell_hashing_columns = sys.argv[2].split()
        self.cell_hashing_columns = [str(x) for x in cell_hashing_columns]

        self.prefix = sys.argv[3]

        self.path_assignment = self.prefix + "_assignment_hashsolo.csv"
        # self.path_plot       = self.prefix + "_hashsolo.jpg"
        self.path_h5ad = self.prefix + "_hashsolo.h5ad"
        self.path_params = self.prefix + "_params_hashsolo.csv"
        self.parse_ext_args(sys.argv[4])

    def parse_ext_args(self, args_string: str) -> None:
        """
        It parses the extended arguments.
        """
        # skip when there are no extended arguments
        if args_string == "null":
            args_string = ""

        # Parse the extended arguments
        args_list = shlex.split(args_string)  # Split the string into a list of arguments
        parser = argparse.ArgumentParser()

        parser.add_argument(
            "--priors",
            type=float,
            nargs=3,
            metavar=("NEGATIVE", "SINGLET", "DOUBLET"),
            help="""List of priors for each hypothesis:
            NEGATIVE = prior for negative hypothesis
            SINGLET  = prior for singlet hypothesis
            DOUBLET  = prior for doublet hypothesis""",
            default=[0.01, 0.8, 0.19],
        )

        parser.add_argument(
            "--pre_existing_clusters",
            help="column in cell_hashing_adata.obs for how to break up demultiplexing",
            type=str,
            default=None,
        )
        parser.add_argument(
            "--clustering_data",
            help="input directory containing transcriptomic data in 10x mtx format.",
            type=str,
            default=None,
        )
        parser.add_argument(
            "--number_of_noise_barcodes",
            help="Number of barcodes to use to create noise distribution",
            type=int,
            default=None,
        )

        parser.add_argument(
            "--round_digits",
            help=(
                "Number of decimal places to round numeric values in cell_hashing_data.obs "
                "before saving. If omitted, no rounding is applied."
            ),
            type=int,
            default=None,
        )

        args = parser.parse_args(args_list)

        # Convert "null" values to default values
        # convert true to True and false to False
        for attr in vars(args):
            value = getattr(args, attr)
            if value == "null":
                setattr(args, attr, parser.get_default(attr))
            elif value == "true":
                setattr(args, attr, True)
            elif value == "false":
                setattr(args, attr, False)

        # Assign args attributes to self attributes
        for attr in vars(args):
            setattr(self, attr, getattr(args, attr))

    def print_args(self) -> None:
        """
        Print the arguments.
        """
        for attr in vars(self):
            print(f"{attr}: {getattr(self, attr)}")


if __name__ == "__main__":
    # ------------------------------ parse and print arguments ------------------------------
    args = Arguments()
    args.print_args()

    # ----------------------------------- read input data -----------------------------------
    if args.use_10x:
        cell_hashing_data = sc.read_10x_mtx(args.data, gex_only=False)

        # Move HTO data from .X to .obs columns
        cell_hashing_data.obs[cell_hashing_data.var_names] = pd.DataFrame(
            # Convert sparse matrix to dense array if needed, otherwise use as-is
            cell_hashing_data.X.toarray() if hasattr(cell_hashing_data.X, "toarray") else cell_hashing_data.X,
            index=cell_hashing_data.obs_names,
            columns=cell_hashing_data.var_names,
        )
        args.cell_hashing_columns = list(cell_hashing_data.obs.columns)
    else:
        cell_hashing_data = sc.read_h5ad(args.data)

    if len(args.cell_hashing_columns) == 2:
        # This edge case issue may be fixed in future versions: https://github.com/calico/solo/issues/91
        args.number_of_noise_barcodes = 1

    # -------------------------------------- call hashsolo -----------------------------------
    if args.clustering_data is not None:
        trans_data = sc.read_10x_mtx(args.clustering_data)
        trans_data.var_names_make_unique()
        sce.pp.hashsolo(
            cell_hashing_data,
            cell_hashing_columns=args.cell_hashing_columns,
            priors=args.priors,
            clustering_data=trans_data,
            pre_existing_clusters=args.pre_existing_clusters,
            number_of_noise_barcodes=args.number_of_noise_barcodes,
        )
    else:
        sce.pp.hashsolo(
            cell_hashing_data,
            cell_hashing_columns=args.cell_hashing_columns,
            priors=args.priors,
            pre_existing_clusters=args.pre_existing_clusters,
            number_of_noise_barcodes=args.number_of_noise_barcodes,
        )

    # ------------------------------------- save results -------------------------------------

    # Round numeric values if requested
    if args.round_digits is not None:
        numeric_cols = cell_hashing_data.obs.select_dtypes(include=["number"]).columns
        cell_hashing_data.obs[numeric_cols] = cell_hashing_data.obs[numeric_cols].round(args.round_digits)

    cell_hashing_data.obs.index.name = "Barcode"

    cell_hashing_data.obs.to_csv(args.path_assignment)

    # plotting does not exist for scanpy.external.pp but we couldn't use the solo package
    # sce.pp.hashsolo.plot_qc_checks_cell_hashing(cell_hashing_data)
    # plt.savefig(args.path_plot, dpi=400)

    cell_hashing_data.write(args.path_h5ad)

    param_list = [[key, getattr(args, key)] for key in sorted(vars(args).keys())]

    param_df = pd.DataFrame(param_list, columns=["Argument", "Value"])
    param_df.fillna("None", inplace=True)
    param_df.to_csv(args.path_params, index=False)

    versions = {
        "SCANPY_HASHSOLO": {
            "python": platform.python_version(),
            "scanpy": sc.__version__,
            "matplotlib": matplotlib.__version__,
            "pandas": pd.__version__,
        }
    }

    with open("versions.yml", "w") as f:
        yaml.dump(versions, f)
'''.strip()

# The script points matplotlib's and numba's caches at the working directory, so run in a scratch one.
scanpy_hashsolo_cmd = shell.create(
    name="scanpy_hashsolo",
    image=HASHSOLO_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "matrix": Dir | None,
        "matrix_name": str,
        "h5ad": File | None,
        "cell_hashing_columns": str,
        "prefix": str,
        "args": str,
    },
    defaults={"matrix_name": "hto", "cell_hashing_columns": "", "args": "null"},
    outputs={"results": Dir},
    script=r"""
        MATRIX={inputs.matrix}
        shopt -s nullglob; H5AD=({inputs.h5ad}); shopt -u nullglob
        P={inputs.prefix}
        OUT={outputs.results}
        W=$(mktemp -d); cd "$W"
        if [ ${#H5AD[@]} -gt 0 ]; then
            DATA=$(basename "${H5AD[0]}"); ln -s "${H5AD[0]}" "$DATA"
        else
            DATA={inputs.matrix_name}; ln -s "$MATRIX" "$DATA"
        fi
        cat > hashsolo.py <<'END_OF_SCRIPT'
HASHSOLO_PY
END_OF_SCRIPT
        python hashsolo.py "$DATA" {inputs.cell_hashing_columns} "$P" {inputs.args} || exit 1
        mv "${P}_assignment_hashsolo.csv" "${P}_hashsolo.h5ad" "${P}_params_hashsolo.csv" versions.yml "$OUT"/
    """.replace("HASHSOLO_PY", HASHSOLO_PY),
)


# One image per environment: HashSolo's differs, so `env` has two.
scanpy_tools_env = flyte.TaskEnvironment.from_task(
    "scanpy_tools",
    mtx_to_h5ad_simpleaf_cmd.as_task(),
    concat_h5ad_cmd.as_task(),
)
scanpy_hashsolo_env = flyte.TaskEnvironment.from_task("scanpy_hashsolo", scanpy_hashsolo_cmd.as_task())
env = flyte.TaskEnvironment(name="scanpy", depends_on=[scanpy_tools_env, scanpy_hashsolo_env])


async def pick(results: Dir, name: str, task: str) -> File:
    f = await results.get_file(name)
    if f is None:
        raise FileNotFoundError(f"{task} wrote no {name}")
    return f


async def mtx_to_h5ad_simpleaf(quant: Dir, sample: str, input_type: str = "raw") -> File:
    """``<sample>_<input_type>_matrix.h5ad`` from a simpleaf ``af_quant`` directory."""
    script = await File.from_local(str(path("mtx_to_h5ad_simpleaf.py")))
    results = await mtx_to_h5ad_simpleaf_cmd(script=script, quant=quant, sample=sample, input_type=input_type)
    return await pick(results, f"{sample}_{input_type}_matrix.h5ad", "mtx_to_h5ad_simpleaf")


async def concat_h5ad(h5ads: list[File], samplesheet: File, prefix: str = "combined_raw") -> File:
    """``<prefix>_matrix.h5ad``: the per-sample ``*_matrix.h5ad`` files concatenated, annotated from the samplesheet.

    Each cell's ``sample`` is its file's name without ``_matrix.h5ad``, as upstream.
    """
    names = [h.path.rstrip("/").rsplit("/", 1)[-1] for h in h5ads]
    if len(set(names)) != len(names):
        raise ValueError(f"concat_h5ad needs distinct file names, got {names}")
    script = await File.from_local(str(path("concat_h5ad.py")))
    results = await concat_h5ad_cmd(script=script, h5ads=h5ads, samplesheet=samplesheet, prefix=prefix)
    return await pick(results, f"{prefix}_matrix.h5ad", "concat_h5ad")


@dataclass
class HashsoloResult:
    assignment: File  # <prefix>_assignment_hashsolo.csv
    h5ad: File  # <prefix>_hashsolo.h5ad
    params: File  # <prefix>_params_hashsolo.csv
    versions: File  # versions.yml


async def scanpy_hashsolo(
    prefix: str,
    matrix: Dir | None = None,
    h5ad: File | None = None,
    cell_hashing_columns: list[str] | None = None,
    matrix_name: str = "hto",
    args: str = "",
) -> HashsoloResult:
    """HashSolo calls for a 10x HTO ``matrix`` directory, or for the ``cell_hashing_columns`` of an ``h5ad``'s ``obs``.

    A matrix directory's hashtags are its features (``cell_hashing_columns`` is ignored, as upstream).
    """
    if (matrix is None) == (h5ad is None):
        raise ValueError("scanpy_hashsolo needs exactly one of matrix and h5ad")
    results = await scanpy_hashsolo_cmd(
        matrix=matrix,
        matrix_name=matrix_name,
        h5ad=h5ad,
        cell_hashing_columns=" ".join(cell_hashing_columns or []),
        prefix=prefix,
        # As upstream's template renders an unset ext.args.
        args=args or "null",
    )
    return HashsoloResult(
        assignment=await pick(results, f"{prefix}_assignment_hashsolo.csv", "scanpy_hashsolo"),
        h5ad=await pick(results, f"{prefix}_hashsolo.h5ad", "scanpy_hashsolo"),
        params=await pick(results, f"{prefix}_params_hashsolo.csv", "scanpy_hashsolo"),
        versions=await pick(results, "versions.yml", "scanpy_hashsolo"),
    )
