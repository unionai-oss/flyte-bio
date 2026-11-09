"""sigprofiler — SigProfiler, de novo extraction and COSMIC assignment of mutational signatures.

Exposes :func:`sigprofiler`, which pools the joint mutation tables in
``tsv_list`` (mutations with ``NV`` 0 are dropped), builds SBS96, DBS78 and
ID83 count matrices with SigProfilerMatrixGenerator against ``genome`` (by
default installing it first, a download of several GB), and runs
SigProfilerExtractor on each matrix that has mutations. As upstream, the
outputs go under ``results/``: ``input/input_data.txt``, the matrices under
``output/`` and one extraction per context (``SBS96/``, ``DBS78/``,
``ID83/``). They are hundreds of files, so the task returns them as
``results.tar.gz``, with a ``versions.yml``.

``args`` is upstream's ``ext.args`` map as Groovy prints it, e.g.
``[minimum_signatures:1, maximum_signatures:2, nmf_replicates:10,
make_decomposition_plots:False]``.

``genome_installed_path`` is a reference volume packed as a ``.tar.gz`` of
``tsb/<genome>/*.txt``, e.g. from :func:`sigprofiler_install_genome`. When it
is given, the task unpacks it and points SigProfiler at it
(``SIGPROFILERMATRIXGENERATOR_VOLUME``), so the template's install finds the
genome already installed and downloads nothing. That is also the way to pin
a reference: on 2026-09-25 the lab replaced its ``GRCh38.tar.gz`` with a new
build that this image's SigProfilerMatrixGenerator (1.3.5) rejects by
checksum, and moved the old one to ``GRCh38_Legacy``, so GRCh38 now needs
``sigprofiler_install_genome("GRCh38_Legacy", "GRCh38")``. The volume is one
archive rather than a ``Dir`` because the platform's input stager, given a
``Dir`` of a few GB, dies before the task starts (the run then fails with
"OutputUploadFailed: /var/outputs/results not found").
``seeds`` is a
SigProfilerExtractor seeds file (``Seeds.txt``); without one the seeds are
random. SigProfilerExtractor runs a process per core the machine reports
(the template passes no ``--cpu``), so the task's CPU request bounds its
speed, not its process count.

The Python code is upstream's template, written out by an unquoted bash
heredoc: Groovy's ``\\$`` and ``\\\\`` escapes mean the same to bash, so only
the placeholders are renamed (``${task.ext.args ?: ''}`` →
``${task_ext_args}``, ``${tsv_list.join(' ')}`` → ``${tsv_list}``). The
tables are staged as ``1.tsv``, ``2.tsv`` … (upstream stages them as
``*.tsv``).
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

SIGPROFILER_IMAGE = (
    "community.wave.seqera.io/library/"
    "pip_sigprofilerassignment_sigprofilerextractor_sigprofilermatrixgenerator_pruned:02a3f95da35d8c9a"
)
# Upstream's process_high; the genome install needs disk.
DEFAULT_RESOURCES = flyte.Resources(cpu=8, memory="16Gi", disk="40Gi")
INSTALL_RESOURCES = flyte.Resources(cpu=1, memory="2Gi", disk="20Gi")
# Where SigProfilerMatrixGenerator downloads its references from.
GENOME_FTP = "ftp://alexandrovlab-ftp.ucsd.edu/pub/tools/SigProfilerMatrixGenerator"

# Upstream's templates/main_script.py, for an unquoted heredoc (see the module docstring).
SIGPROFILER_TEMPLATE = r"""
#!/usr/bin/env python

import os
import shlex
import shutil
import subprocess
import sys
from importlib.metadata import version

import pandas as pd

# Create writable directories in the Nextflow task working directory
os.makedirs(".matplotlib", exist_ok=True)
os.makedirs("tmp", exist_ok=True)

# Tell Matplotlib to use the writable directory
os.environ["MPLCONFIGDIR"] = os.path.abspath(".matplotlib")


def parse_args(x):
    x = x.strip("[]")
    lexer = shlex.shlex(x, posix=True)
    lexer.whitespace = ","
    lexer.whitespace_split = True
    lexer.commenters = ""
    parsed_args = {}
    for item in lexer:
        if ":" not in item:
            continue
        key, val = item.split(":", 1)
        parsed_args[key.strip()] = val.strip().strip("'").strip('"')
    return parsed_args


# Default options

opt = dict()
opt["prefix"] = "${prefix}"

opt.update(
    {
        "genome": "${genome}",
        "input_type": "matrix",
        "context_type": "96,DINUC,ID",
        "minimum_signatures": 1,
        "maximum_signatures": 25,
        "nmf_replicates": 100,
        "min_nmf_iterations": 10000,
        "max_nmf_iterations": 100000,
        "nmf_test_conv": 10000,
        "seeds": "${seeds}",
        "volume": os.path.abspath("tmp"),
        "make_decomposition_plots": True,
        "download_genome_sigprofiler": True,
        "genome_installed_path": "${genome_installed_path}",
    }
)

# Parse extra args
args_opt = parse_args("${task_ext_args}")
for ao_k, ao_v in args_opt.items():
    if ao_k in {
        "minimum_signatures",
        "maximum_signatures",
        "nmf_replicates",
        "min_nmf_iterations",
        "max_nmf_iterations",
        "nmf_test_conv",
    }:
        opt[ao_k] = int(ao_v)
    elif ao_k in {"make_decomposition_plots", "download_genome_sigprofiler"}:
        val = ao_v.lower()
        if val not in {"true", "false"}:
            raise ValueError(f"Invalid value for '{ao_k}': expected 'true' or 'false', got '{ao_v}'")
        opt[ao_k] = val == "true"
    else:
        opt[ao_k] = ao_v

if not opt["seeds"] or opt["seeds"].strip() in {"", "None", "null"}:
    opt["seeds"] = "random"

print(f"Seeds: {opt['seeds']}", flush=True)
print(f"Volume: {opt['volume']}", flush=True)
print(f"MPLCONFIGDIR: {os.environ['MPLCONFIGDIR']}", flush=True)

# Script


def process_tsv_join(tsv_list):
    patients_tsv = tsv_list.split()
    # Read each file into a pandas DataFrame and ensure all columns are of type 'string'
    tables = []
    for p_table in patients_tsv:
        df = pd.read_csv(p_table, sep="\\t", dtype=str)
        tables.append(df)
    multisample_table = pd.concat(tables, ignore_index=True)
    multisample_table_filtered = multisample_table[multisample_table["NV"] != "0"]
    return multisample_table_filtered


def input_processing(data, prefix, genome):
    new_columns = {"Project": prefix, "Genome": genome, "Type": "SOMATIC", "mut_type": "SNP"}
    df = data.assign(**new_columns)
    df["chr"] = df["chr"].astype(str).str[3:]
    df = df.rename(columns={"Indiv": "Sample", "chr": "chrom", "from": "pos_start", "to": "pos_end"})
    df["ID"] = df["Sample"]
    return df.loc[
        :, ["Project", "Sample", "ID", "Genome", "mut_type", "chrom", "pos_start", "pos_end", "ref", "alt", "Type"]
    ]


if __name__ == "__main__":
    prefix = opt["prefix"]
    genome = opt["genome"]
    input_type = opt["input_type"]
    context_types = opt["context_type"].split(",")
    tsv_list = "${tsv_list}"

    if not os.path.exists(prefix):
        os.mkdir(prefix)

    data = process_tsv_join(tsv_list)
    processed = input_processing(data, prefix, genome)
    processed.to_csv(f"{prefix}/input_data.txt", sep="\t", index=False)

    # Conditionally install genome or use provided path
    if opt.get("download_genome_sigprofiler", True):
        print(f"Installing genome {genome} via SigProfilerMatrixGenerator...", flush=True)
        install_genome = f"SigProfilerMatrixGenerator install {genome} -v {opt['volume']}"
        subprocess.run(install_genome, shell=True, check=True)
    else:
        if not opt.get("genome_installed_path"):
            raise ValueError("download_genome_sigprofiler is False but no genome_installed_path was provided.")
        print(f"Using pre-installed genome at: {opt['genome_installed_path']}", flush=True)

    # Mutation counts matrix generation
    generate_count_matrix = (
        f"SigProfilerMatrixGenerator matrix_generator {prefix} {genome} {prefix} --volume {opt['volume']}"
    )
    subprocess.run(generate_count_matrix, shell=True, check=True)

    matrix_files = {
        "SBS96": os.path.join("output", "SBS", f"{prefix}.SBS96.all"),
        "DBS78": os.path.join("output", "DBS", f"{prefix}.DBS78.all"),
        "ID83": os.path.join("output", "ID", f"{prefix}.ID83.all"),
    }

    context_type_map = {"SBS96": "96", "DBS78": "DINUC", "ID83": "ID"}

    # Run SigProfilerExtractor for each mutation type

    for key, matrix_path in matrix_files.items():
        input_data_path = os.path.join(prefix, matrix_path)

        if not os.path.isfile(input_data_path):
            continue  # Skip this mutation type

        print(f"Running SigProfilerExtractor on {key} using {input_data_path}", flush=True)

        output_dir = os.path.join("results", key)
        context_type = context_type_map[key]

        # Convert Python bools to lowercase strings 'true'/'false'
        make_plots = True if opt["make_decomposition_plots"] else False

        sigprofilerextractor_run = (
            f"SigProfilerExtractor sigprofilerextractor "
            f"--reference_genome  {genome} "
            f"--context_type {context_type} "
            f"--minimum_signatures {opt['minimum_signatures']} "
            f"--maximum_signatures {opt['maximum_signatures']} "
            f"--nmf_replicates {opt['nmf_replicates']} "
            f"--seeds {opt['seeds']} "
            f"--volume {opt['volume']} "
            f"--min_nmf_iterations {opt['min_nmf_iterations']} "
            f"--max_nmf_iterations {opt['max_nmf_iterations']} "
            f"--nmf_test_conv {opt['nmf_test_conv']} "
            f"--make_decomposition_plots {make_plots} "
            f"matrix {output_dir} {input_data_path}"
        )
        subprocess.run(sigprofilerextractor_run, shell=True, check=True)

    # save the output results
    source_dir = prefix + "/"
    dest_dir = "results/"
    shutil.copytree(source_dir, "results", dirs_exist_ok=True)

    # Write version
    python_version = ".".join(map(str, sys.version_info[:3]))
    sig_profiler_assignment_version = version("SigProfilerAssignment")
    sig_profiler_extractor_version = version("SigProfilerExtractor")
    sig_profiler_matrix_generator_version = version("SigProfilerMatrixGenerator")
    pandas_version = pd.__version__
    sig_profiler_plotting_version = version("sigProfilerPlotting")

    with open("versions.yml", "a") as f:
        f.write('"${task_process}":' + "\\n")
        f.write("    python: " + python_version + "\\n")
        f.write("    SigProfilerAssignment: " + sig_profiler_assignment_version + "\\n")
        f.write("    SigProfilerExtractor: " + sig_profiler_extractor_version + "\\n")
        f.write("    SigProfilerMatrixGenerator: " + sig_profiler_matrix_generator_version + "\\n")
        f.write("    pandas: " + pandas_version + "\\n")
        f.write("    sigProfilerPlotting: " + sig_profiler_plotting_version + "\\n")
"""

sigprofiler_cmd = shell.create(
    name="sigprofiler",
    image=SIGPROFILER_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "tsv_list": list[File],
        "genome": str,
        "genome_installed_path": File | None,
        "seeds": File | None,
        "prefix": str,
        "args": str,
    },
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        TSV=({inputs.tsv_list})
        shopt -s nullglob; SEEDS=({inputs.seeds}); shopt -u nullglob
        shopt -s nullglob; GENOME_TAR=({inputs.genome_installed_path}); shopt -u nullglob
        genome={inputs.genome}
        prefix={inputs.prefix}
        task_ext_args={inputs.args}
        task_process=SIGPROFILER
        W=$(mktemp -d); cd "$W"
        tsv_list=""
        for I in "${!TSV[@]}"; do
            ln -s "${TSV[$I]}" "$((I + 1)).tsv"
            tsv_list="$tsv_list${tsv_list:+ }$((I + 1)).tsv"
        done
        seeds=""
        if [ ${#SEEDS[@]} -gt 0 ]; then
            ln -s "${SEEDS[0]}" .
            seeds=$(basename "${SEEDS[0]}")
        fi
        genome_installed_path=""
        if [ ${#GENOME_TAR[@]} -gt 0 ]; then
            # SigProfiler looks for references under $SIGPROFILERMATRIXGENERATOR_VOLUME/tsb/<genome>/ first.
            mkdir genome_installed
            tar -xzf "${GENOME_TAR[0]}" -C genome_installed || exit 1
            genome_installed_path="$W/genome_installed"
            export SIGPROFILERMATRIXGENERATOR_VOLUME="$genome_installed_path"
        fi
        cat > main_script.py <<END_OF_TEMPLATE"""
    + SIGPROFILER_TEMPLATE
    + r"""END_OF_TEMPLATE
        python main_script.py || exit 1
        tar -czf {outputs.results}/results.tar.gz results || exit 1
        mv versions.yml {outputs.results}/
    """,
)

# Streams <archive>.tar.gz from the lab's FTP server into tsb/, renamed to <install_as>, checks the files
# against the checksums of the image's SigProfilerMatrixGenerator, and packs tsb/ into one archive.
sigprofiler_install_genome_cmd = shell.create(
    name="sigprofiler_install_genome",
    image=SIGPROFILER_IMAGE,
    resources=INSTALL_RESOURCES,
    inputs={"archive": str, "install_as": str},
    outputs={"volume": File},
    script=rf"""
        ARCHIVE={{inputs.archive}}
        NAME={{inputs.install_as}}
        W=$(mktemp -d)
        mkdir -p "$W/tsb"; cd "$W/tsb"
        python - "{GENOME_FTP}/$ARCHIVE.tar.gz" <<'END_OF_DOWNLOAD' || exit 1
import sys
import tarfile
import urllib.request

with urllib.request.urlopen(sys.argv[1]) as response, tarfile.open(fileobj=response, mode="r|gz") as tar:
    tar.extractall(".", filter="data")
END_OF_DOWNLOAD
        if [ "$ARCHIVE" != "$NAME" ]; then mv "$ARCHIVE" "$NAME" || exit 1; fi
        SIGPROFILERMATRIXGENERATOR_VOLUME="$W" python - "$NAME" <<'END_OF_CHECK' || exit 1
import sys

from SigProfilerMatrixGenerator.scripts.reference_genome_manager import ReferenceGenomeManager

manager = ReferenceGenomeManager()
if not manager.is_genome_installed(sys.argv[1]):
    manager.print_genome_checksum_verification_report(sys.argv[1])
    sys.exit(f"{{sys.argv[1]}} doesn't match this SigProfilerMatrixGenerator's checksums")
END_OF_CHECK
        tar -czf {{outputs.volume}} -C "$W" tsb || exit 1
    """,
)

env = flyte.TaskEnvironment.from_task(
    "sigprofiler", sigprofiler_cmd.as_task(), sigprofiler_install_genome_cmd.as_task()
)


async def sigprofiler_install_genome(archive: str, install_as: str = "") -> File:
    """A reference volume for :func:`sigprofiler`'s ``genome_installed_path``: ``tsb/<install_as>/*.txt``, as .tar.gz.

    ``archive`` names a reference on the lab's FTP server (``GRCh38``, ``GRCh38_Legacy``, ``GRCh37`` …);
    ``install_as`` (default: ``archive``) is the genome name to give :func:`sigprofiler`. The files must match
    the checksums the image's SigProfilerMatrixGenerator (1.3.5) knows for that name.
    """
    return await sigprofiler_install_genome_cmd(archive=archive, install_as=install_as or archive)


@dataclass
class SigprofilerResult:
    results_tar: File  # results.tar.gz: results/{input,output,SBS96,DBS78,ID83}/…
    versions: File
    results: Dir


async def sigprofiler(
    tsv_list: list[File],
    genome: str,
    genome_installed_path: File | None = None,
    seeds: File | None = None,
    prefix: str = "sample",
    args: str = "",
) -> SigprofilerResult:
    """Extract and assign mutational signatures from the joint tables ``tsv_list`` (``genome``: e.g. GRCh38)."""
    results = await sigprofiler_cmd(
        tsv_list=tsv_list,
        genome=genome,
        genome_installed_path=genome_installed_path,
        seeds=seeds,
        prefix=prefix,
        args=args,
    )

    async def pick(name: str) -> File:
        found = await results.get_file(name)
        if found is None:
            raise FileNotFoundError(f"sigprofiler wrote no {name}")
        return found

    return SigprofilerResult(
        results_tar=await pick("results.tar.gz"), versions=await pick("versions.yml"), results=results
    )
