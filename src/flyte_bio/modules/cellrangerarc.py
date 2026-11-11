"""cellrangerarc — 10x Genomics Cell Ranger ARC: single-cell multiome (ATAC + gene expression).

Exposes:

- :func:`cellrangerarc_mkgtf`, which filters a GTF by attributes
  (``cellranger-arc mkgtf``; pass ``--attribute=gene_biotype:...`` in
  ``args``) into ``<prefix>.gtf``.
- :func:`cellrangerarc_mkref`, which builds a Cell Ranger ARC reference
  (``cellranger-arc mkref``: BWA and STAR indexes, ``genes/``, ``regions/``,
  ``reference.json``) from a genome FASTA, a GTF and optional motifs. As
  upstream, the inputs are staged under their own names and named in a
  config (``config``): the one given, or one written from the inputs (genome
  ``<reference_name>``, organism the FASTA's name up to its first dot, and
  ``args`` appended). The reference comes back as
  ``<reference_name>.tar.gz`` (an archive of the ``<reference_name>``
  folder), with the config used. A given config must name the genome
  ``<reference_name>``.
- :func:`cellrangerarc_count`, which analyses one multiome sample
  (``cellranger-arc count``) against such a reference, extracting it first.
  As upstream, it writes the libraries CSV ``<prefix>_lib.csv`` from the
  samples' library types (``gex`` is Gene Expression, anything else
  Chromatin Accessibility) and names; the FASTQs are staged under their own
  names, so name them as ``bcl2fastq`` does
  (``<sample>_S1_L001_R1_001.fastq.gz``). Cell Ranger ARC's ``<prefix>/outs``
  folder comes back as ``<prefix>.outs.tar.gz``, with ``summary.csv``,
  ``web_summary.html`` and the libraries CSV alongside.

count passes the task's cores and memory as ``--localcores`` and
``--localmem``; mkref passes the cores as ``--nthreads``. Telemetry is
disabled (``TENX_DISABLE_TELEMETRY``) and ``HOME`` is the task's working
directory.

The container is upstream's build of 10x's Cell Ranger ARC 2.0.2,
distributed under 10x Genomics' end-user licence: the software may be used
only in connection with 10x Genomics products or data generated with them,
and not redistributed.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

CELLRANGERARC_IMAGE = "quay.io/nf-core/cellranger-arc:2.0.2"
NO_TELEMETRY = {"TENX_DISABLE_TELEMETRY": "1"}
# Upstream: process_low.
MKGTF_RESOURCES = flyte.Resources(cpu=1, memory="4Gi", disk="10Gi")
# Upstream: process_medium. A mouse or human genome needs about 32 GiB for STAR's index.
MKREF_RESOURCES = flyte.Resources(cpu=4, memory="16Gi", disk="50Gi")
# Upstream: process_high (12 CPUs, 72 GB). 10x recommends 64 GB or more for
# real samples: raise with cellrangerarc_count_cmd.override(resources=...).
COUNT_RESOURCES = flyte.Resources(cpu=4, memory="16Gi", disk="50Gi")

# Upstream's inline Python, unchanged but for taking its templated values
# (FASTA, GTF, motifs and config names, args, prefix) as command-line
# arguments, and Groovy's escaping undone.
MKREF_CONFIG_PY = r"""
import sys
from os.path import exists
import shutil

fasta, gtf, motifs, add, reference_config, prefix = sys.argv[1:7]

if ( reference_config == "[]" ):

    config = open("config", "w")
    config.write("{\n")
    config.write('\torganism: "{}"\n'.format(fasta.split(".")[0]))
    config.write('\tgenome: ["' + prefix + '"]\n')
    config.write('\tinput_fasta: ["{}"]\n'.format(fasta))
    config.write('\tinput_gtf: ["{}"]\n'.format(gtf))
    if motifs != "[]":
        config.write('\tinput_motifs: "{}"\n'.format(motifs))
    if add != None:
        config.write(add + "\n")
    config.write("}")
    config.close()

    print("Wrote config file")
else:
    if ( not exists("config") ):
        shutil.move(reference_config, "config")
""".strip()

# Upstream's inline Python, unchanged but for taking its templated values
# (library types and sample names, joined by commas; the libraries CSV's name;
# the FASTQ folder) as command-line arguments, and Groovy's escaping undone.
LIBRARIES_CSV_PY = r"""
import sys

sample_types = sys.argv[1].split(",")
sample_names = sys.argv[2].split(",")
lib_csv_name = sys.argv[3]
fastq_folder = sys.argv[4]
unique_samples_names = set(sample_names)

lib_csv = open(lib_csv_name, "w")
lib_csv.write("fastqs,sample,library_type")

for i in range(0, len(sample_types)):
    if sample_names[i] in unique_samples_names:
        unique_samples_names.remove(
            sample_names[i]
        )  # this has to be done to account for different Lane files (e.g., L002)
        if sample_types[i] == "gex":
            lib_csv.write("\n{},{},{}".format(fastq_folder, sample_names[i], "Gene Expression"))
        else:
            lib_csv.write("\n{},{},{}".format(fastq_folder, sample_names[i], "Chromatin Accessibility"))

lib_csv.close()

print("Wrote lib.csv file to {}".format(lib_csv_name))
""".strip()

cellrangerarc_mkgtf_cmd = shell.create(
    name="cellrangerarc_mkgtf",
    image=CELLRANGERARC_IMAGE,
    resources=MKGTF_RESOURCES,
    env_vars=NO_TELEMETRY,
    inputs={"gtf": File, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        GTF=({inputs.gtf})
        P={inputs.prefix}
        ARGS={inputs.args}
        OUT={outputs.results}
        W=$(mktemp -d); cd "$W"
        export HOME="$W"
        cellranger-arc \
            mkgtf \
            "${GTF[0]}" \
            "$P.gtf" \
            $ARGS > run.log 2>&1 || { tail -c 3000 run.log >&2; exit 1; }
        mv "$P.gtf" "$OUT"/
    """,
)

cellrangerarc_mkref_cmd = shell.create(
    name="cellrangerarc_mkref",
    image=CELLRANGERARC_IMAGE,
    resources=MKREF_RESOURCES,
    env_vars=NO_TELEMETRY,
    inputs={
        "fasta": File,
        "gtf": File,
        "motifs": File | None,
        "reference_config": File | None,
        "reference_name": str,
        "args": str,
    },
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        shopt -s nullglob
        FA=({{inputs.fasta}}); GTF=({{inputs.gtf}}); MOTIFS=({{inputs.motifs}}); CONFIG=({{inputs.reference_config}})
        shopt -u nullglob
        P={{inputs.reference_name}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        export HOME="$W"

        # Every input under its own name in the working directory, as upstream stages them.
        name_or_empty() {{ if [ $# -gt 0 ]; then ln -s "$1" .; basename "$1"; else echo "[]"; fi; }}
        FASTA_NAME=$(name_or_empty "${{FA[@]}}")
        GTF_NAME=$(name_or_empty "${{GTF[@]}}")
        MOTIFS_NAME=$(name_or_empty "${{MOTIFS[@]}}")
        CONFIG_NAME=$(name_or_empty "${{CONFIG[@]}}")

        cat > write_config.py <<'END_OF_SCRIPT'
{MKREF_CONFIG_PY}
END_OF_SCRIPT
        python3 write_config.py "$FASTA_NAME" "$GTF_NAME" "$MOTIFS_NAME" "$ARGS" "$CONFIG_NAME" "$P" || exit 1

        cellranger-arc \
            mkref \
            --config=config \
            --nthreads=$TASK_CPUS \
            $ARGS > run.log 2>&1 || {{ tail -c 3000 run.log >&2; exit 1; }}

        if [ ! -d "$P" ]; then
            echo "cellranger-arc mkref wrote no $P/: the config must name the genome $P" >&2
            exit 1
        fi
        tar -czf "$OUT/$P.tar.gz" "$P"
        cp -L config "$OUT"/
    """,
)

cellrangerarc_count_cmd = shell.create(
    name="cellrangerarc_count",
    image=CELLRANGERARC_IMAGE,
    resources=COUNT_RESOURCES,
    env_vars=NO_TELEMETRY,
    inputs={
        "reads": list[File],
        "reference": File,
        "sample_types": str,
        "sample_names": str,
        "prefix": str,
        "args": str,
    },
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        READS=({{inputs.reads}})
        REF_ARCHIVE=({{inputs.reference}})
        TYPES={{inputs.sample_types}}
        NAMES={{inputs.sample_names}}
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        export HOME="$W"

        mkdir fastqs reference_archive
        for R in "${{READS[@]}}"; do ln -s "$R" fastqs/; done
        tar -xzf "${{REF_ARCHIVE[0]}}" -C reference_archive || exit 1
        if [ -f reference_archive/reference.json ]; then
            mv reference_archive reference
        else
            INNER=(reference_archive/*/); mv "${{INNER[0]%/}}" reference
        fi

        fastq_folder=$(readlink -f fastqs)
        cat > write_libraries.py <<'END_OF_SCRIPT'
{LIBRARIES_CSV_PY}
END_OF_SCRIPT
        python3 write_libraries.py "$TYPES" "$NAMES" "${{P}}_lib.csv" "$fastq_folder" || exit 1

        cellranger-arc \
            count \
            --id="$P" \
            --libraries="${{P}}_lib.csv" \
            --reference=reference \
            --localcores=$TASK_CPUS \
            --localmem=$(( TASK_MEMORY_MB / 1024 )) \
            $ARGS > run.log 2>&1 || {{ tail -c 3000 run.log >&2; exit 1; }}

        tar -czhf "$OUT/$P.outs.tar.gz" "$P/outs"
        cp "${{P}}_lib.csv" "$OUT"/
        for F in summary.csv web_summary.html; do
            if [ -e "$P/outs/$F" ]; then cp -L "$P/outs/$F" "$OUT"/; fi
        done
    """,
)

env = flyte.TaskEnvironment.from_task(
    "cellrangerarc",
    cellrangerarc_mkgtf_cmd.as_task(),
    cellrangerarc_mkref_cmd.as_task(),
    cellrangerarc_count_cmd.as_task(),
)


async def pick_file(results: Dir, name: str, tool: str) -> File:
    found = await results.get_file(name)
    if found is None:
        raise FileNotFoundError(f"{tool} wrote no {name}")
    return found


async def cellrangerarc_mkgtf(gtf: File, prefix: str, args: str = "") -> File:
    """``gtf`` filtered by ``cellranger-arc mkgtf`` (attributes in ``args``) into ``<prefix>.gtf``."""
    results = await cellrangerarc_mkgtf_cmd(gtf=gtf, prefix=prefix, args=args)
    return await pick_file(results, f"{prefix}.gtf", "cellranger-arc mkgtf")


@dataclass
class ArcReference:
    reference: File  # <reference_name>.tar.gz, an archive of the <reference_name> folder
    config: File  # the mkref config used


async def cellrangerarc_mkref(
    fasta: File,
    gtf: File,
    reference_name: str,
    motifs: File | None = None,
    reference_config: File | None = None,
    args: str = "",
) -> ArcReference:
    """A Cell Ranger ARC reference from ``fasta``, ``gtf`` and ``motifs``, or as ``reference_config`` says."""
    results = await cellrangerarc_mkref_cmd(
        fasta=fasta,
        gtf=gtf,
        motifs=motifs,
        reference_config=reference_config,
        reference_name=reference_name,
        args=args,
    )
    tool = "cellranger-arc mkref"
    return ArcReference(
        reference=await pick_file(results, f"{reference_name}.tar.gz", tool),
        config=await pick_file(results, "config", tool),
    )


@dataclass
class ArcCountResult:
    archive: File  # <prefix>.outs.tar.gz, an archive of <prefix>/outs
    libraries: File  # <prefix>_lib.csv
    summary: File  # outs/summary.csv
    web_summary: File  # outs/web_summary.html
    results: Dir  # everything above


async def cellrangerarc_count(
    reads: list[File],
    reference: File,
    sample_types: list[str],
    sample_names: list[str],
    prefix: str,
    args: str = "",
) -> ArcCountResult:
    """Cell Ranger ARC analysis of one multiome sample's ``reads`` against ``reference`` (a ``.tar.gz``).

    ``sample_types`` (``gex`` or e.g. ``atac``) and ``sample_names`` (the
    FASTQs' sample names) pair up, one entry per library.
    """
    if len(sample_types) != len(sample_names) or not sample_types:
        raise ValueError("cellrangerarc_count needs one sample type per sample name")
    if not reads:
        raise ValueError("cellrangerarc_count needs the sample's FASTQs")
    results = await cellrangerarc_count_cmd(
        reads=reads,
        reference=reference,
        sample_types=",".join(sample_types),
        sample_names=",".join(sample_names),
        prefix=prefix,
        args=args,
    )
    tool = "cellranger-arc count"
    return ArcCountResult(
        archive=await pick_file(results, f"{prefix}.outs.tar.gz", tool),
        libraries=await pick_file(results, f"{prefix}_lib.csv", tool),
        summary=await pick_file(results, "summary.csv", tool),
        web_summary=await pick_file(results, "web_summary.html", tool),
        results=results,
    )
