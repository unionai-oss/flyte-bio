"""nanocomp — compare long-read datasets (NanoComp).

Exposes :func:`nanocomp`, which runs NanoComp over one or more files of a
single type, chosen as upstream from the first recognised extension in each
name: FASTQ (``fastq``, ``fq``), FASTA (``fasta``, ``fna``, ``ffn``, ``faa``,
``frn``, ``fa``), BAM (``bam``) or an ONT sequencing summary (``txt``). Mixed
or unrecognised types are an error. The files are passed under their own
names, in the caller's order.

The outputs keep NanoComp's names, prefixed ``<prefix>_``:
``<prefix>_NanoComp-report.html``, the per-plot
``<prefix>_NanoComp_*.html`` files (which plots depends on the input type) and
``<prefix>_NanoStats.txt``. Unlike upstream, which declares no output for
them, the weighted histograms (``*_OverlayHistogram_Weighted.html``,
``*_OverlayLogHistogram_Weighted.html``) are kept too. As upstream, NanoComp's dated log isn't kept.
NanoComp's ``--threads`` is the task's CPUs.
"""

import os
from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

NANOCOMP_IMAGE = "quay.io/biocontainers/nanocomp:1.25.6--pyhdfd78af_0"
DEFAULT_RESOURCES = flyte.Resources(cpu=2, memory="8Gi")

FASTQ_EXTENSIONS = ("fq", "fastq")
FASTA_EXTENSIONS = ("fasta", "fna", "ffn", "faa", "frn", "fa")

# A list[File] is staged under the files' names and globbed alphabetically, so
# ``names`` carries the caller's order (one name per line) into the script.
nanocomp_cmd = shell.create(
    name="nanocomp",
    image=NANOCOMP_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"files": list[File], "names": str, "filetype": str, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        IN=({{inputs.files}})
        NAMES={{inputs.names}}
        TYPE={{inputs.filetype}}
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        set -f
        W=$(mktemp -d); cd "$W"
        export HOME="$W" MPLCONFIGDIR="$W/.matplotlib"
        D=$(dirname "${{IN[0]}}")
        FILES=()
        while IFS= read -r NAME; do
            ln -s "$D/$NAME" "$NAME"
            FILES+=("$NAME")
        done <<< "$NAMES"
        NanoComp --"$TYPE" "${{FILES[@]}}" --threads $TASK_CPUS --prefix "${{P}}_" $ARGS || exit 1
        set +f
        shopt -s nullglob
        mv ./*NanoComp*.html ./*NanoStats.txt "$OUT"/
    """,
)

env = flyte.TaskEnvironment.from_task("nanocomp", nanocomp_cmd.as_task())


@dataclass
class NanocompResult:
    report_html: File  # <prefix>_NanoComp-report.html
    html: list[File]  # every NanoComp html file (report and plots), in name order
    stats_txt: File  # <prefix>_NanoStats.txt
    results: Dir  # everything kept


def file_name(file: File) -> str:
    return file.name or os.path.basename(file.path.rstrip("/"))


def input_type(name: str) -> str:
    """NanoComp's input option for a file named ``name``, as upstream picks it."""
    parts = [p for p in name.split(".") if p]
    if len(parts) < 2:
        raise ValueError("Every input file to nanocomp has to have a file ending.")
    for part in parts[1:]:
        if part in FASTQ_EXTENSIONS:
            return "fastq"
        if part in FASTA_EXTENSIONS:
            return "fasta"
        if part == "bam":
            return "bam"
        if part == "txt":
            return "summary"
    raise ValueError(
        f"There was no suitable filetype found for {name}. NanoComp only accepts fasta (fasta, fna, ffn, faa, frn, "
        "fa), fastq (fastq, fq), bam and Nanopore sequencing summary (txt)."
    )


async def nanocomp(files: list[File], prefix: str = "sample", args: str = "") -> NanocompResult:
    """NanoComp comparison of ``files`` (all FASTQ, FASTA, BAM or sequencing summaries)."""
    if not files:
        raise ValueError("There was no suitable filetype found in NanoComp input.")
    names = [file_name(f) for f in files]
    if len(set(names)) != len(names):
        raise ValueError(f"nanocomp inputs need distinct names, got {names}")
    types = {input_type(name) for name in names}
    if len(types) > 1:
        raise ValueError(
            "You gave different filetypes to NanoComp. Please use only *one* of fasta, fastq, bam or Nanopore "
            "sequencing summary."
        )
    results = await nanocomp_cmd(files=files, names="\n".join(names), filetype=types.pop(), prefix=prefix, args=args)
    found = sorted(await results.list_files(), key=file_name)
    by_name = {file_name(f): f for f in found}
    for name in (f"{prefix}_NanoComp-report.html", f"{prefix}_NanoStats.txt"):
        if name not in by_name:
            raise FileNotFoundError(f"NanoComp wrote no {name}")
    return NanocompResult(
        report_html=by_name[f"{prefix}_NanoComp-report.html"],
        html=[f for f in found if file_name(f).endswith(".html")],
        stats_txt=by_name[f"{prefix}_NanoStats.txt"],
        results=results,
    )
