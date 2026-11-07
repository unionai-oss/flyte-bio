"""bioawk — run a caller-supplied bioawk program over sequence files.

Exposes :func:`bioawk`, awk with bioawk's parsers for FASTA/FASTQ, BED, SAM,
VCF and GFF (``-c fastx`` in ``args`` names the fields ``$name``, ``$seq``,
…). The caller gives the program, as text (upstream's ``ext.args2``) or as a
``program_file`` (which wins when both are set), and the files to run it over,
which bioawk reads gzipped or not.

As upstream:

- the output goes to ``<prefix>.<extension>``, through ``gzip`` when the
  extension ends in ``gz``. With ``redirect_output=False`` nothing is
  redirected and the program writes its own output files
  (``print > "out.tsv"``);
- every ``*.<extension>`` file the step writes is an output;
- an input named ``<prefix>.<extension>`` is an error.

Unlike upstream, ``program`` is the awk program itself, passed to bioawk as one
argument, not shell text spliced into the command: don't wrap it in quotes.
``args`` (bioawk options such as ``-c fastx``) is split on whitespace.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

BIOAWK_IMAGE = "community.wave.seqera.io/library/bioawk:1.0--63ed57d97ec41790"
DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory="6Gi")

# A list[File] is staged under the files' names and globbed alphabetically, so
# ``names`` carries the caller's order (one name per line) into the script.
bioawk_cmd = shell.create(
    name="bioawk",
    image=BIOAWK_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "files": list[File],
        "names": str,
        "program": str,
        "program_file": File | None,
        "prefix": str,
        "extension": str,
        "redirect_output": bool,
        "args": str,
    },
    defaults={"program": "", "redirect_output": True, "args": ""},
    outputs={"results": Dir},
    script=r"""
        set -o pipefail
        shopt -s nullglob; IN=({inputs.files}); PROGRAM_FILE=({inputs.program_file}); shopt -u nullglob
        NAMES={inputs.names}
        PROGRAM={inputs.program}
        P={inputs.prefix}
        S={inputs.extension}
        ARGS={inputs.args}
        OUT={outputs.results}
        set -f
        D=$(dirname "${IN[0]}")
        W=$(mktemp -d); cd "$W"

        FILES=()
        while IFS= read -r NAME; do ln -s "$D/$NAME" "$NAME"; FILES+=("$NAME"); done <<< "$NAMES"
        if [ ${#PROGRAM_FILE[@]} -gt 0 ]; then PROG=(-f "${PROGRAM_FILE[0]}"); else PROG=("$PROGRAM"); fi

        if [ {inputs.redirect_output} = false ]; then
            bioawk $ARGS "${PROG[@]}" "${FILES[@]}"
        elif [[ $S == *gz ]]; then
            bioawk $ARGS "${PROG[@]}" "${FILES[@]}" | gzip > "$P.$S"
        else
            bioawk $ARGS "${PROG[@]}" "${FILES[@]}" > "$P.$S"
        fi

        rm -f "${FILES[@]}"
        set +f; shopt -s nullglob
        WRITTEN=(./*."$S")
        if [ ${#WRITTEN[@]} -gt 0 ]; then mv "${WRITTEN[@]}" "$OUT"/; fi
    """,
)

env = flyte.TaskEnvironment.from_task("bioawk", bioawk_cmd.as_task())


@dataclass
class BioawkResult:
    output: File | None  # <prefix>.<extension>; None if the program wrote only other names
    outputs: list[File]  # every *.<extension> file written, as upstream's output channel
    results: Dir


def file_name(file: File) -> str:
    return file.path.rstrip("/").rsplit("/", 1)[-1]


async def bioawk(
    files: list[File],
    extension: str,
    program: str = "",
    program_file: File | None = None,
    prefix: str = "sample",
    redirect_output: bool = True,
    args: str = "",
) -> BioawkResult:
    """Run bioawk over ``files`` in list order, writing ``<prefix>.<extension>``.

    ``program_file`` is used when given, else ``program``.
    """
    names = [file_name(f) for f in files]
    if not names:
        raise ValueError("bioawk needs at least one input file")
    if len(set(names)) != len(names):
        raise ValueError(f"bioawk inputs need distinct names, got {names}")
    if f"{prefix}.{extension}" in names:
        raise ValueError('Input and output names are the same, use "task.ext.prefix" to disambiguate.')
    if program_file is None and not program:
        raise ValueError("bioawk needs a program or a program_file")
    results = await bioawk_cmd(
        files=files,
        names="\n".join(names),
        program=program,
        program_file=program_file,
        prefix=prefix,
        extension=extension,
        redirect_output=redirect_output,
        args=args,
    )
    outputs = sorted([f for f in await results.list_files() if file_name(f).endswith(f".{extension}")], key=file_name)
    output = next((f for f in outputs if file_name(f) == f"{prefix}.{extension}"), None)
    return BioawkResult(output=output, outputs=outputs, results=results)
