"""gawk — run a caller-supplied awk program over files.

Exposes :func:`gawk`, a generic awk step: the caller gives the program, as
text (upstream's ``ext.args2``) or as a ``program_file`` (which wins when both
are set), and the files to run it over.

As upstream:

- inputs ending in ``.gz`` are decompressed first and awk reads them under
  their names without ``.gz``; other inputs keep their names, so ``FILENAME``
  is what upstream's program would see;
- awk's output goes to ``<prefix>.<suffix>``, through ``gzip`` when the suffix
  ends in ``gz``. With ``redirect_output=False`` nothing is redirected and the
  program writes its own output files (``print > "out.bed"``);
- the suffix defaults to the first input's extension;
- every ``*.<suffix>`` file the step writes is an output;
- an input named ``<prefix>.<suffix>`` is an error.

Unlike upstream, ``program`` is the awk program itself, passed to awk as one
argument, not shell text spliced into the command: don't wrap it in quotes.
``args`` (awk options such as ``-v n=1``) is split on whitespace.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

GAWK_IMAGE = "community.wave.seqera.io/library/gawk:5.3.1--e09efb5dfc4b8156"
DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory="6Gi")

# A list[File] is staged under the files' names and globbed alphabetically, so
# ``names`` carries the caller's order (one name per line) into the script.
gawk_cmd = shell.create(
    name="gawk",
    image=GAWK_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "files": list[File],
        "names": str,
        "program": str,
        "program_file": File | None,
        "prefix": str,
        "suffix": str,
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
        S={inputs.suffix}
        ARGS={inputs.args}
        OUT={outputs.results}
        set -f
        W=$(mktemp -d); cd "$W"

        FILES=()
        if [ ${#IN[@]} -gt 0 ]; then
            D=$(dirname "${IN[0]}")
            while IFS= read -r NAME; do
                if [[ $NAME == *.gz ]]; then
                    gzip -cd "$D/$NAME" > "${NAME%.gz}"
                    FILES+=("${NAME%.gz}")
                else
                    ln -s "$D/$NAME" "$NAME"
                    FILES+=("$NAME")
                fi
            done <<< "$NAMES"
        else
            exec < /dev/null
        fi
        if [ ${#PROGRAM_FILE[@]} -gt 0 ]; then PROG=(-f "${PROGRAM_FILE[0]}"); else PROG=("$PROGRAM"); fi

        if [ {inputs.redirect_output} = false ]; then
            awk $ARGS "${PROG[@]}" "${FILES[@]}"
        elif [[ $S == *gz ]]; then
            awk $ARGS "${PROG[@]}" "${FILES[@]}" | gzip > "$P.$S"
        else
            awk $ARGS "${PROG[@]}" "${FILES[@]}" > "$P.$S"
        fi

        if [ ${#FILES[@]} -gt 0 ]; then rm -f "${FILES[@]}"; fi
        set +f; shopt -s nullglob
        WRITTEN=(./*."$S")
        if [ ${#WRITTEN[@]} -gt 0 ]; then mv "${WRITTEN[@]}" "$OUT"/; fi
    """,
)

env = flyte.TaskEnvironment.from_task("gawk", gawk_cmd.as_task())


@dataclass
class GawkResult:
    output: File | None  # <prefix>.<suffix>; None if the program wrote only other names
    outputs: list[File]  # every *.<suffix> file written, as upstream's output channel
    results: Dir


def file_name(file: File) -> str:
    return file.path.rstrip("/").rsplit("/", 1)[-1]


async def gawk(
    files: list[File],
    program: str = "",
    program_file: File | None = None,
    prefix: str = "sample",
    suffix: str = "",
    redirect_output: bool = True,
    args: str = "",
) -> GawkResult:
    """Run awk over ``files`` in list order, writing ``<prefix>.<suffix>``.

    ``program_file`` is used when given, else ``program``. ``suffix`` defaults
    to the first file's extension.
    """
    names = [file_name(f) for f in files]
    if len(set(names)) != len(names):
        raise ValueError(f"gawk inputs need distinct names, got {names}")
    if not suffix:
        if not files:
            raise ValueError("suffix not set and no input provided")
        suffix = names[0].rsplit(".", 1)[-1] if "." in names[0] else ""
    if f"{prefix}.{suffix}" in names:
        raise ValueError("Input and output names are the same, set prefix in module configuration to disambiguate!")
    if program_file is None and not program:
        raise ValueError("gawk needs a program or a program_file")
    results = await gawk_cmd(
        files=files,
        names="\n".join(names),
        program=program,
        program_file=program_file,
        prefix=prefix,
        suffix=suffix,
        redirect_output=redirect_output,
        args=args,
    )
    outputs = [f for f in await results.list_files() if file_name(f).endswith(f".{suffix}")]
    output = next((f for f in outputs if file_name(f) == f"{prefix}.{suffix}"), None)
    return GawkResult(output=output, outputs=outputs, results=results)
