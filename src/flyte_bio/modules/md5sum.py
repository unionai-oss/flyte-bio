"""md5sum — MD5 checksums of files.

Exposes :func:`md5sum`, which checksums a list of files, as upstream:

- with ``as_separate_files`` (the default): one ``<file name>.md5`` per file;
- otherwise one ``<prefix>.md5`` listing every file, in name order (upstream's
  ``find -L *``).

The files are linked into a working directory under their own names first,
so each line names the file as upstream's does. Files named ``*.md5`` are
skipped, as upstream. ``args`` (e.g. ``--tag`` for BSD-style lines) is passed
to ``md5sum`` and split on whitespace.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

MD5SUM_IMAGE = "community.wave.seqera.io/library/coreutils_grep_gzip_lbzip2_pruned:838ba80435a629f8"
DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory="2Gi")

md5sum_cmd = shell.create(
    name="md5sum",
    image=MD5SUM_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"files": list[File], "as_separate_files": bool, "prefix": str, "args": str},
    defaults={"as_separate_files": True, "args": ""},
    outputs={"results": Dir},
    script=r"""
        F=({inputs.files})
        P={inputs.prefix}
        ARGS={inputs.args}
        W=$(mktemp -d); cd "$W"
        for f in "${F[@]}"; do ln -s "$f" "$(basename "$f")"; done
        if [ {inputs.as_separate_files} = true ]; then
            find -L * -maxdepth 0 -type f ! -name '*.md5' \
                -exec sh -c 'md5sum '"$ARGS"' "$1" > "$1.md5"' _ "{}" \;
        else
            find -L * -type f ! -name '*.md5' -exec md5sum $ARGS "{}" + > "$P.md5"
        fi
        for f in "${F[@]}"; do rm "$(basename "$f")"; done
        mv ./*.md5 {outputs.results}/
    """,
)

env = flyte.TaskEnvironment.from_task("md5sum", md5sum_cmd.as_task())


@dataclass
class Md5sumResult:
    checksums: list[File]  # one <name>.md5 per file in list order, or the one <prefix>.md5
    results: Dir


def file_name(file: File) -> str:
    return file.path.rstrip("/").rsplit("/", 1)[-1]


async def md5sum(
    files: list[File], as_separate_files: bool = True, prefix: str = "sample", args: str = ""
) -> Md5sumResult:
    """MD5 checksums of ``files``: ``<name>.md5`` each, or ``<prefix>.md5`` for all."""
    names = [file_name(f) for f in files]
    if not files:
        raise ValueError("md5sum needs at least one file")
    if len(set(names)) != len(names):
        raise ValueError(f"md5sum inputs need distinct names, got {names}")
    results = await md5sum_cmd(files=files, as_separate_files=as_separate_files, prefix=prefix, args=args)
    wanted = [f"{n}.md5" for n in names if not n.endswith(".md5")] if as_separate_files else [f"{prefix}.md5"]
    checksums = []
    for name in wanted:
        found = await results.get_file(name)
        if found is None:
            raise FileNotFoundError(f"md5sum wrote no {name}")
        checksums.append(found)
    return Md5sumResult(checksums=checksums, results=results)
