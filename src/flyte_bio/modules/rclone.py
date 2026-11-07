"""rclone — move and verify data between storage remotes.

Exposes:

- :func:`rclone_check` — compare a source and a destination by size and hash.
- :func:`rclone_checksum` — check a destination against a sum file (e.g.
  ``md5sum`` output).
- :func:`rclone_copy` — copy a source to a destination, writing rclone's log.

Sources and destinations are rclone paths: ``remote:path``, a local path, or,
as upstream, a URL: ``s3://bucket/key`` becomes ``s3:bucket/key`` and an
``http(s)://host/path`` URL is read with rclone's HTTP backend
(``--http-url https://host :http:path``). Remotes come from an optional rclone
config file (``config``), or rclone's ``RCLONE_CONFIG_*`` environment.

A task runs in its own pod, so local data has to be handed to it: ``staged``,
a Dir, is reachable in the task as the remote ``staged:`` (an rclone alias of
the staged directory; read-only).

:func:`rclone_check` and :func:`rclone_checksum` keep upstream's names for
their reports, and, as upstream, don't write the ones that come out empty:

- ``<prefix>.combined.txt`` — every file, with its status (``=`` match,
  ``*`` differ, ``+`` only on the source, ``-`` only on the destination,
  ``!`` error);
- ``<prefix>.differ.txt``, ``<prefix>.missing_on_dst.txt``,
  ``<prefix>.missing_on_src.txt``, ``<prefix>.match.txt``,
  ``<prefix>.error.txt`` — the paths with each status;
- ``<prefix>.exit_code.txt`` — rclone's exit code (``1`` when anything
  differs). A difference doesn't fail the task.

As upstream, the reports are sorted, ``combined`` by path; sorting is
byte-wise (``LC_ALL=C``), so the order doesn't depend on the image's locale.
rclone checks with as many checkers as the task has CPUs (and copies with
half as many transfers).
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

RCLONE_IMAGE = "community.wave.seqera.io/library/rclone:1.74.3--2ef33c5b9132aa97"
DEFAULT_RESOURCES = flyte.Resources(cpu=2, memory="2Gi")
REPORTS = ("combined", "differ", "missing_on_dst", "missing_on_src", "match", "error", "exit_code")

# Sets CONFIG (the --config option, if a config file was given), exposes the
# staged Dir as the alias remote `staged:`, and defines `rclone_path PATH`,
# which turns a URL into an rclone path as upstream does, setting HTTP_URL to
# the --http-url option an http(s) URL needs.
REMOTES = r"""
        shopt -s nullglob; CONFIG_FILE=({inputs.config}); shopt -u nullglob
        CONFIG=()
        if [ ${#CONFIG_FILE[@]} -gt 0 ]; then CONFIG=(--config "${CONFIG_FILE[0]}"); fi
        if [ -d {inputs.staged} ]; then
            export RCLONE_CONFIG_STAGED_TYPE=alias RCLONE_CONFIG_STAGED_REMOTE={inputs.staged}
        fi
        HTTP_URL=()
        rclone_path() {
            local URL_PATTERN='^(https?://[^/]+)(/.*)?$'
            if [[ $1 =~ $URL_PATTERN ]]; then
                local URL_PATH=${BASH_REMATCH[2]:-/}
                HTTP_URL=(--http-url "${BASH_REMATCH[1]}")
                RCLONE_PATH=":http:${URL_PATH#/}"
            else
                RCLONE_PATH=$(echo "$1" | sed -E 's#^([a-zA-Z][a-zA-Z0-9+.-]*)://#\1:#')
            fi
        }
"""

# Run in the output Dir: write and sort the reports, drop the empty ones.
REPORT_FILES = r"""
        touch "$P.combined.txt" "$P.differ.txt" "$P.missing_on_dst.txt" "$P.missing_on_src.txt" \
            "$P.match.txt" "$P.error.txt"
        REPORT_OPTIONS=(--combined "$P.combined.txt" --differ "$P.differ.txt"
            --missing-on-dst "$P.missing_on_dst.txt" --missing-on-src "$P.missing_on_src.txt"
            --match "$P.match.txt" --error "$P.error.txt" --checkers "$TASK_CPUS")
"""

SORT_REPORTS = r"""
        export LC_ALL=C
        sort -k2 "$P.combined.txt" -o "$P.combined.txt"
        for REPORT in differ missing_on_dst missing_on_src match error; do
            sort "$P.$REPORT.txt" -o "$P.$REPORT.txt"
        done
        for f in *.txt; do [ -s "$f" ] || rm -f "$f"; done
"""

CHECK_INPUTS = {"prefix": str, "config": File | None, "staged": Dir | None, "args": str}

rclone_check_cmd = shell.create(
    name="rclone_check",
    image=RCLONE_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"source": str, "destination": str, **CHECK_INPUTS},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        {REMOTES}
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        rclone_path {{inputs.source}}; SOURCE=$RCLONE_PATH
        cd {{outputs.results}}
        {REPORT_FILES}
        rclone check "${{CONFIG[@]}}" "${{HTTP_URL[@]}}" $ARGS "${{REPORT_OPTIONS[@]}}" \
            "$SOURCE" {{inputs.destination}} \
            && echo 0 > "$P.exit_code.txt" || echo $? > "$P.exit_code.txt"
        {SORT_REPORTS}
    """,
)

# As upstream, symlinks on the destination are followed (--copy-links).
rclone_checksum_cmd = shell.create(
    name="rclone_checksum",
    image=RCLONE_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"sumfile": File, "hash": str, "destination": str, **CHECK_INPUTS},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        {REMOTES}
        SUMFILE=({{inputs.sumfile}})
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        rclone_path {{inputs.destination}}; DESTINATION=$RCLONE_PATH
        cd {{outputs.results}}
        {REPORT_FILES}
        rclone checksum "${{CONFIG[@]}}" --copy-links "${{HTTP_URL[@]}}" $ARGS "${{REPORT_OPTIONS[@]}}" \
            {{inputs.hash}} "${{SUMFILE[0]}}" "$DESTINATION" \
            && echo 0 > "$P.exit_code.txt" || echo $? > "$P.exit_code.txt"
        {SORT_REPORTS}
    """,
)

# Runs in a working directory holding the filter file under its own name, so
# that args can name it as upstream's do (e.g. `--exclude-from exclude.txt`).
rclone_copy_cmd = shell.create(
    name="rclone_copy",
    image=RCLONE_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"source": str, "destination": str, "filter_file": File | None, **CHECK_INPUTS},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        {REMOTES}
        shopt -s nullglob; FILTER=({{inputs.filter_file}}); shopt -u nullglob
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        rclone_path {{inputs.source}}; SOURCE=$RCLONE_PATH
        cd "$(mktemp -d)"
        if [ ${{#FILTER[@]}} -gt 0 ]; then ln -s "${{FILTER[0]}}" "$(basename "${{FILTER[0]}}")"; fi
        TRANSFERS=$(( TASK_CPUS / 2 )); if [ $TRANSFERS -lt 1 ]; then TRANSFERS=1; fi
        rclone "${{CONFIG[@]}}" copy "${{HTTP_URL[@]}}" $ARGS --log-file "$OUT/$P-rclone-copy.log" \
            --transfers $TRANSFERS --checkers $TASK_CPUS "$SOURCE" {{inputs.destination}} || exit 1
    """,
)

env = flyte.TaskEnvironment.from_task(
    "rclone",
    rclone_check_cmd.as_task(),
    rclone_checksum_cmd.as_task(),
    rclone_copy_cmd.as_task(),
)


@dataclass
class RcloneCheckResult:
    """The reports that came out non-empty; ``exit_code`` is always there."""

    combined: File | None  # <prefix>.combined.txt
    differ: File | None  # <prefix>.differ.txt
    missing_on_dst: File | None  # <prefix>.missing_on_dst.txt
    missing_on_src: File | None  # <prefix>.missing_on_src.txt
    match: File | None  # <prefix>.match.txt
    error: File | None  # <prefix>.error.txt
    exit_code: File  # <prefix>.exit_code.txt
    results: Dir  # everything the task wrote


async def check_result(results: Dir, prefix: str, tool: str) -> RcloneCheckResult:
    reports = {name: await results.get_file(f"{prefix}.{name}.txt") for name in REPORTS}
    if reports["exit_code"] is None:
        raise FileNotFoundError(f"rclone {tool} wrote no {prefix}.exit_code.txt")
    return RcloneCheckResult(**reports, results=results)


async def rclone_check(
    source: str,
    destination: str,
    prefix: str,
    config: File | None = None,
    staged: Dir | None = None,
    args: str = "",
) -> RcloneCheckResult:
    """Compare ``source`` with ``destination`` (rclone paths or URLs); ``args`` go to ``rclone check``."""
    results = await rclone_check_cmd(
        source=source, destination=destination, prefix=prefix, config=config, staged=staged, args=args
    )
    return await check_result(results, prefix, "check")


async def rclone_checksum(
    sumfile: File,
    hash: str,
    destination: str,
    prefix: str,
    config: File | None = None,
    staged: Dir | None = None,
    args: str = "",
) -> RcloneCheckResult:
    """Check ``destination`` against ``sumfile``, a list of ``hash`` sums (e.g. ``MD5``) and paths."""
    results = await rclone_checksum_cmd(
        sumfile=sumfile, hash=hash, destination=destination, prefix=prefix, config=config, staged=staged, args=args
    )
    return await check_result(results, prefix, "checksum")


async def rclone_copy(
    source: str,
    destination: str,
    prefix: str,
    filter_file: File | None = None,
    config: File | None = None,
    staged: Dir | None = None,
    args: str = "",
) -> File:
    """Copy ``source`` to ``destination``; returns rclone's log, ``<prefix>-rclone-copy.log``.

    ``args`` go to ``rclone copy`` and may name ``filter_file`` by its file
    name, e.g. ``--exclude-from exclude.txt``.
    """
    results = await rclone_copy_cmd(
        source=source,
        destination=destination,
        filter_file=filter_file,
        prefix=prefix,
        config=config,
        staged=staged,
        args=args,
    )
    log = await results.get_file(f"{prefix}-rclone-copy.log")
    if log is None:
        raise FileNotFoundError(f"rclone copy wrote no {prefix}-rclone-copy.log")
    return log
