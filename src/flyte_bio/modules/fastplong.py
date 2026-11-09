"""fastplong — quality control, filtering and trimming of long reads.

Exposes :func:`fastplong`, which processes one sample's long reads (PacBio or
Nanopore FASTQ) and returns the filtered reads and fastplong's JSON, HTML and
log reports. As upstream, the outputs are:

- ``<prefix>.fastplong.fastq.gz`` (none with ``discard_trimmed_pass``);
- ``<prefix>.fastplong.{json,html,log}``;
- with ``save_trimmed_fail``: ``<prefix>.fail.fastq.gz``.

``--thread`` follows the task's CPUs (the JSON report records it). Unless ``args`` sets
``--report_title``, the report is titled ``<prefix>_fastplong_report``.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

FASTPLONG_IMAGE = "community.wave.seqera.io/library/fastplong:0.4.1--73e8274104613b58"
DEFAULT_RESOURCES = flyte.Resources(cpu=4, memory="8Gi")

fastplong_cmd = shell.create(
    name="fastplong",
    image=FASTPLONG_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "reads": File,
        "adapter_fasta": File | None,
        "prefix": str,
        "discard_trimmed_pass": bool,
        "save_trimmed_fail": bool,
        "args": str,
    },
    defaults={"discard_trimmed_pass": False, "save_trimmed_fail": False, "args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        READS=({{inputs.reads}})
        shopt -s nullglob; ADAPTER=({{inputs.adapter_fasta}}); shopt -u nullglob
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        # As upstream, the reads are read under their own name from the working directory: the JSON report
        # records the command line.
        NAME=$(basename "${{READS[0]}}")
        ln -s "${{READS[0]}}" "$NAME"
        if [[ " $ARGS " != *" --report_title "* ]]; then ARGS="$ARGS --report_title ${{P}}_fastplong_report"; fi
        OUTS=(); FAIL=(); ADAPTERS=()
        if [ {{inputs.discard_trimmed_pass}} != true ]; then OUTS=(--out "$P.fastplong.fastq.gz"); fi
        if [ ${{#ADAPTER[@]}} -gt 0 ]; then
            ln -s "${{ADAPTER[0]}}" "$(basename "${{ADAPTER[0]}}")"
            ADAPTERS=(--adapter_fasta "$(basename "${{ADAPTER[0]}}")")
        fi
        if [ {{inputs.save_trimmed_fail}} = true ]; then FAIL=(--failed_out "$P.fail.fastq.gz"); fi
        fastplong --in "$NAME" "${{OUTS[@]}}" --json "$P.fastplong.json" --html "$P.fastplong.html" \
            "${{ADAPTERS[@]}}" "${{FAIL[@]}}" --thread $TASK_CPUS $ARGS 2> "$P.fastplong.log" || {{
            tail -n 20 "$P.fastplong.log" >&2; exit 1;
        }}
        cat "$P.fastplong.log" >&2
        shopt -s nullglob
        mv ./*.fastplong.fastq.gz ./*.json ./*.html ./*.log ./*.fail.fastq.gz "$OUT"/
    """,
)

env = flyte.TaskEnvironment.from_task("fastplong", fastplong_cmd.as_task())


@dataclass
class FastplongResult:
    reads: File | None  # <prefix>.fastplong.fastq.gz; None with discard_trimmed_pass
    json: File
    html: File
    log: File
    reads_fail: File | None  # <prefix>.fail.fastq.gz, with save_trimmed_fail
    results: Dir  # everything above


async def fastplong(
    reads: File,
    adapter_fasta: File | None = None,
    prefix: str = "sample",
    discard_trimmed_pass: bool = False,
    save_trimmed_fail: bool = False,
    args: str = "",
) -> FastplongResult:
    """Filter and trim one sample's long reads."""
    results = await fastplong_cmd(
        reads=reads,
        adapter_fasta=adapter_fasta,
        prefix=prefix,
        discard_trimmed_pass=discard_trimmed_pass,
        save_trimmed_fail=save_trimmed_fail,
        args=args,
    )

    return await read_fastplong_results(results, prefix, discard_trimmed_pass, save_trimmed_fail)


async def read_fastplong_results(
    results: Dir, prefix: str, discard_trimmed_pass: bool = False, save_trimmed_fail: bool = False
) -> FastplongResult:
    """Pick the outputs out of ``fastplong_cmd``'s results Dir (e.g. from an overridden task)."""

    async def need(name: str) -> File:
        found = await results.get_file(name)
        if found is None:
            raise FileNotFoundError(f"fastplong wrote no {name}")
        return found

    return FastplongResult(
        reads=None if discard_trimmed_pass else await need(f"{prefix}.fastplong.fastq.gz"),
        json=await need(f"{prefix}.fastplong.json"),
        html=await need(f"{prefix}.fastplong.html"),
        log=await need(f"{prefix}.fastplong.log"),
        reads_fail=await need(f"{prefix}.fail.fastq.gz") if save_trimmed_fail else None,
        results=results,
    )
