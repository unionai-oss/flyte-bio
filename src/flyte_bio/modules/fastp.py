"""fastp — adapter and quality trimming of FASTQ reads, with QC reports.

Exposes :func:`fastp`, which trims one sample's reads (single-end, paired-end,
or interleaved with ``--interleaved_in`` in ``args``) and returns the trimmed
reads and fastp's JSON, HTML and log reports.

As upstream, the reads are first symlinked to ``<prefix>.fastq.gz`` or
``<prefix>_R1.fastq.gz`` / ``<prefix>_R2.fastq.gz``, so the reports name the
sample consistently. Outputs keep upstream's names:

- trimmed reads: ``<prefix>.fastp.fastq.gz`` (single-end, interleaved) or
  ``<prefix>_R1.fastp.fastq.gz`` / ``<prefix>_R2.fastp.fastq.gz``;
- reports: ``<prefix>.fastp.{json,html,log}``;
- with ``save_trimmed_fail``: ``<prefix>.fail.fastq.gz`` (single-end), or
  ``<prefix>.paired.fail.fastq.gz``, ``<prefix>_R1.fail.fastq.gz`` and
  ``<prefix>_R2.fail.fastq.gz``;
- with ``save_merged`` (paired-end): ``<prefix>.merged.fastq.gz``.

``discard_trimmed_pass`` runs fastp for its reports only, writing no reads.
fastp's ``--thread`` excludes its reader threads, so as upstream it is the
task's CPUs minus 1 (single-end) or 2 (paired or interleaved), at least 2.
Paired-end runs add ``--detect_adapter_for_pe``, as upstream.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

FASTP_IMAGE = "community.wave.seqera.io/library/fastp:1.4.0--df415eab5e33cb11"
DEFAULT_RESOURCES = flyte.Resources(cpu=4, memory="8Gi")

fastp_cmd = shell.create(
    name="fastp",
    image=FASTP_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "reads_1": File,
        "reads_2": File | None,
        "adapter_fasta": File | None,
        "prefix": str,
        "discard_trimmed_pass": bool,
        "save_trimmed_fail": bool,
        "save_merged": bool,
        "args": str,
    },
    defaults={"discard_trimmed_pass": False, "save_trimmed_fail": False, "save_merged": False, "args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        R1=({{inputs.reads_1}})
        shopt -s nullglob; R2=({{inputs.reads_2}}); ADAPTER=({{inputs.adapter_fasta}}); shopt -u nullglob
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"

        INTERLEAVED=false
        if [[ " $ARGS " == *" --interleaved_in "* ]]; then INTERLEAVED=true; fi
        SINGLE=true
        if [ ${{#R2[@]}} -gt 0 ]; then SINGLE=false; fi
        READERS=1
        if [ $INTERLEAVED = true ] || [ $SINGLE = false ]; then READERS=2; fi
        THREADS=$(( TASK_CPUS - READERS ))
        if [ $THREADS -lt 2 ]; then THREADS=2; fi

        EXTRA=()
        if [ ${{#ADAPTER[@]}} -gt 0 ]; then EXTRA+=(--adapter_fasta "${{ADAPTER[0]}}"); fi
        if [ {{inputs.save_trimmed_fail}} = true ]; then
            if [ $SINGLE = true ]; then
                EXTRA+=(--failed_out "$P.fail.fastq.gz")
            else
                EXTRA+=(--failed_out "$P.paired.fail.fastq.gz")
                EXTRA+=(--unpaired1 "${{P}}_R1.fail.fastq.gz" --unpaired2 "${{P}}_R2.fail.fastq.gz")
            fi
        fi
        DISCARD={{inputs.discard_trimmed_pass}}

        if [ $INTERLEAVED = true ]; then
            ln -s "${{R1[0]}}" "$P.fastq.gz"
            fastp --stdout --in1 "$P.fastq.gz" --thread $THREADS \
                --json "$P.fastp.json" --html "$P.fastp.html" "${{EXTRA[@]}}" $ARGS \
                2> "$P.fastp.log" | gzip -c > "$P.fastp.fastq.gz"
        elif [ $SINGLE = true ]; then
            ln -s "${{R1[0]}}" "$P.fastq.gz"
            OUTS=(--out1 "$P.fastp.fastq.gz")
            if [ $DISCARD = true ]; then OUTS=(); fi
            fastp --in1 "$P.fastq.gz" "${{OUTS[@]}}" --thread $THREADS \
                --json "$P.fastp.json" --html "$P.fastp.html" "${{EXTRA[@]}}" $ARGS \
                2> "$P.fastp.log"
        else
            ln -s "${{R1[0]}}" "${{P}}_R1.fastq.gz"
            ln -s "${{R2[0]}}" "${{P}}_R2.fastq.gz"
            OUTS=(--out1 "${{P}}_R1.fastp.fastq.gz" --out2 "${{P}}_R2.fastp.fastq.gz")
            if [ $DISCARD = true ]; then OUTS=(); fi
            if [ {{inputs.save_merged}} = true ]; then EXTRA+=(-m --merged_out "$P.merged.fastq.gz"); fi
            fastp --in1 "${{P}}_R1.fastq.gz" --in2 "${{P}}_R2.fastq.gz" "${{OUTS[@]}}" \
                --json "$P.fastp.json" --html "$P.fastp.html" "${{EXTRA[@]}}" \
                --thread $THREADS --detect_adapter_for_pe $ARGS \
                2> "$P.fastp.log"
        fi
        cat "$P.fastp.log" >&2
        shopt -s nullglob
        mv ./*.fastp.fastq.gz ./*.json ./*.html ./*.log ./*.fail.fastq.gz ./*.merged.fastq.gz "$OUT"/
    """,
)

env = flyte.TaskEnvironment.from_task("fastp", fastp_cmd.as_task())


@dataclass
class FastpResult:
    reads: list[File]  # the trimmed reads: one file, or R1 and R2; empty with discard_trimmed_pass
    json: File
    html: File
    log: File
    reads_fail: list[File]  # with save_trimmed_fail
    reads_merged: File | None  # with save_merged (paired-end)
    results: Dir  # everything fastp wrote


async def fastp(
    reads_1: File,
    reads_2: File | None = None,
    adapter_fasta: File | None = None,
    prefix: str = "sample",
    discard_trimmed_pass: bool = False,
    save_trimmed_fail: bool = False,
    save_merged: bool = False,
    args: str = "",
) -> FastpResult:
    """Trim one sample's reads; ``reads_2`` is None for single-end and interleaved reads."""
    results = await fastp_cmd(
        reads_1=reads_1,
        reads_2=reads_2,
        adapter_fasta=adapter_fasta,
        prefix=prefix,
        discard_trimmed_pass=discard_trimmed_pass,
        save_trimmed_fail=save_trimmed_fail,
        save_merged=save_merged,
        args=args,
    )

    async def pick(name: str) -> File:
        found = await results.get_file(name)
        if found is None:
            raise FileNotFoundError(f"fastp wrote no {name}")
        return found

    paired = reads_2 is not None
    reads = []
    if not discard_trimmed_pass:
        mates = ["_R1", "_R2"] if paired else [""]
        names = [f"{prefix}{mate}.fastp.fastq.gz" for mate in mates]
        reads = [await pick(n) for n in names]
    reads_fail = []
    if save_trimmed_fail:
        names = (
            [f"{prefix}.paired.fail.fastq.gz", f"{prefix}_R1.fail.fastq.gz", f"{prefix}_R2.fail.fastq.gz"]
            if paired
            else [f"{prefix}.fail.fastq.gz"]
        )
        reads_fail = [await pick(n) for n in names]
    return FastpResult(
        reads=reads,
        json=await pick(f"{prefix}.fastp.json"),
        html=await pick(f"{prefix}.fastp.html"),
        log=await pick(f"{prefix}.fastp.log"),
        reads_fail=reads_fail,
        reads_merged=await pick(f"{prefix}.merged.fastq.gz") if paired and save_merged else None,
        results=results,
    )
