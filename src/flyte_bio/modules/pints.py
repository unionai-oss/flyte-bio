"""pints — call transcriptional regulatory elements (TREs) from nascent-RNA reads with PINTS.

Exposes :func:`pints_caller`, which calls peaks on one sample's BAMs (each BAM
a replicate, numbered in the caller's order) for an assay such as ``GROseq``,
``PROseq`` or ``CoPRO``. Outputs keep PINTS's names, one per replicate ``<n>``:

- ``<prefix>_<n>_divergent_peaks.bed``;
- ``<prefix>_<n>_bidirectional_peaks.bed``;
- ``<prefix>_<n>_unidirectional_peaks.bed``;
- the log ``peakcalling_<timestamp>.log``.

As upstream, PINTS fails (exit 1) when the reads carry no signal. Its
``--thread`` follows the task's CPUs. BigWig inputs aren't supported, as
upstream.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

PINTS_IMAGE = "community.wave.seqera.io/library/pybedtools_bedtools_htslib_pip_pypints:39699b96998ec5f6"
# Upstream's process_medium asks for more; raise with pints_caller_cmd.as_task().override(resources=...).
DEFAULT_RESOURCES = flyte.Resources(cpu=2, memory="8Gi")

# A list[File] is staged under the files' names and globbed alphabetically, so
# the BAM names ("/"-separated) carry the caller's order, which numbers the replicates.
pints_caller_cmd = shell.create(
    name="pints_caller",
    image=PINTS_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "bams": list[File],
        "bam_names": str,
        "bais": list[File],
        "assay_type": str,
        "prefix": str,
        "args": str,
    },
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        BAMS=({{inputs.bams}})
        shopt -s nullglob; BAIS=({{inputs.bais}}); shopt -u nullglob
        BAM_NAMES={{inputs.bam_names}}
        ASSAY={{inputs.assay_type}}
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        mkdir -p .tmp; export TMPDIR="$W/.tmp"

        IFS=/ read -ra NAMES <<< "$BAM_NAMES"
        BAM_DIR=$(dirname "${{BAMS[0]}}")
        for NAME in "${{NAMES[@]}}"; do ln -s "$BAM_DIR/$NAME" "$NAME"; done
        for F in "${{BAIS[@]}}"; do ln -s "$F" "$(basename "$F")"; done

        pints_caller \
            --bam-file "${{NAMES[@]}}" \
            --save-to . \
            --file-prefix "$P" \
            --thread $TASK_CPUS \
            --dont-check-updates \
            --exp-type "$ASSAY" \
            $ARGS || exit 1
        shopt -s nullglob
        mv ./*_divergent_peaks.bed ./*_bidirectional_peaks.bed ./*_unidirectional_peaks.bed ./peakcalling_*.log "$OUT"/
    """,
)

env = flyte.TaskEnvironment.from_task("pints", pints_caller_cmd.as_task())


def file_name(file: File) -> str:
    return file.path.rstrip("/").rsplit("/", 1)[-1]


@dataclass
class PintsCallerResult:
    divergent: list[File]  # <prefix>_<n>_divergent_peaks.bed, by replicate
    bidirectional: list[File]  # <prefix>_<n>_bidirectional_peaks.bed
    unidirectional: list[File]  # <prefix>_<n>_unidirectional_peaks.bed
    log: File  # peakcalling_<timestamp>.log
    results: Dir  # everything above


async def pints_caller(
    bams: list[File],
    assay_type: str,
    bais: list[File] | None = None,
    prefix: str = "sample",
    args: str = "",
) -> PintsCallerResult:
    """Call TREs on ``bams`` (one per replicate) of assay ``assay_type`` (e.g. ``GROseq``, ``PROseq``).

    ``bais`` are optional indexes, linked beside the BAMs under their own names.
    """
    if not bams:
        raise ValueError("pints_caller needs BAMs")
    names = [file_name(f) for f in bams]
    if len(set(names)) != len(names):
        raise ValueError(f"pints_caller needs BAMs with distinct names, got {names}")
    results = await pints_caller_cmd(
        bams=bams,
        bam_names="/".join(names),
        bais=bais or [],
        assay_type=assay_type,
        prefix=prefix,
        args=args,
    )
    files = {file_name(f): f for f in await results.list_files()}

    def replicates(suffix: str) -> list[File]:
        found = [(name, f) for name, f in files.items() if name.endswith(suffix)]
        # <prefix>_<n><suffix>: order by the replicate number n.
        return [f for _, f in sorted(found, key=lambda item: int(item[0][: -len(suffix)].rsplit("_", 1)[-1]))]

    logs = [f for name, f in sorted(files.items()) if name.startswith("peakcalling_") and name.endswith(".log")]
    if not logs:
        raise FileNotFoundError(f"PINTS wrote no peakcalling log, only {sorted(files)}")
    return PintsCallerResult(
        divergent=replicates("_divergent_peaks.bed"),
        bidirectional=replicates("_bidirectional_peaks.bed"),
        unidirectional=replicates("_unidirectional_peaks.bed"),
        log=logs[0],
        results=results,
    )
