"""mgikit — MGI (DNBSEQ) FASTQ demultiplexing.

Exposes :func:`mgikit_demultiplex`, which runs ``mgikit demultiplex`` on an
MGI run folder (a lane's ``*_read_1.fq.gz`` / ``*_read_2.fq.gz`` and
``BioInfo.csv``), given as a ``.tar.gz`` archive, with a sample sheet
(``sample_id``, ``i7``, ``i5``, ``i7_rc``, ``i5_rc``, ``template``). The
archive is extracted into ``<archive name without .tar.gz>/``, dropping its
top-level directory when it has just one, as upstream's unzip step names it.

mgikit writes one FASTQ per sample and read into ``<prefix>/``, which can be
many, so they come back as one archive. The output Dir holds:

- ``fastq.tar.gz``: ``<sample>_S<n>_L0<l>_R1_001.fastq.gz``,
  ``..._R2_001.fastq.gz``;
- the undetermined and ambiguous reads, which upstream moves out of the
  sample folder: ``Undetermined_L0<l>_R1_001.fastq.gz``, ... and
  ``Ambiguous_...`` when mgikit writes them;
- the reports: ``<flowcell>.L0<l>.mgikit.general``, ``.info``,
  ``.sample_stats``, ``.undetermined_barcode`` (and ``.complete``), and
  ``.ambiguous_barcode`` when there are ambiguous reads.

mgikit's ``--threads`` is the task's CPUs unless ``args`` sets it (upstream's
test passes ``--threads 1``).
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

# Upstream's main.nf names the quay.io biocontainer directly.
MGIKIT_IMAGE = "quay.io/biocontainers/mgikit:2.1.0--h3ab6199_0"
# Upstream's process_medium.
DEFAULT_RESOURCES = flyte.Resources(cpu=6, memory="12Gi")

mgikit_demultiplex_cmd = shell.create(
    name="mgikit_demultiplex",
    image=MGIKIT_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"samplesheet": File, "run_dir": File, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        set -e
        SHEET=({{inputs.samplesheet}})
        RUN=({{inputs.run_dir}})
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        NAME=$(basename "${{RUN[0]}}" .tar.gz)
        mkdir -p "$NAME"
        # As upstream's untar: drop the top-level directory only when the archive has just one.
        if [[ $(tar -tzf "${{RUN[0]}}" | grep / | cut -d/ -f1 | sort -u | wc -l) -eq 1 ]]; then
            tar -C "$NAME" --strip-components 1 -xzf "${{RUN[0]}}"
        else
            tar -C "$NAME" -xzf "${{RUN[0]}}"
        fi
        THREADS=(--threads $TASK_CPUS)
        if [[ " $ARGS " == *" --threads "* ]]; then THREADS=(); fi
        # mgikit is a Rust tool: a panic exits 255, which the platform would retry forever.
        mgikit demultiplex -i "$NAME" -s "${{SHEET[0]}}" -o "$P" "${{THREADS[@]}}" $ARGS || exit 1
        shopt -s nullglob
        UNASSIGNED=("$P"/Undetermined*.fastq.gz "$P"/Ambiguous*.fastq.gz)
        if [ ${{#UNASSIGNED[@]}} -gt 0 ]; then mv "${{UNASSIGNED[@]}}" "$OUT"/; fi
        (cd "$P" && ls | grep '\.fastq\.gz$' | tar -czf "$OUT/fastq.tar.gz" -T -)
        cp "$P"/*mgikit.* "$OUT"/
    """,
)

env = flyte.TaskEnvironment.from_task("mgikit", mgikit_demultiplex_cmd.as_task())


@dataclass
class MgikitDemultiplexResult:
    fastq: File  # fastq.tar.gz: the samples' FASTQs
    undetermined: list[File]  # Undetermined_*.fastq.gz
    ambiguous: list[File]  # Ambiguous_*.fastq.gz
    general: File  # *.mgikit.general
    info: File  # *.mgikit.info
    sample_stats: File  # *.mgikit.sample_stats
    results: Dir  # all of the above, and the barcode reports


async def mgikit_demultiplex(
    samplesheet: File, run_dir: File, prefix: str = "out", args: str = ""
) -> MgikitDemultiplexResult:
    """Demultiplex the MGI run folder archive ``run_dir`` (``.tar.gz``) with ``samplesheet``."""
    results = await mgikit_demultiplex_cmd(samplesheet=samplesheet, run_dir=run_dir, prefix=prefix, args=args)
    files = sorted(await results.list_files(), key=lambda f: f.path)

    def named(test) -> list[File]:
        return [f for f in files if test(f.path.rsplit("/", 1)[-1])]

    def one(suffix: str) -> File:
        found = named(lambda name: name.endswith(suffix))
        if len(found) != 1:
            raise FileNotFoundError(f"mgikit wrote {len(found)} *{suffix} reports")
        return found[0]

    fastq = named(lambda name: name == "fastq.tar.gz")
    if not fastq:
        raise FileNotFoundError("mgikit wrote no FASTQs")
    return MgikitDemultiplexResult(
        fastq=fastq[0],
        undetermined=named(lambda name: name.startswith("Undetermined") and name.endswith(".fastq.gz")),
        ambiguous=named(lambda name: name.startswith("Ambiguous") and name.endswith(".fastq.gz")),
        general=one(".mgikit.general"),
        info=one(".mgikit.info"),
        sample_stats=one(".mgikit.sample_stats"),
        results=results,
    )
