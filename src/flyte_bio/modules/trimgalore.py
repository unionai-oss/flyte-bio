"""trimgalore — adapter and quality trimming of FASTQ reads.

Exposes :func:`trimgalore`, which trims one sample's reads (one or two mates)
and returns the trimmed reads, the trimming report(s) and, optionally,
FastQC reports on the trimmed reads.

As upstream, the reads are first symlinked to ``<prefix>.fastq.gz`` or
``<prefix>_1.fastq.gz`` / ``<prefix>_2.fastq.gz``, so outputs are named
after the sample: ``<prefix>_trimmed.fq.gz`` (single-end) or
``<prefix>_1_val_1.fq.gz`` / ``<prefix>_2_val_2.fq.gz`` (paired).
Everything Trim Galore writes lands in one output Dir; the wrapper picks
the files out by name.
"""

import re
from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

TRIMGALORE_IMAGE = "community.wave.seqera.io/library/trim-galore:2.1.0--27e6376b8f6c1872"

# Upstream runs Trim Galore with 8 CPUs; --cores is derived from it below.
TRIMGALORE_CPUS = 8
DEFAULT_RESOURCES = flyte.Resources(cpu=TRIMGALORE_CPUS, memory="8Gi")

# Upstream's --cores rule: CPUs minus Trim Galore's own overhead (4 paired,
# 3 single-end), clamped to 1..8. `--fastqc_args` is set here rather than via
# `args` because its value is a quoted multi-word string, which `args`
# (expanded unquoted, one word per flag value) can't carry.
trimgalore_cmd = shell.create(
    name="trimgalore",
    image=TRIMGALORE_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"reads_1": File, "reads_2": File | None, "prefix": str, "fastqc": bool, "args": str},
    defaults={"fastqc": False, "args": ""},
    outputs={"results": Dir},
    script=rf"""
        shopt -s nullglob; R2=({{inputs.reads_2}}); shopt -u nullglob
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        FQ=()
        if [ {{inputs.fastqc}} = true ]; then FQ=(--fastqc_args "-t {TRIMGALORE_CPUS}"); fi
        W=$(mktemp -d); cd "$W"
        if [ ${{#R2[@]}} -eq 0 ]; then
            ln -s {{inputs.reads_1}} "$P.fastq.gz"
            trim_galore "${{FQ[@]}}" $ARGS --cores {max(1, min(8, TRIMGALORE_CPUS - 3))} --gzip \
                -o {{outputs.results}} "$P.fastq.gz"
        else
            ln -s {{inputs.reads_1}} "${{P}}_1.fastq.gz"
            ln -s "${{R2[0]}}" "${{P}}_2.fastq.gz"
            trim_galore "${{FQ[@]}}" $ARGS --cores {max(1, min(8, TRIMGALORE_CPUS - 4))} --paired --gzip \
                -o {{outputs.results}} "${{P}}_1.fastq.gz" "${{P}}_2.fastq.gz"
        fi
    """,
)


env = flyte.TaskEnvironment.from_task(
    "trimgalore",
    trimgalore_cmd.as_task(),
)


@dataclass
class TrimGaloreResult:
    reads_1: File
    reads_2: File | None  # None for single-end
    reports: list[File]  # one trimming report per mate
    results: Dir  # everything Trim Galore wrote, incl. FastQC html/zip when requested

    async def reads_after_filtering(self) -> float:
        """Reads left after length filtering, parsed from the (last) trimming report as upstream does."""
        async with self.reports[-1].open("rb") as fh:
            text = bytes(await fh.read()).decode()
        total = filtered = 0.0
        for line in text.splitlines():
            if m := re.search(r"([\d\.]+)\ssequences processed in total", line):
                total = float(m.group(1))
            if m := re.search(r"shorter than the length cutoff[^:]+:\s*([\d\.]+)", line):
                filtered = float(m.group(1))
        return total - filtered


async def trimgalore(
    reads_1: File,
    reads_2: File | None = None,
    prefix: str = "trimmed",
    fastqc: bool = False,
    args: str = "",
) -> TrimGaloreResult:
    """Trim one sample's reads; ``reads_2`` is None for single-end."""
    results = await trimgalore_cmd(
        reads_1=reads_1,
        reads_2=reads_2,
        prefix=prefix,
        fastqc=fastqc,
        args=args,
    )

    async def pick(name: str) -> File:
        found = await results.get_file(name)
        if found is None:
            raise FileNotFoundError(f"trimgalore wrote no {name}")
        return found

    if reads_2 is None:
        trimmed = await pick(f"{prefix}_trimmed.fq.gz")
        report = await pick(f"{prefix}.fastq.gz_trimming_report.txt")
        return TrimGaloreResult(reads_1=trimmed, reads_2=None, reports=[report], results=results)
    val_1, val_2, report_1, report_2 = (
        await pick(f"{prefix}_1_val_1.fq.gz"),
        await pick(f"{prefix}_2_val_2.fq.gz"),
        await pick(f"{prefix}_1.fastq.gz_trimming_report.txt"),
        await pick(f"{prefix}_2.fastq.gz_trimming_report.txt"),
    )
    return TrimGaloreResult(reads_1=val_1, reads_2=val_2, reports=[report_1, report_2], results=results)
