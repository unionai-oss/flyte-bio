"""modkit — base-modification calls and pileups from modBAMs (``MM``/``ML`` tags).

Exposes:

- :func:`modkit_callmods` — ``modkit call-mods``: turn each read's
  modification probabilities into calls, writing ``<prefix>.bam`` and
  modkit's log ``<prefix>.log``.
- :func:`modkit_pileup` — ``modkit pileup``: per-reference-position
  modification counts as bgzipped bedMethyl, ``<prefix>.bed.gz``. With
  ``--phased`` in ``args`` modkit writes one file per haplotype and a
  combined one (``<prefix>_hp1.bed.gz``, ``<prefix>_hp2.bed.gz``,
  ``<prefix>_combined.bed.gz``, ...). ``fasta`` (with its ``fai``) is passed
  as ``--ref`` and ``bed`` restricts the pileup (``--include-bed``).
- :func:`modkit_repair` — ``modkit repair``: copy the modification tags of
  the untrimmed reads (``donor_bam``) onto their trimmed copies
  (``acceptor_bam``), writing ``<prefix>.bam`` and ``<prefix>.log``.

As upstream, every tool gets the task's CPUs as threads (pileup also as
bgzip threads) and the inputs are linked into a working directory under their
own names, beside their indexes. ``args`` go to the modkit subcommand.
"""

import os
from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

MODKIT_IMAGE = "community.wave.seqera.io/library/ont-modkit:0.6.4--0c623277f8775a30"
# Upstream labels call-mods and repair process_medium and pileup process_high;
# its tests ran every tool on 4 CPUs.
DEFAULT_RESOURCES = flyte.Resources(cpu=4, memory="8Gi")
PILEUP_RESOURCES = flyte.Resources(cpu=4, memory="16Gi")

modkit_callmods_cmd = shell.create(
    name="modkit_callmods",
    image=MODKIT_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"bam": File, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        BAM=({{inputs.bam}})
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        NAME=$(basename "${{BAM[0]}}"); ln -s "${{BAM[0]}}" "$NAME"
        modkit call-mods $ARGS "$NAME" "$P.bam" --threads $TASK_CPUS --log-filepath "./$P.log" || exit 1
        mv "$P.bam" "$P.log" "$OUT"/
    """,
)

modkit_pileup_cmd = shell.create(
    name="modkit_pileup",
    image=MODKIT_IMAGE,
    resources=PILEUP_RESOURCES,
    inputs={
        "bam": File,
        "bai": File,
        "fasta": File | None,
        "fai": File | None,
        "bed": File | None,
        "prefix": str,
        "args": str,
    },
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        BAM=({{inputs.bam}}); BAI=({{inputs.bai}})
        shopt -s nullglob; FA=({{inputs.fasta}}); FAI=({{inputs.fai}}); BED=({{inputs.bed}}); shopt -u nullglob
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"

        NAME=$(basename "${{BAM[0]}}"); ln -s "${{BAM[0]}}" "$NAME"
        case "${{BAI[0]}}" in *.csi) ln -s "${{BAI[0]}}" "$NAME.csi" ;; *) ln -s "${{BAI[0]}}" "$NAME.bai" ;; esac
        REF=()
        if [ ${{#FA[@]}} -gt 0 ]; then
            FA_NAME=$(basename "${{FA[0]}}"); ln -s "${{FA[0]}}" "$FA_NAME"
            if [ ${{#FAI[@]}} -gt 0 ]; then ln -s "${{FAI[0]}}" "$FA_NAME.fai"; fi
            REF=(--ref "$FA_NAME")
        fi
        INCLUDE=()
        if [ ${{#BED[@]}} -gt 0 ]; then
            BED_NAME=$(basename "${{BED[0]}}"); ln -s "${{BED[0]}}" "$BED_NAME"
            INCLUDE=(--include-bed "$BED_NAME")
        fi

        modkit pileup $ARGS \
            --bgzf --bgzf-threads $TASK_CPUS --threads $TASK_CPUS \
            --prefix "$P" "${{REF[@]}}" "${{INCLUDE[@]}}" "$NAME" "$P.tmp" || exit 1

        # A phased pileup writes a directory of files; an unphased one, a single file.
        if [ -d "$P.tmp" ]; then
            for f in "$P.tmp"/*; do
                if [ -f "$f" ]; then mv "$f" "$OUT"/; fi
            done
        else
            mv "$P.tmp" "$OUT/$P.bed.gz"
        fi
        for f in ./*.log; do
            if [ -f "$f" ]; then mv "$f" "$OUT"/; fi
        done
    """,
)

modkit_repair_cmd = shell.create(
    name="modkit_repair",
    image=MODKIT_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"donor_bam": File, "acceptor_bam": File, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        DONOR=({{inputs.donor_bam}}); ACCEPTOR=({{inputs.acceptor_bam}})
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        DONOR_NAME=$(basename "${{DONOR[0]}}"); ln -s "${{DONOR[0]}}" "$DONOR_NAME"
        ACCEPTOR_NAME=$(basename "${{ACCEPTOR[0]}}"); ln -s "${{ACCEPTOR[0]}}" "$ACCEPTOR_NAME"
        modkit repair $ARGS \
            --threads $TASK_CPUS \
            --donor-bam "$DONOR_NAME" \
            --acceptor-bam "$ACCEPTOR_NAME" \
            --output-bam "$P.bam" \
            --log-filepath "./$P.log" || exit 1
        mv "$P.bam" "$P.log" "$OUT"/
    """,
)

env = flyte.TaskEnvironment.from_task(
    "modkit", modkit_callmods_cmd.as_task(), modkit_pileup_cmd.as_task(), modkit_repair_cmd.as_task()
)


@dataclass
class ModkitBamResult:
    bam: File  # <prefix>.bam
    log: File  # <prefix>.log
    results: Dir  # both


@dataclass
class ModkitPileupResult:
    beds: list[File]  # <prefix>.bed.gz, or with --phased one file per haplotype plus combined, by name
    log: File | None  # a log, when args asked for one (--log-filepath)
    results: Dir  # everything modkit wrote


def file_name(file: File) -> str:
    return file.name or os.path.basename(file.path)


def check_prefix(prefix: str, *bams: File) -> None:
    for bam in bams:
        if file_name(bam) == f"{prefix}.bam":
            raise ValueError(f"input and output are both {prefix}.bam: set another prefix")


async def bam_and_log(results: Dir, prefix: str, tool: str) -> ModkitBamResult:
    bam = await results.get_file(f"{prefix}.bam")
    log = await results.get_file(f"{prefix}.log")
    if bam is None or log is None:
        raise FileNotFoundError(f"modkit {tool} wrote no {prefix}.bam or {prefix}.log")
    return ModkitBamResult(bam=bam, log=log, results=results)


async def modkit_callmods(bam: File, prefix: str = "sample", args: str = "") -> ModkitBamResult:
    """Modification calls for the reads in ``bam``: ``<prefix>.bam`` and ``<prefix>.log``."""
    check_prefix(prefix, bam)
    results = await modkit_callmods_cmd(bam=bam, prefix=prefix, args=args)
    return await bam_and_log(results, prefix, "call-mods")


async def modkit_pileup(
    bam: File,
    bai: File,
    fasta: File | None = None,
    fai: File | None = None,
    bed: File | None = None,
    prefix: str = "sample",
    args: str = "",
) -> ModkitPileupResult:
    """bedMethyl pileup of ``bam`` (with its index ``bai``), optionally against ``fasta`` and within ``bed``."""
    results = await modkit_pileup_cmd(bam=bam, bai=bai, fasta=fasta, fai=fai, bed=bed, prefix=prefix, args=args)
    files = sorted(await results.list_files(), key=file_name)
    beds = [f for f in files if file_name(f).endswith(".bed.gz")]
    if not beds:
        raise FileNotFoundError("modkit pileup wrote no .bed.gz")
    logs = [f for f in files if file_name(f).endswith(".log")]
    return ModkitPileupResult(beds=beds, log=logs[0] if logs else None, results=results)


async def modkit_repair(donor_bam: File, acceptor_bam: File, prefix: str = "sample", args: str = "") -> ModkitBamResult:
    """``acceptor_bam``'s reads with ``donor_bam``'s modification tags: ``<prefix>.bam`` and ``<prefix>.log``."""
    check_prefix(prefix, donor_bam, acceptor_bam)
    results = await modkit_repair_cmd(donor_bam=donor_bam, acceptor_bam=acceptor_bam, prefix=prefix, args=args)
    return await bam_and_log(results, prefix, "repair")
