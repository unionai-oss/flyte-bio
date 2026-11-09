"""rastair — methylation calls from TAPS (C-to-T converted) alignments with rastair.

Exposes:

- :func:`rastair_mbias` — conversion bias along the reads (``rastair mbias``):
  ``<prefix>.rastair_mbias.txt``.
- :func:`rastair_mbiasparser` — plot that bias and pick how many bases to
  trim (rastair's ``plot_mbias.R`` and ``parse_mbias.R``):
  ``<prefix>.rastair_mbias_processed.csv``, its plot
  ``<prefix>.rastair_mbias_processed.pdf``, and the trimming for the original
  top and bottom strands (the CSV's first two lines, e.g. ``0,0,0,0``).
- :func:`rastair_call` — per-position methylation calls (``rastair call``),
  trimming reads by ``trim_ot`` / ``trim_ob`` (``--nOT`` / ``--nOB``):
  ``<prefix>.rastair_call.txt``.
- :func:`rastair_methylkit` — those calls in methylKit's format
  (``rastair_call_to_methylkit.sh``): ``<prefix>.rastair_methylkit.txt.gz``.

As upstream, ``rastair call`` and ``rastair mbias`` get the task's CPUs as
``--threads``; the BAM and its index, and the FASTA and its ``.fai``, are
linked side by side under their own names. ``args`` go to rastair.

rastair 0.8.2 is GPL-3.0 licensed (bioconda's record; its repository's
LICENSE is the AGPL-3.0). Later releases (2.x) moved to a non-commercial
licence.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

RASTAIR_IMAGE = "community.wave.seqera.io/library/rastair:0.8.2--bf70eeab4121509c"

# Upstream labels call and mbias process_medium, the parsers process_low. The
# 4 CPUs match upstream's tests.
MEDIUM_RESOURCES = flyte.Resources(cpu=4, memory="8Gi")
LOW_RESOURCES = flyte.Resources(cpu=1, memory="4Gi")

# Links the BAM, its index, the FASTA and its .fai into the working directory
# under their own names, and sets BAM_NAME and FASTA_NAME.
LINK_INPUTS = r"""
        W=$(mktemp -d); cd "$W"
        for F in "${BAM[0]}" "${BAI[0]}" "${FA[0]}" "${FAI[0]}"; do ln -s "$F" "$(basename "$F")"; done
        BAM_NAME=$(basename "${BAM[0]}")
        FASTA_NAME=$(basename "${FA[0]}")
"""

ALIGNMENT_INPUTS = {"bam": File, "bai": File, "fasta": File, "fai": File, "prefix": str}

rastair_call_cmd = shell.create(
    name="rastair_call",
    image=RASTAIR_IMAGE,
    resources=MEDIUM_RESOURCES,
    inputs={**ALIGNMENT_INPUTS, "trim_ot": str, "trim_ob": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        BAM=({{inputs.bam}}); BAI=({{inputs.bai}}); FA=({{inputs.fasta}}); FAI=({{inputs.fai}})
        P={{inputs.prefix}}
        TRIM_OT={{inputs.trim_ot}}
        TRIM_OB={{inputs.trim_ob}}
        ARGS={{inputs.args}}
        {LINK_INPUTS}
        rastair call $ARGS --threads $TASK_CPUS --nOT "$TRIM_OT" --nOB "$TRIM_OB" --fasta-file "$FASTA_NAME" \
            "$BAM_NAME" > {{outputs.results}}/"$P.rastair_call.txt" || exit 1
    """,
)

rastair_mbias_cmd = shell.create(
    name="rastair_mbias",
    image=RASTAIR_IMAGE,
    resources=MEDIUM_RESOURCES,
    inputs={**ALIGNMENT_INPUTS, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        BAM=({{inputs.bam}}); BAI=({{inputs.bai}}); FA=({{inputs.fasta}}); FAI=({{inputs.fai}})
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        {LINK_INPUTS}
        rastair mbias $ARGS --threads $TASK_CPUS --fasta-file "$FASTA_NAME" \
            "$BAM_NAME" > {{outputs.results}}/"$P.rastair_mbias.txt" || exit 1
    """,
)

rastair_mbiasparser_cmd = shell.create(
    name="rastair_mbiasparser",
    image=RASTAIR_IMAGE,
    resources=LOW_RESOURCES,
    inputs={"mbias": File, "prefix": str},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        MBIAS=({{inputs.mbias}})
        P={{inputs.prefix}}
        W=$(mktemp -d); cd "$W"
        NAME=$(basename "${{MBIAS[0]}}")
        ln -s "${{MBIAS[0]}}" "$NAME"
        plot_mbias.R --pdf -o "$P.rastair_mbias_processed.pdf" "$NAME" > "$P.rastair_mbias_processed.txt" || exit 1
        parse_mbias.R "$P.rastair_mbias_processed.txt" "$P.rastair_mbias_processed.csv" || exit 1
        shopt -s nullglob
        mv "$P.rastair_mbias_processed.csv" ./*.rastair_mbias_processed.pdf {{outputs.results}}/
    """,
)

rastair_methylkit_cmd = shell.create(
    name="rastair_methylkit",
    image=RASTAIR_IMAGE,
    resources=LOW_RESOURCES,
    inputs={"calls": File, "prefix": str},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        set -o pipefail
        CALLS=({{inputs.calls}})
        P={{inputs.prefix}}
        cat "${{CALLS[0]}}" | rastair_call_to_methylkit.sh | gzip -c \
            > {{outputs.results}}/"$P.rastair_methylkit.txt.gz" || exit 1
    """,
)

env = flyte.TaskEnvironment.from_task(
    "rastair",
    rastair_call_cmd.as_task(),
    rastair_mbias_cmd.as_task(),
    rastair_mbiasparser_cmd.as_task(),
    rastair_methylkit_cmd.as_task(),
)


async def pick(results: Dir, name: str, tool: str) -> File:
    found = await results.get_file(name)
    if found is None:
        raise FileNotFoundError(f"{tool} wrote no {name}")
    return found


async def rastair_mbias(bam: File, bai: File, fasta: File, fai: File, prefix: str, args: str = "") -> File:
    """The C-to-T conversion bias along the reads of an indexed BAM: ``<prefix>.rastair_mbias.txt``."""
    results = await rastair_mbias_cmd(bam=bam, bai=bai, fasta=fasta, fai=fai, prefix=prefix, args=args)
    return await pick(results, f"{prefix}.rastair_mbias.txt", "rastair mbias")


@dataclass
class RastairMbiasParserResult:
    csv: File  # <prefix>.rastair_mbias_processed.csv
    pdf: File | None  # <prefix>.rastair_mbias_processed.pdf
    trim_ot: str  # the CSV's first line: --nOT for rastair_call
    trim_ob: str  # its second line: --nOB
    results: Dir


async def rastair_mbiasparser(mbias: File, prefix: str) -> RastairMbiasParserResult:
    """Plot a :func:`rastair_mbias` table and pick the trimming for :func:`rastair_call`."""
    results = await rastair_mbiasparser_cmd(mbias=mbias, prefix=prefix)
    csv = await pick(results, f"{prefix}.rastair_mbias_processed.csv", "rastair mbias parser")
    async with csv.open("rb") as fh:
        lines = bytes(await fh.read()).decode().splitlines()
    if len(lines) < 2:
        raise ValueError(f"{csv.name} has {len(lines)} lines, expected the OT and OB trimming")
    return RastairMbiasParserResult(
        csv=csv,
        pdf=await results.get_file(f"{prefix}.rastair_mbias_processed.pdf"),
        trim_ot=lines[0],
        trim_ob=lines[1],
        results=results,
    )


async def rastair_call(
    bam: File, bai: File, fasta: File, fai: File, prefix: str, trim_ot: str, trim_ob: str, args: str = ""
) -> File:
    """Per-position methylation calls from an indexed BAM: ``<prefix>.rastair_call.txt``."""
    results = await rastair_call_cmd(
        bam=bam, bai=bai, fasta=fasta, fai=fai, prefix=prefix, trim_ot=trim_ot, trim_ob=trim_ob, args=args
    )
    return await pick(results, f"{prefix}.rastair_call.txt", "rastair call")


async def rastair_methylkit(calls: File, prefix: str) -> File:
    """:func:`rastair_call` calls in methylKit's format: ``<prefix>.rastair_methylkit.txt.gz``."""
    results = await rastair_methylkit_cmd(calls=calls, prefix=prefix)
    return await pick(results, f"{prefix}.rastair_methylkit.txt.gz", "rastair methylkit")
