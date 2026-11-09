"""methyldackel — methylation calls and M-bias from bisulfite alignments with MethylDackel.

Exposes:

- :func:`methyldackel_extract` — per-cytosine methylation (``MethylDackel
  extract``): ``<bam name>_CpG.bedGraph`` by default, ``.methylKit`` with
  ``--methylKit``, and ``_CHG`` / ``_CHH`` files with ``--CHG`` / ``--CHH``.
- :func:`methyldackel_mbias` — methylation bias along the reads
  (``MethylDackel mbias --txt``): ``<prefix>.mbias.txt``, with its suggested
  trimming on the first line.

As upstream, the BAM and its index, and the FASTA and its ``.fai``, are linked
side by side under their own names; ``extract`` names its outputs after the
BAM. ``args`` go to MethylDackel (e.g. ``-@ <threads>``).
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

METHYLDACKEL_IMAGE = "quay.io/biocontainers/methyldackel:0.6.1--he4a0461_7"

# Upstream labels extract process_medium and mbias process_low; MethylDackel
# runs one thread unless args set -@.
MEDIUM_RESOURCES = flyte.Resources(cpu=2, memory="8Gi")
LOW_RESOURCES = flyte.Resources(cpu=1, memory="4Gi")

# Links the BAM, its index, the FASTA and its .fai into the working directory
# under their own names, and sets BAM_NAME and FASTA_NAME.
LINK_INPUTS = r"""
        W=$(mktemp -d); cd "$W"
        for F in "${BAM[0]}" "${BAI[0]}" "${FA[0]}" "${FAI[0]}"; do ln -s "$F" "$(basename "$F")"; done
        BAM_NAME=$(basename "${BAM[0]}")
        FASTA_NAME=$(basename "${FA[0]}")
"""

INPUTS = {"bam": File, "bai": File, "fasta": File, "fai": File}

methyldackel_extract_cmd = shell.create(
    name="methyldackel_extract",
    image=METHYLDACKEL_IMAGE,
    resources=MEDIUM_RESOURCES,
    inputs={**INPUTS, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        BAM=({{inputs.bam}}); BAI=({{inputs.bai}}); FA=({{inputs.fasta}}); FAI=({{inputs.fai}})
        ARGS={{inputs.args}}
        {LINK_INPUTS}
        MethylDackel extract $ARGS "$FASTA_NAME" "$BAM_NAME"
        shopt -s nullglob
        mv ./*.bedGraph ./*.methylKit {{outputs.results}}/
    """,
)

methyldackel_mbias_cmd = shell.create(
    name="methyldackel_mbias",
    image=METHYLDACKEL_IMAGE,
    resources=LOW_RESOURCES,
    inputs={**INPUTS, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        BAM=({{inputs.bam}}); BAI=({{inputs.bai}}); FA=({{inputs.fasta}}); FAI=({{inputs.fai}})
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        {LINK_INPUTS}
        MethylDackel mbias $ARGS "$FASTA_NAME" "$BAM_NAME" "$P" --txt > {{outputs.results}}/"$P.mbias.txt"
    """,
)

env = flyte.TaskEnvironment.from_task(
    "methyldackel", methyldackel_extract_cmd.as_task(), methyldackel_mbias_cmd.as_task()
)


@dataclass
class MethyldackelExtractResult:
    bedgraph: list[File]  # *.bedGraph: <bam name>_CpG.bedGraph by default
    methylkit: list[File]  # *.methylKit, with --methylKit
    results: Dir


async def methyldackel_extract(
    bam: File, bai: File, fasta: File, fai: File, args: str = ""
) -> MethyldackelExtractResult:
    """Per-cytosine methylation calls from a coordinate-sorted, indexed bisulfite BAM."""
    results = await methyldackel_extract_cmd(bam=bam, bai=bai, fasta=fasta, fai=fai, args=args)
    files = sorted([f async for f in results.walk()], key=lambda f: f.name)
    return MethyldackelExtractResult(
        bedgraph=[f for f in files if f.name.endswith(".bedGraph")],
        methylkit=[f for f in files if f.name.endswith(".methylKit")],
        results=results,
    )


async def methyldackel_mbias(bam: File, bai: File, fasta: File, fai: File, prefix: str, args: str = "") -> File:
    """The methylation bias along the reads of an indexed bisulfite BAM: ``<prefix>.mbias.txt``."""
    results = await methyldackel_mbias_cmd(bam=bam, bai=bai, fasta=fasta, fai=fai, prefix=prefix, args=args)
    txt = await results.get_file(f"{prefix}.mbias.txt")
    if txt is None:
        raise FileNotFoundError(f"MethylDackel mbias wrote no {prefix}.mbias.txt")
    return txt
