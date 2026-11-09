"""pilon — polish a draft assembly with aligned reads.

Exposes :func:`pilon`, which improves an assembly (FASTA) from a
coordinate-sorted, indexed BAM of reads aligned to it. ``pilon_mode`` says
what the BAM holds, as Pilon's option of that name: ``frags`` (paired-end
fragments), ``jumps`` (mate pairs), ``unpaired`` or ``bam`` (let Pilon
decide). As upstream, the outputs are named by ``--output <prefix>``:

- ``<prefix>.fasta``: the improved assembly;
- with the matching options in ``args``: ``<prefix>.vcf`` (``--vcf``),
  ``<prefix>.changes`` (``--changes``; upstream's ``*.change`` pattern never
  matches it) and the ``.bed`` / ``.wig`` tracks (``--tracks``).

The JVM heap (``-Xmx``) is the task's memory minus 1 GiB, or 80% of it
below 2 GiB, as upstream.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

PILON_IMAGE = "quay.io/biocontainers/pilon:1.24--hdfd78af_0"
DEFAULT_RESOURCES = flyte.Resources(cpu=4, memory="16Gi")
PILON_MODES = ("frags", "jumps", "unpaired", "bam")

pilon_cmd = shell.create(
    name="pilon",
    image=PILON_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"fasta": File, "bam": File, "bai": File, "pilon_mode": str, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        FA=({{inputs.fasta}}); BAM=({{inputs.bam}}); BAI=({{inputs.bai}})
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        # The BAM needs its index beside it; inputs go in their own directory, apart from the outputs.
        mkdir in
        FASTA=in/$(basename "${{FA[0]}}"); ALN=in/$(basename "${{BAM[0]}}")
        ln -s "${{FA[0]}}" "$FASTA"
        ln -s "${{BAM[0]}}" "$ALN"
        ln -s "${{BAI[0]}}" "$ALN.bai"
        MEM_MB=$(( TASK_MEMORY_MB - 1024 ))
        if [ "$TASK_MEMORY_MB" -lt 2048 ]; then MEM_MB=$(( TASK_MEMORY_MB * 8 / 10 )); fi
        PILON_JAR=$(echo "$(dirname "$(which pilon)")"/../share/pilon*/pilon.jar)
        java -Xmx${{MEM_MB}}M -jar $PILON_JAR \
            --genome "$FASTA" \
            --output {{inputs.prefix}} \
            $ARGS \
            --{{inputs.pilon_mode}} "$ALN" || exit 1
        shopt -s nullglob
        mv ./*.fasta ./*.vcf ./*.changes ./*.bed ./*.wig "$OUT"/
    """,
)

env = flyte.TaskEnvironment.from_task("pilon", pilon_cmd.as_task())


@dataclass
class PilonResult:
    improved_assembly: File  # <prefix>.fasta
    vcf: File | None  # <prefix>.vcf, with --vcf
    change_record: File | None  # <prefix>.changes, with --changes
    tracks: list[File]  # *.bed and *.wig, with --tracks
    results: Dir  # everything above


async def pilon(
    fasta: File, bam: File, bai: File, pilon_mode: str = "bam", prefix: str = "sample", args: str = ""
) -> PilonResult:
    """Polish ``fasta`` with the reads in ``bam`` (indexed by ``bai``)."""
    if pilon_mode not in PILON_MODES:
        raise ValueError(f"Unrecognised mode to run Pilon: {pilon_mode!r}. Options: {', '.join(PILON_MODES)}")
    results = await pilon_cmd(fasta=fasta, bam=bam, bai=bai, pilon_mode=pilon_mode, prefix=prefix, args=args)
    files = {f.name: f async for f in results.walk()}
    if f"{prefix}.fasta" not in files:
        raise FileNotFoundError(f"Pilon wrote no {prefix}.fasta")
    return PilonResult(
        improved_assembly=files[f"{prefix}.fasta"],
        vcf=files.get(f"{prefix}.vcf"),
        change_record=files.get(f"{prefix}.changes"),
        tracks=[files[n] for n in sorted(files) if n.endswith((".bed", ".wig"))],
        results=results,
    )
