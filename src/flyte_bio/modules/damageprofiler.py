"""damageprofiler — ancient-DNA damage patterns of aligned reads (DamageProfiler).

Exposes :func:`damageprofiler`, which profiles one sample's BAM: C->T and G->A
misincorporation frequencies at the 5' and 3' read ends, the read-length
distribution and edit distances. A reference FASTA (with its ``.fai``) is
optional, as is a species list (``-sf``: reference names, one per line) that
profiles each listed reference separately.

As upstream, the output is DamageProfiler's output directory (upstream's
``<prefix>/``), returned as a Dir holding its files under DamageProfiler's
names: ``5pCtoT_freq.txt``, ``3pGtoA_freq.txt``,
``{5p,3p}_freq_misincorporations.txt``, ``DNA_comp_genome.txt``,
``DNA_composition_sample.txt``, ``misincorporation.txt``,
``lgdistribution.txt``, ``editDistance.txt``, ``dmgprof.json``,
``DamageProfiler.log``, and the plots (``DamagePlot.pdf``,
``DamagePlot_{five,three}_prime.svg``, ``Length_plot*.{pdf,svg}``,
``edit_distance.{pdf,svg}``). With a species list, each species' files are in
a subdirectory named after it.

DamageProfiler is Java: its heap is 80% of the task's memory (upstream's
bioconda wrapper defaults to 1 GiB).
"""

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import JAVA_HEAP_MB, TASK_LIMITS

DAMAGEPROFILER_IMAGE = "quay.io/biocontainers/damageprofiler:1.1--hdfd78af_2"
DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory="6Gi")

damageprofiler_cmd = shell.create(
    name="damageprofiler",
    image=DAMAGEPROFILER_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "bam": File,
        "fasta": File | None,
        "fai": File | None,
        "species_list": File | None,
        "prefix": str,
        "args": str,
    },
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        BAM=({{inputs.bam}})
        shopt -s nullglob
        FA=({{inputs.fasta}}); FAI=({{inputs.fai}}); SPECIES=({{inputs.species_list}})
        shopt -u nullglob
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        export HOME="$W"

        ln -s "${{BAM[0]}}" "$(basename "${{BAM[0]}}")"
        EXTRA=()
        if [ ${{#FA[@]}} -gt 0 ]; then
            ln -s "${{FA[0]}}" "$(basename "${{FA[0]}}")"
            EXTRA+=(-r "$(basename "${{FA[0]}}")")
        fi
        if [ ${{#FAI[@]}} -gt 0 ]; then ln -s "${{FAI[0]}}" "$(basename "${{FAI[0]}}")"; fi
        if [ ${{#SPECIES[@]}} -gt 0 ]; then EXTRA+=(-sf "${{SPECIES[0]}}"); fi

        damageprofiler -Xmx{JAVA_HEAP_MB}m \
            -i "$(basename "${{BAM[0]}}")" \
            -o "$P/" \
            $ARGS \
            "${{EXTRA[@]}}" || exit 1
        mv "$P"/* "$OUT"/
    """,
)

env = flyte.TaskEnvironment.from_task("damageprofiler", damageprofiler_cmd.as_task())


async def damageprofiler(
    bam: File,
    fasta: File | None = None,
    fai: File | None = None,
    species_list: File | None = None,
    prefix: str = "sample",
    args: str = "",
) -> Dir:
    """DamageProfiler's output directory for ``bam`` (upstream's ``<prefix>/``)."""
    return await damageprofiler_cmd(bam=bam, fasta=fasta, fai=fai, species_list=species_list, prefix=prefix, args=args)
