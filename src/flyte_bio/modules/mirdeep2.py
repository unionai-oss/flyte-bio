"""mirdeep2 — miRDeep2 small-RNA read mapping and miRNA discovery.

Exposes:

- :func:`mirdeep2_mapper` — ``mapper.pl``: preprocess reads (FASTA with
  ``-c``, FASTQ with ``-e``; adapter clipping, length filter, collapsing, as
  ``args`` asks) and map them to a genome with bowtie. Writes, as upstream,
  ``<prefix>_collapsed.fa`` and
  ``<prefix>_reads_collapsed_vs_<genome name>_genome.arf``. The bowtie index
  is a Dir of ``<name>.*.ebwt`` (e.g. from :func:`bowtie1_build`); the genome
  name defaults to the index's ``<name>``. bowtie gets the task's CPUs.
- :func:`mirdeep2_mirdeep2` — ``miRDeep2.pl``: score known and novel miRNAs
  from the mapper's collapsed reads and ``.arf`` against the genome, with
  optional mature, hairpin and other-species mature references (``none`` when
  unset, as upstream). miRDeep2 date-stamps its results; as upstream they are
  renamed to ``result_<prefix>.{bed,csv,html}``.

Both tools are perl, and a perl ``die`` exits 255, so a failure is turned into
exit 1. The inputs are linked into a working directory first, as miRDeep2
writes its run folders next to where it runs.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

MIRDEEP2_IMAGE = "quay.io/biocontainers/mirdeep2:2.0.1.2--0"
MAPPER_RESOURCES = flyte.Resources(cpu=4, memory="8Gi")
MIRDEEP2_RESOURCES = flyte.Resources(cpu=1, memory="6Gi")

mirdeep2_mapper_cmd = shell.create(
    name="mirdeep2_mapper",
    image=MIRDEEP2_IMAGE,
    resources=MAPPER_RESOURCES,
    inputs={"reads": File, "index": Dir, "genome_name": str, "prefix": str, "args": str},
    defaults={"genome_name": "", "args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        READS=({{inputs.reads}})
        GENOME={{inputs.genome_name}}
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}

        INDEX=$(find -L {{inputs.index}} -name "*.3.ebwt" ! -name "*.rev.*" | sed 's/\.3\.ebwt$//' | head -n 1)
        if [ -z "$INDEX" ]; then echo "Bowtie index files not found" >&2; exit 1; fi
        if [ -z "$GENOME" ]; then GENOME=$(basename "$INDEX"); fi

        W=$(mktemp -d); cd "$W"
        NAME=$(basename "${{READS[0]}}")
        ln -s "${{READS[0]}}" "$NAME"
        mapper.pl "$NAME" $ARGS -p "$INDEX" -o $TASK_CPUS \
            -s "${{P}}_collapsed.fa" -t "${{P}}_reads_collapsed_vs_${{GENOME}}_genome.arf" || exit 1
        mv "${{P}}_collapsed.fa" "${{P}}_reads_collapsed_vs_${{GENOME}}_genome.arf" "$OUT"/
    """,
)

mirdeep2_mirdeep2_cmd = shell.create(
    name="mirdeep2_mirdeep2",
    image=MIRDEEP2_IMAGE,
    resources=MIRDEEP2_RESOURCES,
    inputs={
        "reads": File,
        "genome_mappings": File,
        "fasta": File,
        "mature": File | None,
        "hairpin": File | None,
        "mature_other_species": File | None,
        "prefix": str,
        "args": str,
    },
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        READS=({inputs.reads}); ARF=({inputs.genome_mappings}); FA=({inputs.fasta})
        shopt -s nullglob
        MATURE=({inputs.mature}); HAIRPIN=({inputs.hairpin}); OTHER=({inputs.mature_other_species})
        shopt -u nullglob
        P={inputs.prefix}
        ARGS={inputs.args}
        OUT={outputs.results}

        W=$(mktemp -d); cd "$W"
        link() {
            if [ $# -eq 0 ]; then echo none; return; fi
            ln -s "$1" "$(basename "$1")"
            basename "$1"
        }
        R=$(link "${READS[0]}"); G=$(link "${FA[0]}"); M=$(link "${ARF[0]}")
        SPECIES=$(link "${MATURE[@]}"); PRECURSORS=$(link "${HAIRPIN[@]}"); OTHERS=$(link "${OTHER[@]}")

        miRDeep2.pl "$R" "$G" "$M" "$SPECIES" "$OTHERS" "$PRECURSORS" $ARGS || exit 1

        mv result_*.bed "$OUT/result_$P.bed"
        mv result_*.csv "$OUT/result_$P.csv"
        mv result_*.html "$OUT/result_$P.html"
    """,
)

env = flyte.TaskEnvironment.from_task(
    "mirdeep2",
    mirdeep2_mapper_cmd.as_task(),
    mirdeep2_mirdeep2_cmd.as_task(),
)


async def pick(results: Dir, suffix: str, tool: str) -> File:
    found = [f async for f in results.walk() if f.path.endswith(suffix)]
    if len(found) != 1:
        raise FileNotFoundError(f"{tool} wrote {len(found)} files ending in {suffix}, expected 1")
    return found[0]


@dataclass
class Mirdeep2MapperResult:
    collapsed: File  # <prefix>_collapsed.fa
    arf: File  # <prefix>_reads_collapsed_vs_<genome name>_genome.arf
    results: Dir


async def mirdeep2_mapper(
    reads: File, index: Dir, prefix: str = "sample", genome_name: str = "", args: str = ""
) -> Mirdeep2MapperResult:
    """Preprocess and map ``reads`` to the bowtie ``index`` with ``mapper.pl``.

    ``args`` carries the input format and preprocessing (e.g. ``-c -j -k
    <adapter> -l 18 -m``). ``genome_name`` names the ``.arf`` (upstream's
    index id); it defaults to the index's file prefix.
    """
    results = await mirdeep2_mapper_cmd(reads=reads, index=index, genome_name=genome_name, prefix=prefix, args=args)
    return Mirdeep2MapperResult(
        collapsed=await pick(results, f"{prefix}_collapsed.fa", "mapper.pl"),
        arf=await pick(results, "_genome.arf", "mapper.pl"),
        results=results,
    )


@dataclass
class Mirdeep2Result:
    bed: File  # result_<prefix>.bed
    csv: File  # result_<prefix>.csv
    html: File  # result_<prefix>.html
    results: Dir


async def mirdeep2_mirdeep2(
    reads: File,
    genome_mappings: File,
    fasta: File,
    mature: File | None = None,
    hairpin: File | None = None,
    mature_other_species: File | None = None,
    prefix: str = "sample",
    args: str = "",
) -> Mirdeep2Result:
    """Run ``miRDeep2.pl`` on the mapper's collapsed ``reads`` and ``genome_mappings`` (.arf)."""
    results = await mirdeep2_mirdeep2_cmd(
        reads=reads,
        genome_mappings=genome_mappings,
        fasta=fasta,
        mature=mature,
        hairpin=hairpin,
        mature_other_species=mature_other_species,
        prefix=prefix,
        args=args,
    )
    return Mirdeep2Result(
        bed=await pick(results, f"result_{prefix}.bed", "miRDeep2.pl"),
        csv=await pick(results, f"result_{prefix}.csv", "miRDeep2.pl"),
        html=await pick(results, f"result_{prefix}.html", "miRDeep2.pl"),
        results=results,
    )
