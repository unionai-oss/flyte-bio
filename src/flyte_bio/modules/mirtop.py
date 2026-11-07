"""mirtop — miRNA and isomiR annotation in the mirGFF3 format.

Exposes:

- :func:`mirtop_gff` — annotate aligned small-RNA reads (one or more BAMs,
  against the hairpin FASTA) as mirGFF3: ``<prefix>_mirtop.gff``, with one
  sample column per BAM in the caller's order.
- :func:`mirtop_counts` — a count table of a mirGFF3 file (``--add-extra``
  in ``args`` adds the sequence columns): mirtop names it after the GFF
  (``<gff name>.tsv``).
- :func:`mirtop_export` — export a mirGFF3 file to another format (``args``,
  default ``--format isomir``): ``*_rawData.tsv`` (isomir), ``*.fasta`` or
  ``*.vcf``, as mirtop names them.
- :func:`mirtop_stats` — summary statistics of a mirGFF3 file:
  ``<prefix>_stats.txt`` and ``<prefix>_stats.log``.

``gtf`` is the miRBase GFF3 of the species, and ``species`` its miRBase code
(``hsa``).
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

MIRTOP_IMAGE = "community.wave.seqera.io/library/mirtop_pybedtools_pysam_samtools:b9705c2683e775b8"
DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory="6Gi")

# A list[File] is staged under the files' names and globbed alphabetically, so
# ``names`` carries the caller's order (one name per line) into the script.
mirtop_gff_cmd = shell.create(
    name="mirtop_gff",
    image=MIRTOP_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "bams": list[File],
        "names": str,
        "hairpin": File,
        "gtf": File,
        "species": str,
        "prefix": str,
        "args": str,
    },
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        BAMS=({inputs.bams}); HAIRPIN=({inputs.hairpin}); GTF=({inputs.gtf})
        NAMES={inputs.names}
        SPECIES={inputs.species}
        P={inputs.prefix}
        ARGS={inputs.args}
        OUT={outputs.results}
        D=$(dirname "${BAMS[0]}")
        W=$(mktemp -d); cd "$W"
        FILES=()
        while IFS= read -r NAME; do ln -s "$D/$NAME" "$NAME"; FILES+=("$NAME"); done <<< "$NAMES"

        mirtop gff $ARGS --sps "$SPECIES" --hairpin "${HAIRPIN[0]}" --gtf "${GTF[0]}" -o mirtop "${FILES[@]}"
        mv mirtop/mirtop.gff "$OUT/${P}_mirtop.gff"
    """,
)

mirtop_counts_cmd = shell.create(
    name="mirtop_counts",
    image=MIRTOP_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"gff": File, "hairpin": File, "gtf": File, "species": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        GFF=({inputs.gff}); HAIRPIN=({inputs.hairpin}); GTF=({inputs.gtf})
        SPECIES={inputs.species}
        ARGS={inputs.args}
        OUT={outputs.results}
        W=$(mktemp -d); cd "$W"
        G=$(basename "${GFF[0]}"); ln -s "${GFF[0]}" "$G"

        mirtop counts $ARGS --hairpin "${HAIRPIN[0]}" --gtf "${GTF[0]}" --sps "$SPECIES" --gff "$G" -o counts
        mv counts/*.tsv "$OUT"/
    """,
)

mirtop_export_cmd = shell.create(
    name="mirtop_export",
    image=MIRTOP_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"gff": File, "hairpin": File, "gtf": File, "species": str, "args": str},
    defaults={"args": "--format isomir"},
    outputs={"results": Dir},
    script=r"""
        GFF=({inputs.gff}); HAIRPIN=({inputs.hairpin}); GTF=({inputs.gtf})
        SPECIES={inputs.species}
        ARGS={inputs.args}
        OUT={outputs.results}
        W=$(mktemp -d); cd "$W"
        G=$(basename "${GFF[0]}"); ln -s "${GFF[0]}" "$G"

        mirtop export $ARGS --hairpin "${HAIRPIN[0]}" --gtf "${GTF[0]}" --sps "$SPECIES" -o export "$G"
        shopt -s nullglob
        WRITTEN=(export/*_rawData.tsv export/*.fasta export/*.vcf*)
        if [ ${#WRITTEN[@]} -gt 0 ]; then mv "${WRITTEN[@]}" "$OUT"/; fi
    """,
)

mirtop_stats_cmd = shell.create(
    name="mirtop_stats",
    image=MIRTOP_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"gff": File, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        GFF=({inputs.gff})
        P={inputs.prefix}
        ARGS={inputs.args}
        OUT={outputs.results}
        W=$(mktemp -d); cd "$W"
        G=$(basename "${GFF[0]}"); ln -s "${GFF[0]}" "$G"

        mirtop stats $ARGS --out stats "$G"
        mv stats/mirtop_stats.log "$OUT/${P}_stats.log"
        mv stats/mirtop_stats.txt "$OUT/${P}_stats.txt"
    """,
)

env = flyte.TaskEnvironment.from_task(
    "mirtop",
    mirtop_gff_cmd.as_task(),
    mirtop_counts_cmd.as_task(),
    mirtop_export_cmd.as_task(),
    mirtop_stats_cmd.as_task(),
)


def file_name(file: File) -> str:
    return file.path.rstrip("/").rsplit("/", 1)[-1]


async def pick(results: Dir, name: str, tool: str) -> File:
    found = await results.get_file(name)
    if found is None:
        raise FileNotFoundError(f"{tool} wrote no {name}")
    return found


async def mirtop_gff(
    bams: list[File], hairpin: File, gtf: File, species: str, prefix: str = "sample", args: str = ""
) -> File:
    """mirGFF3 annotation of ``bams`` (in list order): ``<prefix>_mirtop.gff``."""
    names = [file_name(b) for b in bams]
    if not names:
        raise ValueError("mirtop gff needs at least one BAM")
    if len(set(names)) != len(names):
        raise ValueError(f"mirtop gff inputs need distinct names, got {names}")
    results = await mirtop_gff_cmd(
        bams=bams, names="\n".join(names), hairpin=hairpin, gtf=gtf, species=species, prefix=prefix, args=args
    )
    return await pick(results, f"{prefix}_mirtop.gff", "mirtop gff")


async def mirtop_counts(gff: File, hairpin: File, gtf: File, species: str, args: str = "") -> File:
    """The count table of ``gff``, named after it (``<gff name>.tsv``)."""
    results = await mirtop_counts_cmd(gff=gff, hairpin=hairpin, gtf=gtf, species=species, args=args)
    tables = await results.list_files()
    if len(tables) != 1:
        raise FileNotFoundError(f"mirtop counts wrote {len(tables)} tables, expected 1")
    return tables[0]


@dataclass
class MirtopExportResult:
    tsv: list[File]  # *_rawData.tsv (isomir format)
    fasta: list[File]  # *.fasta (fasta format)
    vcf: list[File]  # *.vcf* (vcf format)
    results: Dir


async def mirtop_export(
    gff: File, hairpin: File, gtf: File, species: str, args: str = "--format isomir"
) -> MirtopExportResult:
    """Export ``gff`` to the format ``args`` names."""
    results = await mirtop_export_cmd(gff=gff, hairpin=hairpin, gtf=gtf, species=species, args=args)
    files = sorted(await results.list_files(), key=file_name)
    return MirtopExportResult(
        tsv=[f for f in files if file_name(f).endswith("_rawData.tsv")],
        fasta=[f for f in files if file_name(f).endswith(".fasta")],
        vcf=[f for f in files if ".vcf" in file_name(f)],
        results=results,
    )


@dataclass
class MirtopStatsResult:
    txt: File  # <prefix>_stats.txt
    log: File  # <prefix>_stats.log
    results: Dir


async def mirtop_stats(gff: File, prefix: str = "sample", args: str = "") -> MirtopStatsResult:
    """Summary statistics of ``gff``."""
    results = await mirtop_stats_cmd(gff=gff, prefix=prefix, args=args)
    return MirtopStatsResult(
        txt=await pick(results, f"{prefix}_stats.txt", "mirtop stats"),
        log=await pick(results, f"{prefix}_stats.log", "mirtop stats"),
        results=results,
    )
