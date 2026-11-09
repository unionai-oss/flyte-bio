"""links — LINKS, scaffolding a draft assembly with long reads (or long sequences).

Exposes :func:`links`, which scaffolds an assembly (FASTA, plain or gzipped)
with one or more read files (gzipped ones are decompressed first), as
upstream. The outputs keep LINKS' names under ``-b <prefix>``:

- ``<prefix>.scaffolds.fa`` and the scaffold layout ``<prefix>.scaffolds``;
- ``<prefix>.assembly_correspondence.tsv``, ``<prefix>.gv`` (scaffold graph),
  ``<prefix>.bloom`` (the reads' Bloom filter);
- ``<prefix>.log``, ``<prefix>.pairing_issues``,
  ``<prefix>.tigpair_checkpoint.tsv`` and, when written,
  ``<prefix>.pairing_distribution.csv`` and
  ``<prefix>.simplepair_checkpoint.tsv``.

Pass LINKS' options (``-d`` distances, ``-k``, ``-l``, …) in ``args``. As
upstream, ``-j`` is the task's CPUs, at most 4 (LINKS supports no more).
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

LINKS_IMAGE = "quay.io/biocontainers/links:2.0.1--h4ac6f70_5"
DEFAULT_RESOURCES = flyte.Resources(cpu=4, memory="16Gi")

links_cmd = shell.create(
    name="links",
    image=LINKS_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"assembly": File, "reads": list[File], "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        ASM=({{inputs.assembly}})
        READS=({{inputs.reads}})
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        if [[ "${{ASM[0]}}" == *.gz ]]; then
            gzip -dc "${{ASM[0]}}" > assembly.fa
        else
            ln -s "${{ASM[0]}}" assembly.fa
        fi
        for READ in "${{READS[@]}}"; do
            NAME=$(basename "$READ")
            if [[ "$NAME" == *.gz ]]; then
                gzip -dc "$READ" > "${{NAME%.gz}}"
                echo "${{NAME%.gz}}" >> readfile.fof
            else
                ln -s "$READ" "$NAME"
                echo "$NAME" >> readfile.fof
            fi
        done
        THREADS=$(( TASK_CPUS < 4 ? TASK_CPUS : 4 ))
        LINKS -f assembly.fa -s readfile.fof -j $THREADS -b "$P" $ARGS || exit 1
        shopt -s nullglob
        for F in "$P".log "$P".pairing_distribution.csv "$P".pairing_issues "$P".scaffolds "$P".scaffolds.fa \
                "$P".bloom "$P".gv "$P".assembly_correspondence.tsv "$P".simplepair_checkpoint.tsv \
                "$P".tigpair_checkpoint.tsv; do
            if [ -e "$F" ]; then cp -L "$F" "$OUT"/; fi
        done
    """,
)

env = flyte.TaskEnvironment.from_task("links", links_cmd.as_task())


@dataclass
class LinksResult:
    scaffolds_fasta: File  # <prefix>.scaffolds.fa
    scaffolds_csv: File  # <prefix>.scaffolds
    assembly_correspondence: File  # <prefix>.assembly_correspondence.tsv
    scaffolds_graph: File  # <prefix>.gv
    bloom: File  # <prefix>.bloom
    log: File  # <prefix>.log
    pairing_issues: File  # <prefix>.pairing_issues
    tigpair_checkpoint: File  # <prefix>.tigpair_checkpoint.tsv
    pairing_distribution: File | None  # <prefix>.pairing_distribution.csv
    simplepair_checkpoint: File | None  # <prefix>.simplepair_checkpoint.tsv
    results: Dir  # everything above


async def links(assembly: File, reads: list[File], prefix: str = "sample", args: str = "") -> LinksResult:
    """Scaffold ``assembly`` with ``reads``."""
    if not reads:
        raise ValueError("LINKS needs reads")
    results = await links_cmd(assembly=assembly, reads=reads, prefix=prefix, args=args)

    async def need(suffix: str) -> File:
        found = await results.get_file(f"{prefix}.{suffix}")
        if found is None:
            raise FileNotFoundError(f"LINKS wrote no {prefix}.{suffix}")
        return found

    return LinksResult(
        scaffolds_fasta=await need("scaffolds.fa"),
        scaffolds_csv=await need("scaffolds"),
        assembly_correspondence=await need("assembly_correspondence.tsv"),
        scaffolds_graph=await need("gv"),
        bloom=await need("bloom"),
        log=await need("log"),
        pairing_issues=await need("pairing_issues"),
        tigpair_checkpoint=await need("tigpair_checkpoint.tsv"),
        pairing_distribution=await results.get_file(f"{prefix}.pairing_distribution.csv"),
        simplepair_checkpoint=await results.get_file(f"{prefix}.simplepair_checkpoint.tsv"),
        results=results,
    )
