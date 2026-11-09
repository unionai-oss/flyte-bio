"""arriba — gene-fusion detection from STAR's chimeric alignments, and fusion plots.

Exposes:

- :func:`arriba_arriba`, which calls fusions from a STAR BAM (aligned with
  chimeric reads ``WithinBAM``) against the genome FASTA and GTF, writing
  ``<prefix>.fusions.tsv`` and the discarded candidates,
  ``<prefix>.fusions.discarded.tsv``. Without a blacklist it runs with
  ``-f blacklist`` (the blacklist filter off), as upstream; known fusions,
  cytobands and protein domains are optional.
- :func:`arriba_visualisation`, which draws the fusions of an Arriba table
  (``draw_fusions.R``) into ``<prefix>.pdf``, with the reads' coverage when
  given the BAM (and its index), and the cytobands and protein domains when
  given.

Both share one image (Arriba 2.5.0). The database files (blacklists, known
fusions, cytobands, protein domains) are Arriba's release ``database/`` files.
"""

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

ARRIBA_IMAGE = "community.wave.seqera.io/library/arriba_wget:a3e48cf793a0b654"
# Upstream's process_medium, capped: Arriba holds the genome and annotation in memory.
DEFAULT_RESOURCES = flyte.Resources(cpu=4, memory="16Gi")

arriba_arriba_cmd = shell.create(
    name="arriba_arriba",
    image=ARRIBA_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "bam": File,
        "fasta": File,
        "gtf": File,
        "blacklist": File | None,
        "known_fusions": File | None,
        "cytobands": File | None,
        "protein_domains": File | None,
        "prefix": str,
        "args": str,
    },
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        BAM=({{inputs.bam}}); FA=({{inputs.fasta}}); GTF=({{inputs.gtf}})
        shopt -s nullglob
        BLACKLIST=({{inputs.blacklist}}); KNOWN=({{inputs.known_fusions}})
        CYTOBANDS=({{inputs.cytobands}}); DOMAINS=({{inputs.protein_domains}})
        shopt -u nullglob
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        W=$(mktemp -d); cd "$W"
        OPTIONAL=(-f blacklist)
        if [ ${{#BLACKLIST[@]}} -gt 0 ]; then OPTIONAL=(-b "${{BLACKLIST[0]}}"); fi
        if [ ${{#KNOWN[@]}} -gt 0 ]; then OPTIONAL+=(-k "${{KNOWN[0]}}"); fi
        if [ ${{#CYTOBANDS[@]}} -gt 0 ]; then OPTIONAL+=(-d "${{CYTOBANDS[0]}}"); fi
        if [ ${{#DOMAINS[@]}} -gt 0 ]; then OPTIONAL+=(-p "${{DOMAINS[0]}}"); fi
        arriba \
            -x "${{BAM[0]}}" \
            -a "${{FA[0]}}" \
            -g "${{GTF[0]}}" \
            -o "$P.fusions.tsv" \
            -O "$P.fusions.discarded.tsv" \
            "${{OPTIONAL[@]}}" \
            $ARGS || exit 1
        mv "$P.fusions.tsv" "$P.fusions.discarded.tsv" {{outputs.results}}/
    """,
)

# As upstream: the alignments (with the index beside them), cytobands and
# protein domains are optional.
arriba_visualisation_cmd = shell.create(
    name="arriba_visualisation",
    image=ARRIBA_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "fusions": File,
        "gtf": File,
        "bam": File | None,
        "bai": File | None,
        "protein_domains": File | None,
        "cytobands": File | None,
        "prefix": str,
        "args": str,
    },
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        FUSIONS=({{inputs.fusions}}); GTF=({{inputs.gtf}})
        shopt -s nullglob
        BAM=({{inputs.bam}}); BAI=({{inputs.bai}})
        DOMAINS=({{inputs.protein_domains}}); CYTOBANDS=({{inputs.cytobands}})
        shopt -u nullglob
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        W=$(mktemp -d); cd "$W"
        OPTIONAL=()
        if [ ${{#BAM[@]}} -gt 0 ]; then
            BAM_NAME=$(basename "${{BAM[0]}}")
            ln -s "${{BAM[0]}}" "$BAM_NAME"
            if [ ${{#BAI[@]}} -gt 0 ]; then ln -s "${{BAI[0]}}" "$BAM_NAME.bai"; fi
            OPTIONAL+=(--alignments="$BAM_NAME")
        fi
        if [ ${{#CYTOBANDS[@]}} -gt 0 ]; then OPTIONAL+=(--cytobands="${{CYTOBANDS[0]}}"); fi
        if [ ${{#DOMAINS[@]}} -gt 0 ]; then OPTIONAL+=(--proteinDomains="${{DOMAINS[0]}}"); fi
        draw_fusions.R \
            --fusions="${{FUSIONS[0]}}" \
            --output="$P.pdf" \
            --annotation="${{GTF[0]}}" \
            "${{OPTIONAL[@]}}" \
            $ARGS || exit 1
        mv "$P.pdf" {{outputs.results}}/
    """,
)

env = flyte.TaskEnvironment.from_task(
    "arriba",
    arriba_arriba_cmd.as_task(),
    arriba_visualisation_cmd.as_task(),
)


async def pick(results: Dir, name: str, tool: str) -> File:
    found = await results.get_file(name)
    if found is None:
        raise FileNotFoundError(f"{tool} wrote no {name}")
    return found


async def arriba_arriba(
    bam: File,
    fasta: File,
    gtf: File,
    prefix: str,
    blacklist: File | None = None,
    known_fusions: File | None = None,
    cytobands: File | None = None,
    protein_domains: File | None = None,
    args: str = "",
) -> tuple[File, File]:
    """Fusions in a STAR ``bam``: (``<prefix>.fusions.tsv``, ``<prefix>.fusions.discarded.tsv``).

    Without ``blacklist``, the blacklist filter is off (``-f blacklist``), as upstream.
    """
    results = await arriba_arriba_cmd(
        bam=bam,
        fasta=fasta,
        gtf=gtf,
        blacklist=blacklist,
        known_fusions=known_fusions,
        cytobands=cytobands,
        protein_domains=protein_domains,
        prefix=prefix,
        args=args,
    )
    return (
        await pick(results, f"{prefix}.fusions.tsv", "arriba"),
        await pick(results, f"{prefix}.fusions.discarded.tsv", "arriba"),
    )


async def arriba_visualisation(
    fusions: File,
    gtf: File,
    prefix: str,
    bam: File | None = None,
    bai: File | None = None,
    protein_domains: File | None = None,
    cytobands: File | None = None,
    args: str = "",
) -> File:
    """Plots of the fusions in an Arriba ``fusions`` table: ``<prefix>.pdf``.

    ``bam`` (coordinate-sorted, with its ``bai``) adds the read coverage.
    """
    results = await arriba_visualisation_cmd(
        fusions=fusions,
        gtf=gtf,
        bam=bam,
        bai=bai,
        protein_domains=protein_domains,
        cytobands=cytobands,
        prefix=prefix,
        args=args,
    )
    return await pick(results, f"{prefix}.pdf", "draw_fusions.R")
