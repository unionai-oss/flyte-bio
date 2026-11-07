"""quast — QUAST assembly quality assessment.

Exposes :func:`quast`, which evaluates one or more assemblies, optionally
against a reference FASTA (``-r``) and its gene annotation (``--features``,
GFF), and returns QUAST's report directory and its main tables.

The output Dir keeps upstream's layout and names:

- ``<prefix>/``: QUAST's whole report (``report.{tsv,txt,html,pdf,tex}``,
  ``quast.log``, ``basic_stats/``, ``icarus_viewers/``, and with a reference
  ``contigs_reports/``, ``aligned_stats/``, ``genome_stats/``);
- ``<prefix>.tsv``: ``<prefix>/report.tsv``;
- with a reference: ``<prefix>_misassemblies.tsv``, ``<prefix>_unaligned.tsv``
  and, for an assembly named ``transcriptome.*``, ``<prefix>_transcriptome.tsv``
  (copies of the matching ``contigs_reports/`` tables).

The assemblies keep their names, which QUAST uses as the report's column
headers. ``--threads`` follows the task's CPUs.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

QUAST_IMAGE = "community.wave.seqera.io/library/quast:5.3.0--755a216045b6dbdd"
DEFAULT_RESOURCES = flyte.Resources(cpu=4, memory="8Gi")

quast_cmd = shell.create(
    name="quast",
    image=QUAST_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "assemblies": list[File],
        "fasta": File | None,
        "gff": File | None,
        "prefix": str,
        "args": str,
    },
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        ASSEMBLIES=({{inputs.assemblies}})
        shopt -s nullglob; FASTA=({{inputs.fasta}}); GFF=({{inputs.gff}}); shopt -u nullglob
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"

        NAMES=()
        for F in "${{ASSEMBLIES[@]}}"; do
            ln -s "$F" "$(basename "$F")"
            NAMES+=("$(basename "$F")")
        done
        REFERENCE=()
        if [ ${{#FASTA[@]}} -gt 0 ]; then
            ln -s "${{FASTA[0]}}" "$(basename "${{FASTA[0]}}")"
            REFERENCE=(-r "$(basename "${{FASTA[0]}}")")
        fi
        FEATURES=()
        if [ ${{#GFF[@]}} -gt 0 ]; then
            ln -s "${{GFF[0]}}" "$(basename "${{GFF[0]}}")"
            FEATURES=(--features "$(basename "${{GFF[0]}}")")
        fi

        quast.py \
            --output-dir "$P" \
            "${{REFERENCE[@]}}" \
            "${{FEATURES[@]}}" \
            --threads $TASK_CPUS \
            $ARGS \
            "${{NAMES[@]}}"

        cp "$P/report.tsv" "$OUT/$P.tsv"
        for PAIR in all_alignments_transcriptome.tsv:transcriptome misassemblies_report.tsv:misassemblies \
                unaligned_report.tsv:unaligned; do
            SRC=$P/contigs_reports/${{PAIR%%:*}}
            if [ -f "$SRC" ]; then cp "$SRC" "$OUT/${{P}}_${{PAIR#*:}}.tsv"; fi
        done
        mv "$P" "$OUT"/
        true
    """,
)

env = flyte.TaskEnvironment.from_task("quast", quast_cmd.as_task())


@dataclass
class QuastResult:
    report: Dir  # QUAST's report directory, <prefix>/
    tsv: File  # <prefix>.tsv, the summary table
    misassemblies: File | None  # with a reference
    unaligned: File | None  # with a reference
    transcriptome: File | None  # with a reference, for an assembly named transcriptome.*
    results: Dir  # everything above


async def quast(
    assemblies: list[File],
    fasta: File | None = None,
    gff: File | None = None,
    prefix: str = "sample",
    args: str = "",
) -> QuastResult:
    """Assess ``assemblies``, against reference ``fasta`` (and its ``gff``) when given."""
    if not assemblies:
        raise ValueError("quast needs at least one assembly")
    results = await quast_cmd(assemblies=assemblies, fasta=fasta, gff=gff, prefix=prefix, args=args)
    tsv = await results.get_file(f"{prefix}.tsv")
    if tsv is None:
        raise FileNotFoundError(f"quast wrote no {prefix}.tsv")
    return QuastResult(
        report=Dir.from_existing_remote(results.path.rstrip("/") + "/" + prefix),
        tsv=tsv,
        misassemblies=await results.get_file(f"{prefix}_misassemblies.tsv"),
        unaligned=await results.get_file(f"{prefix}_unaligned.tsv"),
        transcriptome=await results.get_file(f"{prefix}_transcriptome.tsv"),
        results=results,
    )
