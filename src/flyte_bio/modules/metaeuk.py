"""metaeuk — eukaryotic gene prediction in contigs with MetaEuk.

Exposes :func:`metaeuk_easypredict`, ``metaeuk easy-predict`` of a
nucleotide FASTA (contigs, a genome) against reference proteins: either a
protein FASTA or an MMseqs2 database (a Dir as :func:`flyte_bio.modules.mmseqs.mmseqs_databases`
or ``mmseqs_createdb`` writes it, found as upstream by its ``*.version`` file,
else its ``*.dbtype`` file). Outputs keep MetaEuk's names under ``<prefix>``:

- ``<prefix>.fas``: the predicted proteins;
- ``<prefix>.codon.fas``: their coding sequences;
- ``<prefix>.headersMap.tsv``: each prediction's header and its sources;
- ``<prefix>.gff``: the predictions as GFF.

MetaEuk runs with ``--threads`` of the task's CPUs; ``args`` passes its
options (``--min-length``, ``-e``, ``--metaeuk-eval`` and so on).
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

METAEUK_IMAGE = "quay.io/biocontainers/metaeuk:6.a5d39d9--pl5321hf1761c0_2"
# Upstream's process_medium (at the 4 CPUs upstream's tests run with). Large
# reference databases need more memory and disk for MMseqs2's temporary files.
DEFAULT_RESOURCES = flyte.Resources(cpu=4, memory="8Gi", disk="20Gi")

metaeuk_easypredict_cmd = shell.create(
    name="metaeuk_easypredict",
    image=METAEUK_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"fasta": File, "database": File | None, "database_dir": Dir | None, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        FA=({{inputs.fasta}})
        shopt -s nullglob; DATABASE=({{inputs.database}}); shopt -u nullglob
        DATABASE_DIR={{inputs.database_dir}}
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        if [ ${{#DATABASE[@]}} -gt 0 ]; then
            DB="${{DATABASE[0]}}"
        else
            # MetaEuk needs the database's base name, as upstream finds it.
            MARKER=$(find -L "$DATABASE_DIR" -maxdepth 1 -name '*.version' | sort | head -n 1)
            if [ -z "$MARKER" ]; then
                MARKER=$(find -L "$DATABASE_DIR" -maxdepth 1 -name '*.dbtype' | sort | head -n 1)
            fi
            if [ -z "$MARKER" ]; then echo "no MMseqs2 database in $DATABASE_DIR" >&2; exit 1; fi
            BASE=$(basename "$MARKER")
            DB="$DATABASE_DIR/${{BASE%%.*}}"
        fi
        metaeuk easy-predict \
            "${{FA[0]}}" \
            "$DB" \
            "$P" \
            tmp/ \
            $ARGS \
            --threads $TASK_CPUS || exit 1
        mv "$P.fas" "$P.codon.fas" ./*.tsv ./*.gff "$OUT/"
    """,
)

env = flyte.TaskEnvironment.from_task("metaeuk", metaeuk_easypredict_cmd.as_task())


@dataclass
class MetaeukResult:
    faa: File  # <prefix>.fas
    codon: File  # <prefix>.codon.fas
    tsv: File  # <prefix>.headersMap.tsv
    gff: File  # <prefix>.gff


async def metaeuk_easypredict(
    fasta: File, database: File | Dir, prefix: str = "sample", args: str = ""
) -> MetaeukResult:
    """Genes predicted in ``fasta`` against ``database`` (a protein FASTA or an MMseqs2 database Dir)."""
    results = await metaeuk_easypredict_cmd(
        fasta=fasta,
        database=database if isinstance(database, File) else None,
        database_dir=database if isinstance(database, Dir) else None,
        prefix=prefix,
        args=args,
    )

    async def picked(name: str) -> File:
        found = await results.get_file(name)
        if found is None:
            raise FileNotFoundError(f"metaeuk easy-predict wrote no {name}")
        return found

    return MetaeukResult(
        faa=await picked(f"{prefix}.fas"),
        codon=await picked(f"{prefix}.codon.fas"),
        tsv=await picked(f"{prefix}.headersMap.tsv"),
        gff=await picked(f"{prefix}.gff"),
    )
