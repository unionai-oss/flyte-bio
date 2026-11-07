"""prokka — rapid annotation of prokaryotic (and viral) genome assemblies.

Exposes :func:`prokka`, which annotates one assembly FASTA (plain or
gzipped; as upstream, a ``.gz`` input is decompressed first) and returns
prokka's outputs, named ``<prefix>.<ext>`` as upstream:

- ``gff`` (annotation with the sequences), ``gbk`` (GenBank), ``fna``
  (contigs), ``faa`` (proteins), ``ffn`` (feature nucleotide sequences),
  ``sqn`` (Sequin), ``fsa`` (contigs for tbl2asn), ``tbl`` (feature table);
- ``err`` (tbl2asn's validation errors), ``log``, ``txt`` (statistics) and
  ``tsv`` (feature table as TSV).

``proteins`` is a trusted protein FASTA or GenBank file to annotate from
first (``--proteins``); ``prodigal_tf`` a Prodigal training file
(``--prodigaltf``). The ``gbk``, ``sqn`` and ``log`` files record the date
and the run's paths.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

PROKKA_IMAGE = "community.wave.seqera.io/library/prokka_openjdk_parallel:f21b98bcef4c3579"
DEFAULT_RESOURCES = flyte.Resources(cpu=2, memory="8Gi")

OUTPUT_EXTENSIONS = ("gff", "gbk", "fna", "faa", "ffn", "sqn", "fsa", "tbl", "err", "log", "txt", "tsv")

prokka_cmd = shell.create(
    name="prokka",
    image=PROKKA_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "fasta": File,
        "proteins": File | None,
        "prodigal_tf": File | None,
        "prefix": str,
        "args": str,
    },
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        FA=({{inputs.fasta}})
        shopt -s nullglob; PROTEINS=({{inputs.proteins}}); TF=({{inputs.prodigal_tf}}); shopt -u nullglob
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"

        INPUT="${{FA[0]}}"
        NAME=$(basename "$INPUT")
        if [[ "$NAME" == *.gz ]]; then
            INPUT="$W/${{NAME%.gz}}"
            gunzip -c "${{FA[0]}}" > "$INPUT"
        fi
        EXTRA=()
        if [ ${{#PROTEINS[@]}} -gt 0 ]; then EXTRA+=(--proteins "${{PROTEINS[0]}}"); fi
        if [ ${{#TF[@]}} -gt 0 ]; then EXTRA+=(--prodigaltf "${{TF[0]}}"); fi

        prokka \
            $ARGS \
            --cpus $TASK_CPUS \
            --prefix "$P" \
            "${{EXTRA[@]}}" \
            "$INPUT"

        mv "$P"/* "$OUT"/
    """,
)

env = flyte.TaskEnvironment.from_task("prokka", prokka_cmd.as_task())


@dataclass
class ProkkaResult:
    gff: File
    gbk: File
    fna: File
    faa: File
    ffn: File
    sqn: File
    fsa: File
    tbl: File
    err: File
    log: File
    txt: File
    tsv: File
    results: Dir  # everything prokka wrote


async def prokka(
    fasta: File,
    proteins: File | None = None,
    prodigal_tf: File | None = None,
    prefix: str = "sample",
    args: str = "",
) -> ProkkaResult:
    """Annotate the assembly ``fasta``."""
    results = await prokka_cmd(fasta=fasta, proteins=proteins, prodigal_tf=prodigal_tf, prefix=prefix, args=args)

    async def pick(ext: str) -> File:
        found = await results.get_file(f"{prefix}.{ext}")
        if found is None:
            raise FileNotFoundError(f"prokka wrote no {prefix}.{ext}")
        return found

    files = {ext: await pick(ext) for ext in OUTPUT_EXTENSIONS}
    return ProkkaResult(**files, results=results)
