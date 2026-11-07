"""spades — SPAdes de novo genome (and RNA, viral, metagenome) assembly.

Exposes :func:`spades`, which assembles one sample from Illumina reads
(single-end, or an R1/R2 pair), optionally with PacBio or Nanopore long reads,
or from a SPAdes dataset YAML (``--dataset``), with optional custom HMMs
(``--custom-hmms``). Pass the assembly mode (``--rnaviral``, ``--meta``,
``--isolate``, …) in ``args``.

As upstream, every input is staged into the working directory under its
original name, so a dataset YAML can name the reads by their basenames. The
outputs keep upstream's names, each gzipped with ``gzip -n``:

- ``<prefix>.scaffolds.fa.gz``, ``<prefix>.contigs.fa.gz``,
  ``<prefix>.transcripts.fa.gz``, ``<prefix>.gene_clusters.fa.gz``,
  ``<prefix>.assembly.gfa.gz``: whichever ones the mode writes;
- ``<prefix>.spades.log`` and, when SPAdes warned, ``<prefix>.warnings.log``.

``--threads`` and ``--memory`` (whole GiB) follow the task's CPUs and memory.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

SPADES_IMAGE = "community.wave.seqera.io/library/spades:4.1.0--77799c52e1d1054a"

# Upstream's label asks for 12 CPUs and 72 GB. These defaults fit viral and
# small bacterial assemblies; raise them with
# spades_cmd.as_task().override(resources=...) for larger genomes or
# metagenomes. The 4 CPUs match upstream's tests.
DEFAULT_RESOURCES = flyte.Resources(cpu=4, memory="16Gi")

spades_cmd = shell.create(
    name="spades",
    image=SPADES_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "reads_1": File | None,
        "reads_2": File | None,
        "pacbio": File | None,
        "nanopore": File | None,
        "dataset_yml": File | None,
        "hmm": File | None,
        "prefix": str,
        "args": str,
    },
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        shopt -s nullglob
        ILLUMINA=({{inputs.reads_1}} {{inputs.reads_2}}); PACBIO=({{inputs.pacbio}}); NANOPORE=({{inputs.nanopore}})
        YML=({{inputs.dataset_yml}}); HMM=({{inputs.hmm}})
        shopt -u nullglob
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"

        name() {{ basename "$1"; }}
        for F in "${{ILLUMINA[@]}}" "${{PACBIO[@]}}" "${{NANOPORE[@]}}" "${{YML[@]}}" "${{HMM[@]}}"; do
            ln -s "$F" "$(name "$F")"
        done

        READS=()
        if [ ${{#YML[@]}} -gt 0 ]; then
            READS=(--dataset "$(name "${{YML[0]}}")")
        else
            if [ ${{#ILLUMINA[@]}} -eq 1 ]; then
                READS+=(-s "$(name "${{ILLUMINA[0]}}")")
            elif [ ${{#ILLUMINA[@]}} -eq 2 ]; then
                READS+=(-1 "$(name "${{ILLUMINA[0]}}")" -2 "$(name "${{ILLUMINA[1]}}")")
            fi
            if [ ${{#PACBIO[@]}} -gt 0 ]; then READS+=(--pacbio "$(name "${{PACBIO[0]}}")"); fi
            if [ ${{#NANOPORE[@]}} -gt 0 ]; then READS+=(--nanopore "$(name "${{NANOPORE[0]}}")"); fi
        fi
        HMMS=()
        if [ ${{#HMM[@]}} -gt 0 ]; then HMMS=(--custom-hmms "$(name "${{HMM[0]}}")"); fi

        spades.py \
            $ARGS \
            --threads $TASK_CPUS \
            --memory $(( TASK_MEMORY_MB / 1024 )) \
            "${{HMMS[@]}}" \
            "${{READS[@]}}" \
            -o ./
        mv spades.log "$OUT/$P.spades.log"

        for PAIR in scaffolds.fasta:scaffolds.fa contigs.fasta:contigs.fa transcripts.fasta:transcripts.fa \
                assembly_graph_with_scaffolds.gfa:assembly.gfa gene_clusters.fasta:gene_clusters.fa; do
            SRC=${{PAIR%%:*}}; DEST=$P.${{PAIR#*:}}
            if [ -f "$SRC" ]; then
                mv "$SRC" "$DEST"
                gzip -n "$DEST"
                mv "$DEST.gz" "$OUT"/
            fi
        done
        if [ -f warnings.log ]; then mv warnings.log "$OUT/$P.warnings.log"; fi
        true
    """,
)

env = flyte.TaskEnvironment.from_task("spades", spades_cmd.as_task())


@dataclass
class SpadesResult:
    scaffolds: File | None
    contigs: File | None
    transcripts: File | None
    gene_clusters: File | None
    gfa: File | None  # the assembly graph with scaffolds
    warnings: File | None  # only when SPAdes warned
    log: File
    results: Dir  # everything above


async def spades(
    illumina: list[File] | None = None,
    pacbio: File | None = None,
    nanopore: File | None = None,
    dataset_yml: File | None = None,
    hmm: File | None = None,
    prefix: str = "sample",
    args: str = "",
) -> SpadesResult:
    """Assemble one sample.

    ``illumina`` is one file (single-end, ``-s``) or an R1/R2 pair. With
    ``dataset_yml`` SPAdes reads only the YAML, which names the reads (also
    passed here, so they're staged next to it) by their basenames.
    """
    illumina = illumina or []
    if len(illumina) > 2:
        raise ValueError(f"spades takes one or two Illumina files, got {len(illumina)}")
    if not (illumina or pacbio or nanopore or dataset_yml):
        raise ValueError("spades needs reads or a dataset YAML")
    results = await spades_cmd(
        reads_1=illumina[0] if illumina else None,
        reads_2=illumina[1] if len(illumina) == 2 else None,
        pacbio=pacbio,
        nanopore=nanopore,
        dataset_yml=dataset_yml,
        hmm=hmm,
        prefix=prefix,
        args=args,
    )

    return await read_spades_results(results, prefix)


async def read_spades_results(results: Dir, prefix: str) -> SpadesResult:
    """Pick the outputs out of ``spades_cmd``'s results Dir (e.g. from an overridden task)."""

    async def find(suffix: str) -> File | None:
        return await results.get_file(f"{prefix}.{suffix}")

    log = await find("spades.log")
    if log is None:
        raise FileNotFoundError(f"spades wrote no {prefix}.spades.log")
    return SpadesResult(
        scaffolds=await find("scaffolds.fa.gz"),
        contigs=await find("contigs.fa.gz"),
        transcripts=await find("transcripts.fa.gz"),
        gene_clusters=await find("gene_clusters.fa.gz"),
        gfa=await find("assembly.gfa.gz"),
        warnings=await find("warnings.log"),
        log=log,
        results=results,
    )
