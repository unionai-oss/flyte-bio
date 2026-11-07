"""mgnifam — build protein families from sequence clusters with mgnifam.

Exposes :func:`mgnifam_generatefamilies`, which runs ``mgnifam
generate_families`` on one chunk of a clustering (``clustering``, a TSV of
representative and member sequence names) against the sequences in
``fasta``. Each cluster that survives becomes a family: a seed alignment, a
profile HMM built from it, a full alignment of the sequences the HMM
recruits, and the family's reference annotation line. ``fasta_index`` is an
existing SSI index of ``fasta`` (e.g. from
:func:`flyte_bio.modules.hmmer.esl_sfetch_index`); without one, mgnifam
builds its own.

``prefix`` names the chunk (``--chunk_id``) and its files. The result Dir is
mgnifam's output directory, holding upstream's outputs under their names:

- per family ``<prefix>_<n>``: ``seed_msa/<prefix>_<n>.sto.gz``,
  ``full_msa/<prefix>_<n>.sto.gz``, ``hmm/<prefix>_<n>.hmm.gz`` and
  ``rf/<prefix>_<n>.txt``;
- for the chunk: ``<prefix>_families.tsv``, ``<prefix>_metadata.csv``,
  ``<prefix>.log`` (wall-clock timestamps), ``<prefix>_reps.fasta.gz``,
  ``<prefix>_successful.txt``, ``<prefix>_discarded.csv``,
  ``<prefix>_converged.txt`` and ``<prefix>_mgnifam_stats.json``.

Every output is optional, as upstream: a chunk whose clusters are all
discarded has no families. mgnifam's options (``--prefetch_targets``,
``--max_seq_identity``, ``--batch_size``) go in ``args``.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

MGNIFAM_IMAGE = "quay.io/biocontainers/mgnifam:4.0.0--pyhdfd78af_0"
# Upstream labels this process_high; its tests run with 4 CPUs.
DEFAULT_RESOURCES = flyte.Resources(cpu=4, memory="16Gi")

mgnifam_generatefamilies_cmd = shell.create(
    name="mgnifam_generatefamilies",
    image=MGNIFAM_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"clustering": File, "fasta": File, "fasta_index": File | None, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        CLUSTERS=({{inputs.clustering}})
        FA=({{inputs.fasta}})
        shopt -s nullglob; INDEX=({{inputs.fasta_index}}); shopt -u nullglob
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"

        # The inputs side by side under their own names, as upstream stages them.
        ln -s "${{CLUSTERS[0]}}" "$(basename "${{CLUSTERS[0]}}")"
        ln -s "${{FA[0]}}" "$(basename "${{FA[0]}}")"
        INDEXOPT=()
        if [ ${{#INDEX[@]}} -gt 0 ]; then
            ln -s "${{INDEX[0]}}" "$(basename "${{INDEX[0]}}")"
            INDEXOPT=(--fasta_index "$(basename "${{INDEX[0]}}")")
        fi
        mgnifam generate_families \
            --clusters_chunk "$(basename "${{CLUSTERS[0]}}")" \
            --fasta_file "$(basename "${{FA[0]}}")" \
            --output_dir "$P" \
            --cpus $TASK_CPUS \
            --chunk_id "$P" \
            "${{INDEXOPT[@]}}" \
            $ARGS || exit 1
        shopt -s dotglob
        mv "$P"/* "$OUT"/
    """,
)

env = flyte.TaskEnvironment.from_task("mgnifam", mgnifam_generatefamilies_cmd.as_task())


@dataclass
class MgnifamFamilies:
    seed_msa: list[File]  # seed_msa/<prefix>_<n>.sto.gz, by name
    full_msa: list[File]  # full_msa/<prefix>_<n>.sto.gz, by name
    hmm: list[File]  # hmm/<prefix>_<n>.hmm.gz, by name
    rf: list[File]  # rf/<prefix>_<n>.txt, by name
    tsv: File | None  # <prefix>_families.tsv
    csv: File | None  # <prefix>_metadata.csv
    log: File | None  # <prefix>.log
    reps_fasta: File | None  # <prefix>_reps.fasta.gz
    successful: File | None  # <prefix>_successful.txt
    discarded: File | None  # <prefix>_discarded.csv
    converged: File | None  # <prefix>_converged.txt
    stats: File | None  # <prefix>_mgnifam_stats.json
    results: Dir  # mgnifam's output directory


async def mgnifam_generatefamilies(
    clustering: File, fasta: File, fasta_index: File | None = None, prefix: str = "sample", args: str = ""
) -> MgnifamFamilies:
    """Protein families from one clustering chunk of ``fasta``'s sequences."""
    results = await mgnifam_generatefamilies_cmd(
        clustering=clustering, fasta=fasta, fasta_index=fasta_index, prefix=prefix, args=args
    )
    base = results.path.rstrip("/")
    paths = sorted([f.path[len(base) :].lstrip("/") async for f in results.walk()])

    def family_files(folder: str, suffix: str) -> list[File]:
        return [
            File.from_existing_remote(f"{base}/{p}")
            for p in paths
            if p.startswith(f"{folder}/") and "/" not in p[len(folder) + 1 :] and p.endswith(suffix)
        ]

    def single(name: str) -> File | None:
        return File.from_existing_remote(f"{base}/{name}") if name in paths else None

    return MgnifamFamilies(
        seed_msa=family_files("seed_msa", ".sto.gz"),
        full_msa=family_files("full_msa", ".sto.gz"),
        hmm=family_files("hmm", ".hmm.gz"),
        rf=family_files("rf", ".txt"),
        tsv=single(f"{prefix}_families.tsv"),
        csv=single(f"{prefix}_metadata.csv"),
        log=single(f"{prefix}.log"),
        reps_fasta=single(f"{prefix}_reps.fasta.gz"),
        successful=single(f"{prefix}_successful.txt"),
        discarded=single(f"{prefix}_discarded.csv"),
        converged=single(f"{prefix}_converged.txt"),
        stats=single(f"{prefix}_mgnifam_stats.json"),
        results=results,
    )
