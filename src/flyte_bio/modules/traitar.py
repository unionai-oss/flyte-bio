"""traitar — Traitar, prediction of microbial phenotypic traits from Pfam annotations.

Exposes:

- :func:`traitar_pfamget` — download the Pfam HMMs Traitar's models were
  trained on (``traitar pfam``): a Dir holding ``Pfam-A.hmm.gz`` and
  ``Pfam-A.hmm``.
- :func:`traitar_run` — predict phenotypes for one sample (``traitar
  phenotype``) from nucleotides (``from_nucleotides``: genes are called with
  Prodigal first), proteins (``from_genes``) or a Pfam annotation summary
  (``from_annotation_summary``). A gzipped input is decompressed first. As
  upstream, the sample file lists it under the sample name ``prefix``.

  The predictions keep upstream's names: ``predictions_majority-vote_combined.txt``,
  ``predictions_single-votes_combined.txt`` and ``predictions_flat_*.txt``
  (from ``phenotype_prediction/``), and ``predictions_raw.txt``,
  ``predictions_majority-vote.txt``, ``predictions_single-votes.txt`` and
  ``predictions_conservative-vote.txt``. Traitar's whole output folder
  (gene prediction, Pfam annotation, heatmaps) is archived as
  ``<prefix>.tar.gz``.

As upstream, Traitar gets the task's CPUs as parallel jobs (``-c``).
"""

import os
from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

TRAITAR_IMAGE = "community.wave.seqera.io/library/hmmer_prodigal_pandas_parallel_pruned:ccae2eabc2a54ac8"
# Pfam-A is about 1.5 GiB unpacked, and the download is kept beside it.
PFAMGET_RESOURCES = flyte.Resources(cpu=1, memory="4Gi", disk="20Gi")
# Upstream's process_medium, scaled down.
RUN_RESOURCES = flyte.Resources(cpu=8, memory="16Gi", disk="20Gi")
INPUT_TYPES = ("from_genes", "from_nucleotides", "from_annotation_summary")

traitar_pfamget_cmd = shell.create(
    name="traitar_pfamget",
    image=TRAITAR_IMAGE,
    resources=PFAMGET_RESOURCES,
    inputs={"args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        ARGS={inputs.args}
        OUT={outputs.results}
        W=$(mktemp -d); cd "$W"
        export HOME="$W"
        traitar pfam $ARGS pfam_data || exit 1
        mv pfam_data/* "$OUT"/
    """,
)

traitar_run_cmd = shell.create(
    name="traitar_run",
    image=TRAITAR_IMAGE,
    resources=RUN_RESOURCES,
    inputs={"fasta": File, "input_type": str, "pfam_db": Dir, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        FA=({{inputs.fasta}})
        INPUT_TYPE={{inputs.input_type}}
        PFAM={{inputs.pfam_db}}
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        export HOME="$W"
        mkdir -p input_dir
        NAME=$(basename "${{FA[0]}}")
        if [[ "$NAME" == *.gz ]]; then
            NAME="${{NAME%.gz}}"
            gunzip -c "${{FA[0]}}" > "input_dir/$NAME"
        else
            ln -s "${{FA[0]}}" "input_dir/$NAME"
        fi
        printf 'sample_file_name\tsample_name\n%s\t%s\n' "$NAME" "$P" > samples.txt
        traitar phenotype "$PFAM" input_dir samples.txt "$INPUT_TYPE" "$P" -c $TASK_CPUS --overwrite $ARGS || exit 1
        shopt -s nullglob
        cp "$P"/phenotype_prediction/predictions_*.txt "$P"/predictions_*.txt "$OUT"/
        tar -czf "$OUT/$P.tar.gz" "$P"
    """,
)

env = flyte.TaskEnvironment.from_task("traitar", traitar_pfamget_cmd.as_task(), traitar_run_cmd.as_task())


async def traitar_pfamget(args: str = "") -> Dir:
    """The Pfam HMMs Traitar uses: a Dir of ``Pfam-A.hmm.gz`` and ``Pfam-A.hmm``."""
    return await traitar_pfamget_cmd(args=args)


@dataclass
class TraitarResult:
    predictions_combined: File  # predictions_majority-vote_combined.txt
    predictions_single_votes: File  # predictions_single-votes_combined.txt
    predictions_flat: list[File]  # predictions_flat_*.txt
    predictions_raw: list[File]  # predictions_{raw,majority-vote,single-votes,conservative-vote}.txt
    archive: File  # <prefix>.tar.gz, Traitar's whole output folder
    results: Dir


async def traitar_run(
    fasta: File, input_type: str, pfam_db: Dir, prefix: str = "sample", args: str = ""
) -> TraitarResult:
    """Phenotypes predicted for one sample; ``pfam_db`` comes from :func:`traitar_pfamget`."""
    if input_type not in INPUT_TYPES:
        raise ValueError(f"input_type must be one of {', '.join(INPUT_TYPES)}, not {input_type!r}")
    results = await traitar_run_cmd(fasta=fasta, input_type=input_type, pfam_db=pfam_db, prefix=prefix, args=args)

    async def pick(name: str) -> File:
        found = await results.get_file(name)
        if found is None:
            raise FileNotFoundError(f"traitar wrote no {name}")
        return found

    names = sorted(os.path.basename(f.path) for f in await results.list_files())
    combined = ("predictions_majority-vote_combined.txt", "predictions_single-votes_combined.txt")
    flat = [n for n in names if n.startswith("predictions_flat_")]
    raw = [n for n in names if n.startswith("predictions_") and n not in combined and n not in flat]
    return TraitarResult(
        predictions_combined=await pick(combined[0]),
        predictions_single_votes=await pick(combined[1]),
        predictions_flat=[await pick(n) for n in flat],
        predictions_raw=[await pick(n) for n in raw],
        archive=await pick(f"{prefix}.tar.gz"),
        results=results,
    )
