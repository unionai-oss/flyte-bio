"""orthofinder — orthogroups, orthologues, gene trees and the rooted species tree with OrthoFinder.

Exposes :func:`orthofinder`, which runs OrthoFinder on one proteome FASTA per
species (each species is named after its file, so the names must differ) or,
with ``prior_run``, adds the species to a previous run (``-b``). The prior
run is a ``.tar.gz`` of its ``WorkingDirectory``; as upstream it is extracted
(under its single top-level directory's name, if it has one) and OrthoFinder
writes the new results inside it.

The results directory, ``Results_<prefix>``, holds thousands of files on real
data, so it is returned as ``Results_<prefix>.tar.gz``, beside the pieces
upstream emits on their own:

- ``Orthogroups.tsv``: the orthogroups;
- ``SpeciesTree_rooted_node_labels.txt``: the rooted species tree;
- ``WorkingDirectory.tar.gz``: the working directory, to extend the run later;
- ``Single_Copy_Orthologue_Sequences.tar.gz``: when OrthoFinder wrote that
  directory.

``-t`` and ``-a`` follow the task's CPUs.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

ORTHOFINDER_IMAGE = "quay.io/biocontainers/orthofinder:3.1.3--hdfd78af_0"
# Upstream's label is process_high. The 4 CPUs match upstream's tests (the
# gene trees' tie-breaks can follow the thread count); raise them with
# orthofinder_cmd.override(resources=...) for many or large proteomes.
DEFAULT_RESOURCES = flyte.Resources(cpu=4, memory="16Gi")

orthofinder_cmd = shell.create(
    name="orthofinder",
    image=ORTHOFINDER_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"fastas": list[File], "prior_run": File | None, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        FASTAS=({{inputs.fastas}})
        shopt -s nullglob; PRIOR=({{inputs.prior_run}}); shopt -u nullglob
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"

        mkdir input
        for F in "${{FASTAS[@]}}"; do cp "$F" input/; done

        INCLUDE=()
        RESULTS=input/OrthoFinder/Results_$P
        if [ ${{#PRIOR[@]}} -gt 0 ]; then
            mkdir prior
            tar -xzf "${{PRIOR[0]}}" -C prior || exit 1
            TOP=(prior/*)
            if [ ${{#TOP[@]}} -eq 1 ] && [ -d "${{TOP[0]}}" ]; then
                PRIOR_DIR=$(basename "${{TOP[0]}}")
                mv "${{TOP[0]}}" "$PRIOR_DIR"
            else
                PRIOR_DIR=prior_run
                mv prior "$PRIOR_DIR"
            fi
            INCLUDE=(-b "$PRIOR_DIR")
            RESULTS=$PRIOR_DIR/OrthoFinder/Results_$P
        fi

        orthofinder -t $TASK_CPUS -a $TASK_CPUS -f input -n "$P" "${{INCLUDE[@]}}" $ARGS || exit 1

        if [ ! -d "$RESULTS" ]; then echo "OrthoFinder wrote no $RESULTS" >&2; exit 1; fi
        cp "$RESULTS/Orthogroups/Orthogroups.tsv" "$OUT"/ || exit 1
        cp "$RESULTS/Species_Tree/SpeciesTree_rooted_node_labels.txt" "$OUT"/ || exit 1
        tar -czf "$OUT/WorkingDirectory.tar.gz" -C "$RESULTS" WorkingDirectory || exit 1
        if [ -d "$RESULTS/Single_Copy_Orthologue_Sequences" ]; then
            tar -czf "$OUT/Single_Copy_Orthologue_Sequences.tar.gz" -C "$RESULTS" Single_Copy_Orthologue_Sequences
        fi
        tar -czf "$OUT/Results_$P.tar.gz" -C "$(dirname "$RESULTS")" "Results_$P" || exit 1
    """,
)

env = flyte.TaskEnvironment.from_task("orthofinder", orthofinder_cmd.as_task())


@dataclass
class OrthofinderResult:
    orthofinder: File  # Results_<prefix>.tar.gz: the whole results directory
    working: File  # WorkingDirectory.tar.gz
    single_copy_seqs: File | None  # Single_Copy_Orthologue_Sequences.tar.gz, when written
    orthogroups: File  # Orthogroups.tsv
    sptree: File  # SpeciesTree_rooted_node_labels.txt


async def orthofinder(
    fastas: list[File], prior_run: File | None = None, prefix: str = "sample", args: str = ""
) -> OrthofinderResult:
    """Orthology of the proteomes ``fastas`` (one per species), optionally added to ``prior_run``."""
    names = [f.path.rstrip("/").rsplit("/", 1)[-1] for f in fastas]
    if len(set(names)) != len(names):
        raise ValueError(f"orthofinder names species by file name, so they must differ: {names}")
    results = await orthofinder_cmd(fastas=fastas, prior_run=prior_run, prefix=prefix, args=args)

    async def pick(name: str) -> File:
        found = await results.get_file(name)
        if found is None:
            raise FileNotFoundError(f"orthofinder wrote no {name}")
        return found

    return OrthofinderResult(
        orthofinder=await pick(f"Results_{prefix}.tar.gz"),
        working=await pick("WorkingDirectory.tar.gz"),
        single_copy_seqs=await results.get_file("Single_Copy_Orthologue_Sequences.tar.gz"),
        orthogroups=await pick("Orthogroups.tsv"),
        sptree=await pick("SpeciesTree_rooted_node_labels.txt"),
    )
