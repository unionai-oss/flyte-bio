"""vrhyme — binning of viral genomes from metagenomes with vRhyme.

Exposes :func:`vrhyme`, which bins a metagenome's viral contigs (a FASTA,
may be gzipped) by their coverage in one sample's reads (single- or
paired-end; vRhyme maps them itself) and their sequence features, and
returns, as upstream:

- the best binning's bins, the Dir ``vRhyme_best_bins_fasta`` (each bin's
  ``vRhyme_bin_<n>.fasta`` with its proteins ``.faa`` and genes ``.ffn``);
- its membership table, ``vRhyme_best_bins.<n>.membership.tsv``;
- its summary, ``vRhyme_best_bins.<n>.summary.tsv``.

As upstream, ``-t`` is the task's CPUs; options go in ``args``.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

VRHYME_IMAGE = "community.wave.seqera.io/library/vrhyme_numpy_python_scikit-learn_setuptools:02d8cf67999994aa"
# Upstream: process_high.
DEFAULT_RESOURCES = flyte.Resources(cpu=4, memory="16Gi")

vrhyme_cmd = shell.create(
    name="vrhyme",
    image=VRHYME_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"reads_1": File, "reads_2": File | None, "fasta": File, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        R1=({{inputs.reads_1}})
        shopt -s nullglob; R2=({{inputs.reads_2}}); shopt -u nullglob
        FA=({{inputs.fasta}})
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        NAME=$(basename "${{FA[0]}}")
        if [ "${{NAME##*.}}" = gz ]; then
            NAME=${{NAME%.gz}}
            gunzip -c "${{FA[0]}}" > "$NAME" || exit 1
        else
            ln -s "${{FA[0]}}" "$NAME"
        fi
        vRhyme \
            -i "$NAME" \
            -r "${{R1[0]}}" "${{R2[@]}}" \
            -o "$P" \
            -t $TASK_CPUS \
            $ARGS \
            || exit 1
        mv "$P/vRhyme_best_bins_fasta" "$OUT"/ || exit 1
        find "$P" -name 'vRhyme_best_bins.*.membership.tsv' -exec mv {{}} "$OUT"/ \;
        find "$P" -name 'vRhyme_best_bins.*.summary.tsv' -exec mv {{}} "$OUT"/ \;
        true
    """,
)

env = flyte.TaskEnvironment.from_task("vrhyme", vrhyme_cmd.as_task())


@dataclass
class VrhymeResult:
    bins: Dir  # vRhyme_best_bins_fasta: each bin's .fasta, .faa and .ffn
    membership: File  # vRhyme_best_bins.<n>.membership.tsv
    summary: File  # vRhyme_best_bins.<n>.summary.tsv


async def vrhyme(
    reads_1: File,
    fasta: File,
    reads_2: File | None = None,
    prefix: str = "sample",
    args: str = "",
) -> VrhymeResult:
    """Bin the viral contigs in ``fasta`` using one sample's reads; ``reads_2`` is None for single-end reads."""
    results = await vrhyme_cmd(reads_1=reads_1, reads_2=reads_2, fasta=fasta, prefix=prefix, args=args)
    top = await results.list_files()

    def ending(suffix: str) -> File:
        found = [f for f in top if f.name.startswith("vRhyme_best_bins.") and f.name.endswith(suffix)]
        if not found:
            raise FileNotFoundError(f"vRhyme wrote no vRhyme_best_bins.*{suffix}")
        return found[0]

    return VrhymeResult(
        bins=Dir.from_existing_remote(results.path.rstrip("/") + "/vRhyme_best_bins_fasta"),
        membership=ending(".membership.tsv"),
        summary=ending(".summary.tsv"),
    )
