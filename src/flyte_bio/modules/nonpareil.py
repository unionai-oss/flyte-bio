"""nonpareil — metagenome coverage and sequence diversity estimates.

- :func:`nonpareil_nonpareil` estimates one sample's coverage from its reads
  (``nonpareil``; ``format`` ``fastq`` or ``fasta``, ``mode`` ``kmer`` or
  ``alignment``): the redundancy summary ``<prefix>.npo`` (Nonpareil's main
  output), the redundancy values ``<prefix>.npa``, the mates' distribution
  ``<prefix>.npc`` and the log ``<prefix>.npl``. ``-t`` is the task's CPUs
  and ``-R`` its memory in MiB; options such as ``-X`` go in ``args``.
- :func:`nonpareil_curve` plots one sample's Nonpareil curve
  (``Nonpareil.curve`` in R), :func:`nonpareil_set` several samples' curves
  together (``Nonpareil.set``): ``<prefix>.png``. ``r_args`` are extra R
  arguments to the plotting function, e.g. ``col = 'green'``.
- :func:`nonpareil_nonpareilcurvesr` summarises several samples' curves
  (``NonpareilCurves.R``): ``<prefix>.json``, ``<prefix>.tsv``,
  ``<prefix>.csv`` and the plot ``<prefix>.pdf``.

The curves name each sample after its ``.npo`` file.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

NONPAREIL_IMAGE = "quay.io/biocontainers/nonpareil:3.5.5--r43hdcf5f25_0"
# Upstream labels nonpareil process_medium and the plots process_single.
DEFAULT_RESOURCES = flyte.Resources(cpu=4, memory="8Gi")
PLOT_RESOURCES = flyte.Resources(cpu=1, memory="2Gi")

nonpareil_nonpareil_cmd = shell.create(
    name="nonpareil_nonpareil",
    image=NONPAREIL_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"reads": File, "format": str, "mode": str, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        READS=({{inputs.reads}})
        FORMAT={{inputs.format}}
        MODE={{inputs.mode}}
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        nonpareil \
            -s "${{READS[0]}}" \
            -f "$FORMAT" \
            -T "$MODE" \
            -t $TASK_CPUS \
            -R $TASK_MEMORY_MB \
            -b "$P" \
            $ARGS \
            || exit 1
        mv "$P".npa "$P".npc "$P".npl "$P".npo "$OUT"/
    """,
)

# Plot the curves of the .npo files in the working directory with an R function.
PLOT_SCRIPT = r"""
        NPO=({inputs.npo})
        P={inputs.prefix}
        R_ARGS={inputs.r_args}
        OUT={outputs.results}
        W=$(mktemp -d); cd "$W"
        for F in "${NPO[@]}"; do ln -s "$F" "$(basename "$F")"; done
        FIRST=$(basename "${NPO[0]}")
        EXTRA=""
        if [ -n "$R_ARGS" ]; then EXTRA=", $R_ARGS"; fi
        Rscript - <<EOF || exit 1
library(Nonpareil)
png(file='$P.png')
PLOT_CALL
dev.off()
EOF
        mv "$P.png" "$OUT"/
"""

nonpareil_curve_cmd = shell.create(
    name="nonpareil_curve",
    image=NONPAREIL_IMAGE,
    resources=PLOT_RESOURCES,
    inputs={"npo": File, "prefix": str, "r_args": str},
    defaults={"r_args": ""},
    outputs={"results": Dir},
    script=PLOT_SCRIPT.replace("PLOT_CALL", "Nonpareil.curve('$FIRST'$EXTRA)"),
)

nonpareil_set_cmd = shell.create(
    name="nonpareil_set",
    image=NONPAREIL_IMAGE,
    resources=PLOT_RESOURCES,
    inputs={"npo": list[File], "prefix": str, "r_args": str},
    defaults={"r_args": ""},
    outputs={"results": Dir},
    script=PLOT_SCRIPT.replace("PLOT_CALL", "Nonpareil.set(list.files(pattern='*.npo')$EXTRA)"),
)

nonpareil_nonpareilcurvesr_cmd = shell.create(
    name="nonpareil_nonpareilcurvesr",
    image=NONPAREIL_IMAGE,
    resources=PLOT_RESOURCES,
    inputs={"npo": list[File], "npo_order": str, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        NPO=({inputs.npo})
        # Staged list inputs expand in name order; restore the order they were given in.
        read -ra NPO_NAMES <<< {inputs.npo_order}; NPO_DIR=$(dirname "${NPO[0]}")
        NPO=(); for N in "${NPO_NAMES[@]}"; do NPO+=("$NPO_DIR/$N"); done
        P={inputs.prefix}
        ARGS={inputs.args}
        OUT={outputs.results}
        W=$(mktemp -d); cd "$W"
        NAMES=()
        for F in "${NPO[@]}"; do NAMES+=("$(basename "$F")"); ln -s "$F" "$(basename "$F")"; done
        NonpareilCurves.R \
            $ARGS \
            --json "$P.json" \
            --tsv "$P.tsv" \
            --csv "$P.csv" \
            --pdf "$P.pdf" \
            "${NAMES[@]}" \
            || exit 1
        shopt -s nullglob
        mv ./*.json ./*.tsv ./*.csv ./*.pdf "$OUT"/
    """,
)

env = flyte.TaskEnvironment.from_task(
    "nonpareil",
    nonpareil_nonpareil_cmd.as_task(),
    nonpareil_curve_cmd.as_task(),
    nonpareil_set_cmd.as_task(),
    nonpareil_nonpareilcurvesr_cmd.as_task(),
)


@dataclass
class NonpareilResult:
    npo: File  # <prefix>.npo, the redundancy summary
    npa: File  # <prefix>.npa, the redundancy values
    npc: File  # <prefix>.npc, the mates' distribution
    npl: File  # <prefix>.npl, the log
    results: Dir  # all of the above


async def nonpareil_nonpareil(
    reads: File, format: str = "fastq", mode: str = "kmer", prefix: str = "sample", args: str = ""
) -> NonpareilResult:
    """Estimate one sample's coverage from its ``reads`` (``format`` fastq or fasta; ``mode`` kmer or alignment)."""
    results = await nonpareil_nonpareil_cmd(reads=reads, format=format, mode=mode, prefix=prefix, args=args)

    async def pick(suffix: str) -> File:
        found = await results.get_file(f"{prefix}.{suffix}")
        if found is None:
            raise FileNotFoundError(f"nonpareil wrote no {prefix}.{suffix}")
        return found

    return NonpareilResult(
        npo=await pick("npo"), npa=await pick("npa"), npc=await pick("npc"), npl=await pick("npl"), results=results
    )


async def nonpareil_curve(npo: File, prefix: str = "sample", r_args: str = "") -> File:
    """One sample's Nonpareil curve: ``<prefix>.png``."""
    results = await nonpareil_curve_cmd(npo=npo, prefix=prefix, r_args=r_args)
    png = await results.get_file(f"{prefix}.png")
    if png is None:
        raise FileNotFoundError(f"Nonpareil.curve wrote no {prefix}.png")
    return png


async def nonpareil_set(npo: list[File], prefix: str = "sample", r_args: str = "") -> File:
    """Several samples' Nonpareil curves in one plot: ``<prefix>.png``."""
    results = await nonpareil_set_cmd(npo=npo, prefix=prefix, r_args=r_args)
    png = await results.get_file(f"{prefix}.png")
    if png is None:
        raise FileNotFoundError(f"Nonpareil.set wrote no {prefix}.png")
    return png


@dataclass
class NonpareilCurvesResult:
    json: File  # <prefix>.json
    tsv: File  # <prefix>.tsv
    csv: File  # <prefix>.csv
    pdf: File  # <prefix>.pdf
    results: Dir  # all of the above


async def nonpareil_nonpareilcurvesr(npo: list[File], prefix: str = "sample", args: str = "") -> NonpareilCurvesResult:
    """Summaries of several samples' Nonpareil curves (``NonpareilCurves.R``)."""
    results = await nonpareil_nonpareilcurvesr_cmd(
        npo=npo, npo_order=" ".join(f.name for f in npo), prefix=prefix, args=args
    )

    async def pick(suffix: str) -> File:
        found = await results.get_file(f"{prefix}.{suffix}")
        if found is None:
            raise FileNotFoundError(f"NonpareilCurves.R wrote no {prefix}.{suffix}")
        return found

    return NonpareilCurvesResult(
        json=await pick("json"), tsv=await pick("tsv"), csv=await pick("csv"), pdf=await pick("pdf"), results=results
    )
