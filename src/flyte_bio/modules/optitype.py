"""optitype — HLA class I genotyping with OptiType.

Exposes:

- :func:`optitype` — call the HLA-A, -B and -C alleles of one sample from its
  reads aligned to the HLA reference (one or more BAMs with their indexes).

Outputs keep upstream's names, under one output Dir:

- the calls: ``<prefix>_result.tsv``;
- the coverage plot: ``<prefix>_coverage_plot.pdf``.

As upstream, OptiType gets a ``config.ini`` that runs razers3 on the task's
CPUs, solves the integer program with ``solver`` (GLPK, single-threaded, by
default; CBC uses the task's CPUs), deletes its intermediate BAMs, and sets
``unpaired_weight`` and ``use_discordant``.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

OPTITYPE_IMAGE = "community.wave.seqera.io/library/optitype_coincbc:722e52fd4495c161"
# Upstream's process_medium asks for 6 CPUs and 36 GB; HLA-region BAMs need far less.
DEFAULT_RESOURCES = flyte.Resources(cpu=2, memory="8Gi")

optitype_cmd = shell.create(
    name="optitype",
    image=OPTITYPE_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "bams": list[File],
        "bais": list[File],
        "seq_type": str,
        "prefix": str,
        "solver": str,
        "unpaired_weight": str,
        "use_discordant": bool,
        "args": str,
    },
    defaults={"solver": "glpk", "unpaired_weight": "0", "use_discordant": False, "args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        BAMS=({{inputs.bams}})
        shopt -s nullglob; BAIS=({{inputs.bais}}); shopt -u nullglob
        P={{inputs.prefix}}
        SOLVER={{inputs.solver}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"

        INPUTS=()
        for f in "${{BAMS[@]}}"; do
            ln -s "$f" "$(basename "$f")"
            INPUTS+=(--input "$(basename "$f")")
        done
        for f in "${{BAIS[@]}}"; do ln -s "$f" "$(basename "$f")"; done

        SOLVER_THREADS=$TASK_CPUS
        if [ "$SOLVER" = glpk ]; then SOLVER_THREADS=1; fi
        {{
            echo "[mapping]"
            echo "razers3=razers3"
            echo "threads=$TASK_CPUS"
            echo "[ilp]"
            echo "solver=$SOLVER"
            echo "threads=$SOLVER_THREADS"
            echo "[behavior]"
            echo "deletebam=true"
            echo "unpaired_weight={{inputs.unpaired_weight}}"
            echo "use_discordant={{inputs.use_discordant}}"
        }} > config.ini

        optitype run --{{inputs.seq_type}} --config config.ini $ARGS --prefix "$P" --outdir "$OUT" "${{INPUTS[@]}}" \
            || exit 1
    """,
)

env = flyte.TaskEnvironment.from_task("optitype", optitype_cmd.as_task())


@dataclass
class OptitypeResult:
    hla_type: File  # <prefix>_result.tsv
    coverage_plot: File  # <prefix>_coverage_plot.pdf
    results: Dir  # everything the task wrote


async def optitype(
    bams: list[File],
    bais: list[File],
    prefix: str,
    seq_type: str = "dna",
    solver: str = "glpk",
    unpaired_weight: float = 0,
    use_discordant: bool = False,
    args: str = "",
) -> OptitypeResult:
    """Type the HLA class I alleles of one sample's ``bams`` (``seq_type`` ``dna`` or ``rna``).

    ``bais`` are the BAMs' indexes. ``solver``, ``unpaired_weight`` and
    ``use_discordant`` are upstream's ``ext.args2`` settings; ``args`` go to
    ``optitype run`` (upstream's tests use ``--enumerate 1 --beta 0.009``).
    """
    if seq_type not in ("dna", "rna"):
        raise ValueError(f"seq_type must be 'dna' or 'rna', got {seq_type!r}")
    weight = int(unpaired_weight) if float(unpaired_weight).is_integer() else unpaired_weight
    results = await optitype_cmd(
        bams=bams,
        bais=bais,
        seq_type=seq_type,
        prefix=prefix,
        solver=solver,
        unpaired_weight=str(weight),
        use_discordant=use_discordant,
        args=args,
    )

    async def pick(name: str) -> File:
        file = await results.get_file(name)
        if file is None:
            raise FileNotFoundError(f"optitype wrote no {name}")
        return file

    return OptitypeResult(
        hla_type=await pick(f"{prefix}_result.tsv"),
        coverage_plot=await pick(f"{prefix}_coverage_plot.pdf"),
        results=results,
    )
