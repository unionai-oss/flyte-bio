"""fibertoolsrs — fibertools-rs (``ft``): m6A and nucleosome calling on single-molecule chromatin fibers.

Exposes:

- :func:`fibertoolsrs_predictm6a` — ``ft predict-m6a``: predict m6A
  (adenine) methylation from PacBio kinetics, writing ``<prefix>.bam``.
- :func:`fibertoolsrs_addnucleosomes` — ``ft add-nucleosomes``: call
  nucleosomes and methyltransferase-sensitive patches (MSPs) from the m6A
  calls, writing ``<prefix>.bam``.
- :func:`fibertoolsrs_extract` — ``ft extract``: the fibers' ``m6a``,
  ``cpg``, ``msp``, ``nuc`` or ``all`` features as a BED, ``<prefix>.bed.gz``.

Upstream's prefixes default to ``<id>_m6a`` and ``<id>_nucleosomes`` for the
BAMs and ``<id>`` for the BED. As upstream, ``ft`` gets the task's CPUs as
threads and the input is linked into a working directory under its own name,
so the BAM header's ``@PG`` command line names it as upstream's does.
``args`` go to the ``ft`` subcommand.
"""

import os

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

FIBERTOOLSRS_IMAGE = "quay.io/biocontainers/fibertools-rs:0.8.2--h3b373d1_0"
# Upstream labels every tool process_low (2 CPUs). Its test snapshots are of
# 2-thread runs: the output BAMs' @PG command lines record --threads.
DEFAULT_RESOURCES = flyte.Resources(cpu=2, memory="8Gi")

EXTRACT_TYPES = ("m6a", "cpg", "msp", "nuc", "all")


def bam_to_bam(name: str, subcommand: str):
    return shell.create(
        name=name,
        image=FIBERTOOLSRS_IMAGE,
        resources=DEFAULT_RESOURCES,
        inputs={"bam": File, "prefix": str, "args": str},
        defaults={"args": ""},
        outputs={"results": Dir},
        script=rf"""
            {TASK_LIMITS}
            BAM=({{inputs.bam}})
            P={{inputs.prefix}}
            ARGS={{inputs.args}}
            OUT={{outputs.results}}
            W=$(mktemp -d); cd "$W"
            NAME=$(basename "${{BAM[0]}}"); ln -s "${{BAM[0]}}" "$NAME"
            # The output's @PG line records this command line, thread count included.
            ft {subcommand} $ARGS --threads $TASK_CPUS "$NAME" "$P.bam" || exit 1
            mv "$P.bam" "$OUT"/
        """,
    )


fibertoolsrs_predictm6a_cmd = bam_to_bam("fibertoolsrs_predictm6a", "predict-m6a")
fibertoolsrs_addnucleosomes_cmd = bam_to_bam("fibertoolsrs_addnucleosomes", "add-nucleosomes")

fibertoolsrs_extract_cmd = shell.create(
    name="fibertoolsrs_extract",
    image=FIBERTOOLSRS_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"bam": File, "extract_type": str, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        BAM=({{inputs.bam}})
        TYPE={{inputs.extract_type}}
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        NAME=$(basename "${{BAM[0]}}"); ln -s "${{BAM[0]}}" "$NAME"
        ft extract $ARGS --threads $TASK_CPUS "$NAME" "--$TYPE" "$P.bed.gz" || exit 1
        mv "$P.bed.gz" "$OUT"/
    """,
)

env = flyte.TaskEnvironment.from_task(
    "fibertoolsrs",
    fibertoolsrs_predictm6a_cmd.as_task(),
    fibertoolsrs_addnucleosomes_cmd.as_task(),
    fibertoolsrs_extract_cmd.as_task(),
)


def file_name(file: File) -> str:
    return file.name or os.path.basename(file.path)


async def pick(results: Dir, name: str, tool: str) -> File:
    found = await results.get_file(name)
    if found is None:
        raise FileNotFoundError(f"ft {tool} wrote no {name}")
    return found


async def fibertoolsrs_predictm6a(bam: File, prefix: str = "sample_m6a", args: str = "") -> File:
    """``bam`` with m6A calls predicted from its kinetics: ``<prefix>.bam``."""
    if file_name(bam) == f"{prefix}.bam":
        raise ValueError(f"input and output are both {prefix}.bam: set another prefix")
    results = await fibertoolsrs_predictm6a_cmd(bam=bam, prefix=prefix, args=args)
    return await pick(results, f"{prefix}.bam", "predict-m6a")


async def fibertoolsrs_addnucleosomes(bam: File, prefix: str = "sample_nucleosomes", args: str = "") -> File:
    """``bam`` with nucleosome and MSP calls: ``<prefix>.bam``."""
    if file_name(bam) == f"{prefix}.bam":
        raise ValueError(f"input and output are both {prefix}.bam: set another prefix")
    results = await fibertoolsrs_addnucleosomes_cmd(bam=bam, prefix=prefix, args=args)
    return await pick(results, f"{prefix}.bam", "add-nucleosomes")


async def fibertoolsrs_extract(bam: File, extract_type: str, prefix: str = "sample", args: str = "") -> File:
    """``bam``'s ``extract_type`` features (one of m6a, cpg, msp, nuc, all) as ``<prefix>.bed.gz``."""
    if extract_type not in EXTRACT_TYPES:
        raise ValueError(f"Invalid extract_type='{extract_type}'. Allowed: {', '.join(EXTRACT_TYPES)}")
    results = await fibertoolsrs_extract_cmd(bam=bam, extract_type=extract_type, prefix=prefix, args=args)
    return await pick(results, f"{prefix}.bed.gz", "extract")
