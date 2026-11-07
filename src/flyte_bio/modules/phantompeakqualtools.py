"""phantompeakqualtools — ChIP-seq strand cross-correlation quality metrics.

Exposes :func:`phantompeakqualtools`, which runs the image's ``run_spp.R``
(SPP's cross-correlation analysis) on one BAM and returns, under upstream's
names:

- the metrics line (fragment length estimates, NSC, RSC, quality tag),
  ``<prefix>.spp.out``;
- the cross-correlation plot, ``<prefix>.spp.pdf``;
- the R session data, ``<prefix>.spp.Rdata``.

As upstream, ``args`` go to ``Rscript`` and ``args2`` to ``run_spp.R`` (e.g.
``-p=<threads>`` or ``-x=<exclusion range>``). The BAM is linked into the
working directory under its own name, which the metrics line records.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

PHANTOMPEAKQUALTOOLS_IMAGE = "community.wave.seqera.io/library/phantompeakqualtools:1.2.2--f8026fe2526a5e18"
DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory="8Gi")

phantompeakqualtools_cmd = shell.create(
    name="phantompeakqualtools",
    image=PHANTOMPEAKQUALTOOLS_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"bam": File, "prefix": str, "args": str, "args2": str},
    defaults={"args": "", "args2": ""},
    outputs={"results": Dir},
    script=r"""
        BAM=({inputs.bam})
        P={inputs.prefix}
        ARGS={inputs.args}
        ARGS2={inputs.args2}
        OUT={outputs.results}
        W=$(mktemp -d); cd "$W"
        NAME=$(basename "${BAM[0]}")
        ln -s "${BAM[0]}" "$NAME"
        RUN_SPP=$(which run_spp.R)
        Rscript $ARGS -e "library(caTools); source(\"$RUN_SPP\")" -c="$NAME" \
            -savp="$P.spp.pdf" -savd="$P.spp.Rdata" -out="$P.spp.out" $ARGS2
        mv "$P.spp.pdf" "$P.spp.Rdata" "$P.spp.out" "$OUT"/
    """,
)

env = flyte.TaskEnvironment.from_task("phantompeakqualtools", phantompeakqualtools_cmd.as_task())


@dataclass
class PhantompeakqualtoolsResult:
    spp: File  # <prefix>.spp.out: the metrics line
    pdf: File  # <prefix>.spp.pdf: the cross-correlation plot
    rdata: File  # <prefix>.spp.Rdata
    results: Dir  # everything the task wrote


async def phantompeakqualtools(
    bam: File, prefix: str = "sample", args: str = "", args2: str = ""
) -> PhantompeakqualtoolsResult:
    """Strand cross-correlation metrics of ``bam``; ``args`` go to Rscript, ``args2`` to run_spp.R."""
    results = await phantompeakqualtools_cmd(bam=bam, prefix=prefix, args=args, args2=args2)

    async def pick(name: str) -> File:
        file = await results.get_file(name)
        if file is None:
            raise FileNotFoundError(f"phantompeakqualtools wrote no {name}")
        return file

    return PhantompeakqualtoolsResult(
        spp=await pick(f"{prefix}.spp.out"),
        pdf=await pick(f"{prefix}.spp.pdf"),
        rdata=await pick(f"{prefix}.spp.Rdata"),
        results=results,
    )
