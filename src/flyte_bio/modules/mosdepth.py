"""mosdepth — read-depth summaries of BAM / CRAM files.

Exposes :func:`mosdepth`, which writes mosdepth's ``<prefix>.mosdepth.*``
reports (global and per-region coverage distributions, the summary) and,
depending on the options, per-base / per-region depth BEDs, into one output
Dir. The alignment and its index are linked side by side first, as mosdepth
expects; a CRAM also needs its ``fasta``.
"""

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

MOSDEPTH_IMAGE = "community.wave.seqera.io/library/htslib_mosdepth_gzip:4108dd38be84e40a"

MOSDEPTH_CPUS = 4
DEFAULT_RESOURCES = flyte.Resources(cpu=MOSDEPTH_CPUS, memory="8Gi")

mosdepth_cmd = shell.create(
    name="mosdepth",
    image=MOSDEPTH_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"bam": File, "index": File, "bed": File | None, "fasta": File | None, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        BAM=({{inputs.bam}}); IDX=({{inputs.index}})
        shopt -s nullglob; BED=({{inputs.bed}}); FA=({{inputs.fasta}}); shopt -u nullglob
        ARGS={{inputs.args}}
        W=$(mktemp -d); NAME=$(basename "${{BAM[0]}}")
        EXT=bai; [[ "$NAME" == *.cram ]] && EXT=crai
        ln -s "${{BAM[0]}}" "$W/$NAME"
        ln -s "${{IDX[0]}}" "$W/$NAME.$EXT"
        BY=(); [ ${{#BED[@]}} -gt 0 ] && BY=(--by "${{BED[0]}}")
        # htslib indexes the FASTA beside it if needed, so give it a writable home.
        REF=()
        if [ ${{#FA[@]}} -gt 0 ]; then ln -s "${{FA[0]}}" "$W/genome.fa"; REF=(--fasta "$W/genome.fa"); fi
        mosdepth --threads {MOSDEPTH_CPUS} "${{BY[@]}}" "${{REF[@]}}" $ARGS \
            {{outputs.results}}/{{inputs.prefix}} "$W/$NAME"
    """,
)


env = flyte.TaskEnvironment.from_task(
    "mosdepth",
    mosdepth_cmd.as_task(),
)


async def mosdepth(
    bam: File,
    index: File,
    prefix: str,
    bed: File | None = None,
    fasta: File | None = None,
    args: str = "",
) -> Dir:
    """Depth reports for ``bam`` (BAM or CRAM, with its ``.bai`` / ``.crai``), as ``<prefix>.*``.

    ``bed`` reports depth per region (``--by``); ``args`` may give ``--by
    <window>`` instead, but not both, as upstream enforces.
    """
    if bed is not None and ("--by" in args or "-b " in args):
        raise ValueError("--by can only be given once: either a bed or --by in args")
    return await mosdepth_cmd(bam=bam, index=index, bed=bed, fasta=fasta, prefix=prefix, args=args)
