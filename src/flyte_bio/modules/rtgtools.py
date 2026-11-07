"""rtgtools — RTG Tools.

Exposes :func:`rtgtools_format`, which converts sequences into RTG's SDF
format (``rtg format``): a FASTA reference (as ``rtg vcfeval`` takes), FASTQ
reads or SAM/BAM reads, as ``args``' ``--format`` says (FASTA by default).
With ``single_end`` (the default, as upstream) ``input1`` is formatted on its
own; otherwise ``input1`` and ``input2`` are the left and right reads.
``sam_rg`` is a file holding the read group header line to record
(``--sam-rg``).

The SDF is upstream's ``<prefix>.sdf`` directory. As upstream, RTG's memory
(``RTG_MEM``) is 80% of the task's memory.
"""

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import JAVA_HEAP_MB, TASK_LIMITS

RTGTOOLS_IMAGE = "community.wave.seqera.io/library/rtg-tools:3.13--3465421f1b0be0ce"
# Upstream's process_low.
DEFAULT_RESOURCES = flyte.Resources(cpu=2, memory="4Gi")

rtgtools_format_cmd = shell.create(
    name="rtgtools_format",
    image=RTGTOOLS_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "input1": File,
        "input2": File | None,
        "sam_rg": File | None,
        "single_end": bool,
        "prefix": str,
        "args": str,
    },
    defaults={"single_end": True, "args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        IN1=({{inputs.input1}})
        shopt -s nullglob; IN2=({{inputs.input2}}); RG=({{inputs.sam_rg}}); shopt -u nullglob
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"

        if [ {{inputs.single_end}} = true ]; then
            INPUT=("${{IN1[0]}}")
        else
            if [ ${{#IN2[@]}} -eq 0 ]; then
                echo "paired-end formatting needs input2" >&2
                exit 1
            fi
            INPUT=(--left "${{IN1[0]}}" --right "${{IN2[0]}}")
        fi
        RG_OPTION=()
        if [ ${{#RG[@]}} -gt 0 ]; then RG_OPTION=(--sam-rg "${{RG[0]}}"); fi

        rtg RTG_MEM={JAVA_HEAP_MB}M format \
            $ARGS \
            "${{RG_OPTION[@]}}" \
            --output "$P.sdf" \
            "${{INPUT[@]}}"
        mv "$P.sdf" "$OUT"/
    """,
)

env = flyte.TaskEnvironment.from_task("rtgtools", rtgtools_format_cmd.as_task())


async def rtgtools_format(
    input1: File,
    input2: File | None = None,
    sam_rg: File | None = None,
    single_end: bool = True,
    prefix: str = "sample",
    args: str = "",
) -> Dir:
    """``input1`` (and, paired, ``input2``) in RTG's SDF format: the ``<prefix>.sdf`` Dir."""
    results = await rtgtools_format_cmd(
        input1=input1, input2=input2, sam_rg=sam_rg, single_end=single_end, prefix=prefix, args=args
    )
    return Dir.from_existing_remote(results.path.rstrip("/") + f"/{prefix}.sdf")
