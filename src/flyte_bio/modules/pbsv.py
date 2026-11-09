"""pbsv — structural-variant calling from PacBio reads with pbsv.

Exposes:

- :func:`pbsv_discover` — find structural-variant signatures in one
  sample's aligned BAM: ``<prefix>.svsig.gz``.
- :func:`pbsv_call` — call structural variants from one or more samples'
  signatures against the reference: ``<prefix>.vcf``. ``-j`` follows the
  task's CPUs.

``args`` go to pbsv and aren't shell-quoted. Upstream's discover module also
takes the reference but doesn't pass it to pbsv, so :func:`pbsv_discover`
doesn't take one. The inputs are linked into the working directory under
their own names.
"""

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

PBSV_IMAGE = "community.wave.seqera.io/library/pbsv:2.11.0--c85e7f17330a07c9"
# Upstream labels discover process_single and call process_medium; its tests run with 4 CPUs.
DISCOVER_RESOURCES = flyte.Resources(cpu=1, memory="4Gi")
CALL_RESOURCES = flyte.Resources(cpu=4, memory="8Gi")

pbsv_discover_cmd = shell.create(
    name="pbsv_discover",
    image=PBSV_IMAGE,
    resources=DISCOVER_RESOURCES,
    inputs={"bam": File, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        BAM=({inputs.bam})
        P={inputs.prefix}
        ARGS={inputs.args}
        W=$(mktemp -d); cd "$W"
        BAM_NAME=$(basename "${BAM[0]}"); ln -s "${BAM[0]}" "$BAM_NAME"
        pbsv discover $ARGS "$BAM_NAME" "$P.svsig.gz" || exit 1
        mv "$P.svsig.gz" {outputs.results}/
    """,
)

pbsv_call_cmd = shell.create(
    name="pbsv_call",
    image=PBSV_IMAGE,
    resources=CALL_RESOURCES,
    inputs={"svsig": list[File], "fasta": File, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        SVSIG=({{inputs.svsig}}); FA=({{inputs.fasta}})
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        W=$(mktemp -d); cd "$W"
        FA_NAME=$(basename "${{FA[0]}}"); ln -s "${{FA[0]}}" "$FA_NAME"
        NAMES=()
        for F in "${{SVSIG[@]}}"; do NAMES+=("$(basename "$F")"); ln -s "$F" "${{NAMES[-1]}}"; done
        pbsv call $ARGS -j $TASK_CPUS "$FA_NAME" "${{NAMES[@]}}" "$P.vcf" || exit 1
        mv "$P.vcf" {{outputs.results}}/
    """,
)

env = flyte.TaskEnvironment.from_task("pbsv", pbsv_discover_cmd.as_task(), pbsv_call_cmd.as_task())


async def pbsv_discover(bam: File, prefix: str = "sample", args: str = "") -> File:
    """Structural-variant signatures of the reads in ``bam``: ``<prefix>.svsig.gz``."""
    results = await pbsv_discover_cmd(bam=bam, prefix=prefix, args=args)
    svsig = await results.get_file(f"{prefix}.svsig.gz")
    if svsig is None:
        raise FileNotFoundError(f"pbsv discover wrote no {prefix}.svsig.gz")
    return svsig


async def pbsv_call(svsig: list[File], fasta: File, prefix: str = "sample", args: str = "") -> File:
    """Structural variants called from the signatures in ``svsig`` against ``fasta``: ``<prefix>.vcf``."""
    results = await pbsv_call_cmd(svsig=svsig, fasta=fasta, prefix=prefix, args=args)
    vcf = await results.get_file(f"{prefix}.vcf")
    if vcf is None:
        raise FileNotFoundError(f"pbsv call wrote no {prefix}.vcf")
    return vcf
