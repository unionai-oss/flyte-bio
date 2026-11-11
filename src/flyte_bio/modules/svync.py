"""svync — svync, harmonise structural variant VCFs from different callers.

Exposes :func:`svync`, which rewrites a VCF's INFO and FORMAT fields as a
YAML ``config`` describes (adding, renaming or computing them) and writes
``<prefix>.vcf.gz`` with its ``.tbi``. As upstream, ``args`` goes to svync,
``bgzip_args`` to ``bgzip`` and ``tabix_args`` to ``tabix``, and ``bgzip``
gets the task's CPUs as threads.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

# Upstream's main.nf names the bare biocontainer (Docker Hub); this is its quay.io mirror.
SVYNC_IMAGE = "quay.io/biocontainers/svync:0.3.0--h9ee0642_0"
# Upstream's process_low.
DEFAULT_RESOURCES = flyte.Resources(cpu=2, memory="4Gi")

svync_cmd = shell.create(
    name="svync",
    image=SVYNC_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "vcf": File,
        "tbi": File | None,
        "config": File,
        "prefix": str,
        "args": str,
        "bgzip_args": str,
        "tabix_args": str,
    },
    defaults={"args": "", "bgzip_args": "", "tabix_args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        set -o pipefail
        VCF=({{inputs.vcf}}); CONFIG=({{inputs.config}})
        shopt -s nullglob; TBI=({{inputs.tbi}}); shopt -u nullglob
        P={{inputs.prefix}}
        ARGS={{inputs.args}}; BGZIP_ARGS={{inputs.bgzip_args}}; TABIX_ARGS={{inputs.tabix_args}}
        W=$(mktemp -d); cd "$W"
        for f in "${{VCF[0]}}" "${{TBI[@]}}"; do ln -s "$f" "$(basename "$f")"; done
        # svync is a Go tool: exit codes like 255 would be retried forever, so any failure exits 1.
        svync $ARGS --config "${{CONFIG[0]}}" --input "$(basename "${{VCF[0]}}")" \
            | bgzip --threads $TASK_CPUS $BGZIP_ARGS > "$P.vcf.gz" || exit 1
        tabix $TABIX_ARGS "$P.vcf.gz" || exit 1
        mv "$P.vcf.gz" "$P.vcf.gz.tbi" {{outputs.results}}/
    """,
)

env = flyte.TaskEnvironment.from_task("svync", svync_cmd.as_task())


@dataclass
class SvyncResult:
    vcf: File  # <prefix>.vcf.gz
    tbi: File  # <prefix>.vcf.gz.tbi


async def svync(
    vcf: File,
    config: File,
    prefix: str,
    tbi: File | None = None,
    args: str = "",
    bgzip_args: str = "",
    tabix_args: str = "",
) -> SvyncResult:
    """``vcf`` with its fields rewritten as the YAML ``config`` says: ``<prefix>.vcf.gz`` and its ``.tbi``."""
    if vcf.name == f"{prefix}.vcf.gz":
        raise ValueError("input and output names are the same; choose another prefix")
    results = await svync_cmd(
        vcf=vcf, tbi=tbi, config=config, prefix=prefix, args=args, bgzip_args=bgzip_args, tabix_args=tabix_args
    )
    out = await results.get_file(f"{prefix}.vcf.gz")
    index = await results.get_file(f"{prefix}.vcf.gz.tbi")
    if out is None or index is None:
        raise FileNotFoundError(f"svync wrote no {prefix}.vcf.gz or its .tbi")
    return SvyncResult(vcf=out, tbi=index)
