"""glnexus — GLnexus joint genotyping of gVCFs into a multi-sample BCF.

Exposes :func:`glnexus`, which merges and jointly genotypes several samples'
gVCFs (``glnexus_cli``), optionally restricted to the regions of a BED file
(``--bed``). The output is ``<prefix>.bcf``.

As upstream:

- ``--threads`` and ``--mem-gbytes`` follow the task's CPUs and memory
  (whole GiB).
- The configuration preset is chosen in ``args`` (``--config gatk``,
  ``DeepVariant``, …); GLnexus's default is ``gatk``. A ``custom_config``
  YAML is staged in the working directory under its own name, so ``args`` can
  name it (``--config <name>.yml``); like upstream, it isn't passed otherwise.
- GLnexus keeps its working database (``GLnexus.DB``) in a scratch working
  directory.
"""

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

GLNEXUS_IMAGE = "community.wave.seqera.io/library/bcftools_glnexus:cf380f1a6410f606"

# Upstream labels it process_high; these fit test-sized and modest cohorts.
DEFAULT_RESOURCES = flyte.Resources(cpu=4, memory="16Gi")

glnexus_cmd = shell.create(
    name="glnexus",
    image=GLNEXUS_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"gvcfs": list[File], "custom_config": File | None, "bed": File | None, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        GVCFS=({{inputs.gvcfs}})
        shopt -s nullglob; CONFIG=({{inputs.custom_config}}); BED=({{inputs.bed}}); shopt -u nullglob
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        W=$(mktemp -d); cd "$W"
        NAMES=()
        for F in "${{GVCFS[@]}}"; do ln -s "$F" "$(basename "$F")"; NAMES+=("$(basename "$F")"); done
        for F in "${{CONFIG[@]}}" "${{BED[@]}}"; do ln -s "$F" "$(basename "$F")"; done
        REGIONS=()
        if [ ${{#BED[@]}} -gt 0 ]; then REGIONS=(--bed "$(basename "${{BED[0]}}")"); fi
        glnexus_cli \
            --threads $TASK_CPUS \
            --mem-gbytes $(( TASK_MEMORY_MB / 1024 )) \
            "${{REGIONS[@]}}" \
            $ARGS \
            "${{NAMES[@]}}" \
            > "$P.bcf" || exit 1
        mv "$P.bcf" {{outputs.results}}/
    """,
)

env = flyte.TaskEnvironment.from_task("glnexus", glnexus_cmd.as_task())


async def glnexus(
    gvcfs: list[File],
    prefix: str,
    custom_config: File | None = None,
    bed: File | None = None,
    args: str = "",
) -> File:
    """Jointly genotype ``gvcfs`` (in the given order) into ``<prefix>.bcf``."""
    results = await glnexus_cmd(gvcfs=gvcfs, custom_config=custom_config, bed=bed, prefix=prefix, args=args)
    bcf = await results.get_file(f"{prefix}.bcf")
    if bcf is None:
        raise FileNotFoundError(f"glnexus wrote no {prefix}.bcf")
    return bcf
