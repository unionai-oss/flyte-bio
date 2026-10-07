"""multiqc — aggregate QC reports into one HTML report.

Exposes :func:`multiqc`, which runs MultiQC over one Dir of QC outputs with
the given config files (applied in order) and an optional sample-name
replacement table, writing ``<prefix>.html``, ``<prefix>_data/`` and
``<prefix>_plots/`` into an output Dir.

MultiQC derives sample names from the file (and directory) names it finds,
so callers lay the Dir out with the names the producing tools would use.
"""

import asyncio
import tempfile
from dataclasses import dataclass
from pathlib import Path

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

MULTIQC_IMAGE = "community.wave.seqera.io/library/multiqc:1.33--ee7739d47738383b"

DEFAULT_RESOURCES = flyte.Resources(cpu=2, memory="8Gi")

# Config files are passed as explicit `--config` options in list order;
# `configs` are named `<index>_<name>` by the wrapper so a glob sorts them in
# that order. `replace_names` is `list[File]` (0 or 1 item) rather than
# `File | None` until flyteorg/flyte#8118 is deployed.
multiqc_cmd = shell.create(
    name="multiqc",
    image=MULTIQC_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"data": Dir, "configs": list[File], "replace_names": list[File], "prefix": str, "args": str},
    defaults={"configs": [], "replace_names": [], "args": ""},
    outputs={"results": Dir},
    script=r"""
        shopt -s nullglob; CONFIGS=({inputs.configs}); REPLACE=({inputs.replace_names}); shopt -u nullglob
        OPTS=()
        for c in "${CONFIGS[@]}"; do OPTS+=(--config "$c"); done
        [ ${#REPLACE[@]} -gt 0 ] && OPTS+=(--replace-names "${REPLACE[0]}")
        ARGS={inputs.args}
        multiqc \
            --force \
            $ARGS \
            "${OPTS[@]}" \
            --filename {inputs.prefix}.html \
            --outdir {outputs.results} \
            {inputs.data}
    """,
)


env = flyte.TaskEnvironment.from_task(
    "multiqc",
    multiqc_cmd.as_task(),
)


async def ordered(configs: list[File]) -> list[File]:
    """Re-stage ``configs`` as ``<index>_<basename>`` so the task's glob keeps list order.

    A list[File] is staged under original basenames and globbed alphabetically;
    config files are tiny, so they're copied through the caller's pod.
    """
    root = Path(tempfile.mkdtemp(prefix="multiqc_configs_"))

    async def restage(i: int, cfg: File) -> File:
        name = f"{i:02d}_{cfg.path.rstrip('/').rsplit('/', 1)[-1]}"
        await cfg.download(str(root / name))
        return await File.from_local(root / name)

    return list(await asyncio.gather(*(restage(i, c) for i, c in enumerate(configs))))


@dataclass
class MultiQCResult:
    report: File  # <prefix>.html
    results: Dir  # report, <prefix>_data/ and <prefix>_plots/


async def multiqc(
    data: Dir,
    configs: list[File] | None = None,
    replace_names: File | None = None,
    prefix: str = "multiqc_report",
    args: str = "",
) -> MultiQCResult:
    """Run MultiQC over ``data``; ``configs`` apply in list order (later ones override)."""
    results = await multiqc_cmd(
        data=data,
        configs=await ordered(configs or []),
        replace_names=[replace_names] if replace_names is not None else [],
        prefix=prefix,
        args=args,
    )
    report = await results.get_file(f"{prefix}.html")
    if report is None:
        raise FileNotFoundError(f"multiqc wrote no {prefix}.html")
    return MultiQCResult(report=report, results=results)
