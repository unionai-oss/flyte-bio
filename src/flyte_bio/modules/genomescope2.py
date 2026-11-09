"""genomescope2 — genome size, heterozygosity and repeat content from a k-mer histogram.

Exposes :func:`genomescope2`, which fits GenomeScope's model to a k-mer
count histogram (e.g. from :func:`flyte_bio.modules.meryl.meryl_histogram`).
Pass the k-mer length and ploidy in ``args`` (``-k 21 -p 2``). As upstream,
the outputs are named ``<prefix>_<name>``:

- ``<prefix>_{linear,transformed_linear,log,transformed_log}_plot.png``;
- ``<prefix>_model.txt`` and ``<prefix>_summary.txt``;
- ``<prefix>_lookup_table.txt``, ``<prefix>_fitted_hist.png`` (with
  ``--fitted_hist``) and a JSON report, when GenomeScope writes them.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

GENOMESCOPE2_IMAGE = "community.wave.seqera.io/library/genomescope2:2.1.0--a4c756d0a4552c53"
DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory="4Gi")

genomescope2_cmd = shell.create(
    name="genomescope2",
    image=GENOMESCOPE2_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"histogram": File, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        HIST=({inputs.histogram})
        P={inputs.prefix}
        ARGS={inputs.args}
        # As upstream, the histogram is read under its own name from the working directory: the summary
        # records the input path.
        W=$(mktemp -d); cd "$W"
        NAME=$(basename "${HIST[0]}")
        ln -s "${HIST[0]}" "$NAME"
        genomescope2 --input "$NAME" $ARGS --output . --name_prefix "$P" || exit 1
        rm "$NAME"
        if [ -f fitted_hist.png ]; then mv fitted_hist.png "${P}_fitted_hist.png"; fi
        if [ -f lookup_table.txt ]; then mv lookup_table.txt "${P}_lookup_table.txt"; fi
        mv ./* {outputs.results}/
    """,
)

env = flyte.TaskEnvironment.from_task("genomescope2", genomescope2_cmd.as_task())


@dataclass
class Genomescope2Result:
    linear_plot_png: File
    transformed_linear_plot_png: File
    log_plot_png: File
    transformed_log_plot_png: File
    model: File
    summary: File
    lookup_table: File | None
    fitted_histogram_png: File | None  # with --fitted_hist
    json_report: File | None
    results: Dir  # everything above


async def genomescope2(histogram: File, prefix: str = "sample", args: str = "") -> Genomescope2Result:
    """Fit GenomeScope's model to ``histogram``."""
    results = await genomescope2_cmd(histogram=histogram, prefix=prefix, args=args)

    async def find(suffix: str) -> File | None:
        return await results.get_file(f"{prefix}_{suffix}")

    async def need(suffix: str) -> File:
        found = await find(suffix)
        if found is None:
            raise FileNotFoundError(f"genomescope2 wrote no {prefix}_{suffix}")
        return found

    reports = [f async for f in results.walk() if f.name.endswith(".json")]
    return Genomescope2Result(
        linear_plot_png=await need("linear_plot.png"),
        transformed_linear_plot_png=await need("transformed_linear_plot.png"),
        log_plot_png=await need("log_plot.png"),
        transformed_log_plot_png=await need("transformed_log_plot.png"),
        model=await need("model.txt"),
        summary=await need("summary.txt"),
        lookup_table=await find("lookup_table.txt"),
        fitted_histogram_png=await find("fitted_hist.png"),
        json_report=reports[0] if reports else None,
        results=results,
    )
