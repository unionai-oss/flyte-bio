"""scimap — spatial single-cell analysis of multiplexed imaging quantification.

Exposes :func:`scimap_mcmicro`, which runs scimap's standard workflow on a
cell-by-feature table (e.g. from :func:`flyte_bio.modules.mcquant.mcquant`):
it clusters the cells by marker expression (k-means, Leiden, PhenoGraph) and
by spatial neighbourhood (spatial k-means). It writes the
AnnData ``<table stem>.h5ad`` and the same as a table, ``<table stem>.h5ad.csv``,
whose columns are the input's plus ``imageid`` and one per clustering.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

SCIMAP_IMAGE = "docker.io/labsyspharm/scimap:2.1.3"
DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory="6Gi")

scimap_mcmicro_cmd = shell.create(
    name="scimap_mcmicro",
    image=SCIMAP_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"table": File},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        TABLE=({{inputs.table}})
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        mkdir numba_cache_dir
        export HOME="$W" MPLCONFIGDIR="$W/.matplotlib" NUMBA_CACHE_DIR="$W/numba_cache_dir"
        scimap-mcmicro "${{TABLE[0]}}" -o . || exit 1
        shopt -s nullglob
        mv ./*.csv ./*.h5ad "$OUT"/
    """,
)

env = flyte.TaskEnvironment.from_task("scimap", scimap_mcmicro_cmd.as_task())


@dataclass
class ScimapResult:
    csv: File | None  # <table stem>.h5ad.csv
    h5ad: File | None  # <table stem>.h5ad
    results: Dir  # everything scimap wrote


async def scimap_mcmicro(table: File) -> ScimapResult:
    """Cluster the cells of a cell-by-feature ``table``; as upstream, either output may be missing."""
    results = await scimap_mcmicro_cmd(table=table)
    stem = table.name.rsplit(".", 1)[0]
    return ScimapResult(
        csv=await results.get_file(f"{stem}.h5ad.csv"),
        h5ad=await results.get_file(f"{stem}.h5ad"),
        results=results,
    )
