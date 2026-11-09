"""fmhfunprofiler — fmh-funprofiler, functional (KEGG ortholog) profiling of metagenomes with FracMinHash sketches.

Exposes :func:`fmhfunprofiler`, which sketches one sample's reads at k-mer size
``ksize`` and scale ``scaled``, and profiles them against a sourmash sketch of
KEGG ortholog sequences (``ko_sketch``), writing the KO abundance table
``<prefix>.csv``. ``args`` are funcprofiler's options, passed before its
positional arguments, as upstream.
"""

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

FMHFUNPROFILER_IMAGE = "quay.io/biocontainers/fmh-funprofiler:1.1.1--pyh106432d_0"
# Upstream labels it process_medium; funcprofiler is single-threaded.
DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory="8Gi")

fmhfunprofiler_cmd = shell.create(
    name="fmhfunprofiler",
    image=FMHFUNPROFILER_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"reads": File, "ko_sketch": File, "ksize": int, "scaled": int, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        READS=({inputs.reads})
        SKETCH=({inputs.ko_sketch})
        P={inputs.prefix}
        ARGS={inputs.args}
        OUT={outputs.results}
        W=$(mktemp -d); cd "$W"
        ln -s "${READS[0]}" .
        ln -s "${SKETCH[0]}" .
        funcprofiler \
            $ARGS \
            "$(basename "${READS[0]}")" \
            "$(basename "${SKETCH[0]}")" \
            {inputs.ksize} \
            {inputs.scaled} \
            "$P.csv" || exit 1
        mv "$P.csv" "$OUT"/
    """,
)

env = flyte.TaskEnvironment.from_task("fmhfunprofiler", fmhfunprofiler_cmd.as_task())


async def fmhfunprofiler(
    reads: File, ko_sketch: File, ksize: int, scaled: int, prefix: str = "sample", args: str = ""
) -> File:
    """The KEGG ortholog profile ``<prefix>.csv`` of ``reads`` against ``ko_sketch``."""
    results = await fmhfunprofiler_cmd(
        reads=reads, ko_sketch=ko_sketch, ksize=ksize, scaled=scaled, prefix=prefix, args=args
    )
    table = await results.get_file(f"{prefix}.csv")
    if table is None:
        raise FileNotFoundError(f"funcprofiler wrote no {prefix}.csv")
    return table
