"""cooler — Hi-C contact matrices in the cooler (``.cool`` / ``.mcool``) format.

Exposes:

- :func:`cooler_cload` — bin a contact list (pairs) into ``<prefix>.cool``, in
  one of cooler's ``cload`` modes (``pairix``, ``tabix``, ``pairs``, ...), at
  ``bin_size`` over a chrom sizes file. The ``pairix`` and ``tabix`` modes need
  the contacts' index (``.px2`` / ``.tbi``), which is linked next to them, and
  run on the task's CPUs (``--nproc``), as upstream.
- :func:`cooler_digest` — the restriction fragments of a genome (FASTA plus
  chrom sizes) for an ``enzyme``, as ``<FASTA name without extension>_<enzyme>.bed``
  (runs of characters other than letters and digits in the enzyme become ``_``).
- :func:`cooler_makebins` — the fixed-size genomic bins of a chrom sizes file,
  as ``<prefix>.bed``.
- :func:`cooler_balance` — matrix-balance a cooler (one resolution of an
  ``.mcool`` with ``resolution``), writing ``<prefix>.<cool extension>`` with
  the ``weight`` column added, on the task's CPUs.
- :func:`cooler_dump` — dump a cooler's table (pixels by default) as
  ``<prefix>.bedpe``.
- :func:`cooler_zoomify` — aggregate a cooler into a multi-resolution
  ``<prefix>.mcool`` (``args`` picks the resolutions, e.g. ``-r 2,4,8``), on
  the task's CPUs.

``args`` is upstream's ``ext.args`` and is split on whitespace.

Unlike upstream, :func:`cooler_balance` balances a copy of its input, not a
symlink to it: inputs are read-only here, and upstream's symlink writes the
weights into the input file itself.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

COOLER_IMAGE = "quay.io/biocontainers/cooler:0.10.4--pyhdfd78af_0"
# Upstream pins dump to the previous release.
COOLER_DUMP_IMAGE = "quay.io/biocontainers/cooler:0.10.3--pyhdfd78af_0"

# Upstream's process_high is 12 CPUs and 72 GB; this is enough for the test data.
# Genome-wide matrices at fine resolutions need more memory.
DEFAULT_RESOURCES = flyte.Resources(cpu=4, memory="16Gi")
MAKEBINS_RESOURCES = flyte.Resources(cpu=1, memory="2Gi")

cooler_cload_cmd = shell.create(
    name="cooler_cload",
    image=COOLER_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "contacts": File,
        "index": File | None,
        "chromsizes": File,
        "mode": str,
        "bin_size": int,
        "prefix": str,
        "args": str,
    },
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        CONTACTS=({{inputs.contacts}}); SIZES=({{inputs.chromsizes}})
        shopt -s nullglob; INDEX=({{inputs.index}}); shopt -u nullglob
        MODE={{inputs.mode}}
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        W=$(mktemp -d); cd "$W"
        NAME=$(basename "${{CONTACTS[0]}}")
        ln -s "${{CONTACTS[0]}}" "$NAME"
        if [ ${{#INDEX[@]}} -gt 0 ]; then ln -s "${{INDEX[0]}}" "$(basename "${{INDEX[0]}}")"; fi
        NPROC=()
        if [ "$MODE" = pairix ] || [ "$MODE" = tabix ]; then NPROC=(--nproc $TASK_CPUS); fi
        cooler cload "$MODE" $ARGS "${{NPROC[@]}}" "${{SIZES[0]}}":{{inputs.bin_size}} "$NAME" "$P.cool" || exit 1
        mv "$P.cool" {{outputs.results}}/
    """,
)

# The FASTA is linked into a working directory: its index (.fai) is written
# next to it, and inputs are read-only.
cooler_digest_cmd = shell.create(
    name="cooler_digest",
    image=COOLER_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"fasta": File, "chromsizes": File, "enzyme": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        FA=({inputs.fasta}); SIZES=({inputs.chromsizes})
        ENZYME={inputs.enzyme}
        ARGS={inputs.args}
        W=$(mktemp -d); cd "$W"
        NAME=$(basename "${FA[0]}")
        ln -s "${FA[0]}" "$NAME"
        LABEL=$(printf '%s' "$ENZYME" | sed -E 's/[^0-9a-zA-Z]+/_/g')
        cooler digest $ARGS -o {outputs.results}/"${NAME%.*}_$LABEL.bed" "${SIZES[0]}" "$NAME" "$ENZYME" || exit 1
    """,
)

cooler_makebins_cmd = shell.create(
    name="cooler_makebins",
    image=COOLER_IMAGE,
    resources=MAKEBINS_RESOURCES,
    inputs={"chromsizes": File, "bin_size": int, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        SIZES=({inputs.chromsizes})
        ARGS={inputs.args}
        cooler makebins $ARGS "${SIZES[0]}" {inputs.bin_size} > {outputs.results}/{inputs.prefix}.bed || exit 1
    """,
)

# `resolution` is empty for a single-resolution cooler.
cooler_balance_cmd = shell.create(
    name="cooler_balance",
    image=COOLER_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"cool": File, "resolution": str, "prefix": str, "args": str},
    defaults={"resolution": "", "args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        COOL=({{inputs.cool}})
        RESOLUTION={{inputs.resolution}}
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        NAME=$(basename "${{COOL[0]}}")
        OUT={{outputs.results}}/"$P.${{NAME##*.}}"
        cp "${{COOL[0]}}" "$OUT" || exit 1
        chmod u+w "$OUT"
        URI="$OUT"
        if [ -n "$RESOLUTION" ]; then URI="$OUT::/resolutions/$RESOLUTION"; fi
        cooler balance $ARGS -p $TASK_CPUS "$URI" || exit 1
    """,
)

cooler_dump_cmd = shell.create(
    name="cooler_dump",
    image=COOLER_DUMP_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"cool": File, "resolution": str, "prefix": str, "args": str},
    defaults={"resolution": "", "args": ""},
    outputs={"results": Dir},
    script=r"""
        COOL=({inputs.cool})
        RESOLUTION={inputs.resolution}
        ARGS={inputs.args}
        URI="${COOL[0]}"
        if [ -n "$RESOLUTION" ]; then URI="$URI::/resolutions/$RESOLUTION"; fi
        cooler dump $ARGS -o {outputs.results}/{inputs.prefix}.bedpe "$URI" || exit 1
    """,
)

cooler_zoomify_cmd = shell.create(
    name="cooler_zoomify",
    image=COOLER_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"cool": File, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        COOL=({{inputs.cool}})
        ARGS={{inputs.args}}
        cooler zoomify $ARGS -n $TASK_CPUS -o {{outputs.results}}/{{inputs.prefix}}.mcool "${{COOL[0]}}" || exit 1
    """,
)

cooler_tools_env = flyte.TaskEnvironment.from_task(
    "cooler_tools",
    cooler_cload_cmd.as_task(),
    cooler_digest_cmd.as_task(),
    cooler_makebins_cmd.as_task(),
    cooler_balance_cmd.as_task(),
    cooler_zoomify_cmd.as_task(),
)
cooler_dump_env = flyte.TaskEnvironment.from_task("cooler_dump", cooler_dump_cmd.as_task())
env = flyte.TaskEnvironment(name="cooler", depends_on=[cooler_tools_env, cooler_dump_env])


@dataclass
class CoolerResult:
    output: File  # the file the tool wrote, named as upstream
    results: Dir


def file_name(file: File) -> str:
    return file.path.rstrip("/").rsplit("/", 1)[-1]


async def pick(results: Dir, name: str) -> CoolerResult:
    found = await results.get_file(name)
    if found is None:
        raise FileNotFoundError(f"cooler wrote no {name}")
    return CoolerResult(output=found, results=results)


async def cooler_cload(
    contacts: File,
    chromsizes: File,
    mode: str,
    bin_size: int,
    index: File | None = None,
    prefix: str = "sample",
    args: str = "",
) -> CoolerResult:
    """Bin ``contacts`` into ``<prefix>.cool`` with ``cooler cload <mode>`` at ``bin_size``."""
    results = await cooler_cload_cmd(
        contacts=contacts, index=index, chromsizes=chromsizes, mode=mode, bin_size=bin_size, prefix=prefix, args=args
    )
    return await pick(results, f"{prefix}.cool")


async def cooler_digest(fasta: File, chromsizes: File, enzyme: str, args: str = "") -> CoolerResult:
    """``<FASTA name without extension>_<enzyme>.bed``: the ``enzyme`` fragments of the genome."""
    results = await cooler_digest_cmd(fasta=fasta, chromsizes=chromsizes, enzyme=enzyme, args=args)
    beds = [file for file in await results.list_files() if file_name(file).endswith(".bed")]
    if len(beds) != 1:
        raise FileNotFoundError(f"cooler digest wrote {len(beds)} BED files, expected 1")
    return CoolerResult(output=beds[0], results=results)


async def cooler_makebins(chromsizes: File, bin_size: int, prefix: str = "sample", args: str = "") -> CoolerResult:
    """``<prefix>.bed``: the ``bin_size`` bins of the genome in ``chromsizes``."""
    results = await cooler_makebins_cmd(chromsizes=chromsizes, bin_size=bin_size, prefix=prefix, args=args)
    return await pick(results, f"{prefix}.bed")


async def cooler_balance(
    cool: File, resolution: int | None = None, prefix: str = "sample", args: str = ""
) -> CoolerResult:
    """A balanced copy of ``cool`` (its ``resolution`` for an ``.mcool``), named ``<prefix>.<extension>``."""
    name = file_name(cool)
    extension = name.rsplit(".", 1)[-1]
    if name == f"{prefix}.{extension}":
        raise ValueError('Input and output names are the same, use "task.ext.prefix" to disambiguate!')
    results = await cooler_balance_cmd(
        cool=cool, resolution="" if resolution is None else str(resolution), prefix=prefix, args=args
    )
    return await pick(results, f"{prefix}.{extension}")


async def cooler_dump(
    cool: File, resolution: int | None = None, prefix: str = "sample", args: str = ""
) -> CoolerResult:
    """``<prefix>.bedpe``: a table of ``cool`` (its ``resolution`` for an ``.mcool``)."""
    results = await cooler_dump_cmd(
        cool=cool, resolution="" if resolution is None else str(resolution), prefix=prefix, args=args
    )
    return await pick(results, f"{prefix}.bedpe")


async def cooler_zoomify(cool: File, prefix: str = "sample", args: str = "") -> CoolerResult:
    """``<prefix>.mcool``: ``cool`` aggregated to coarser resolutions."""
    results = await cooler_zoomify_cmd(cool=cool, prefix=prefix, args=args)
    return await pick(results, f"{prefix}.mcool")
