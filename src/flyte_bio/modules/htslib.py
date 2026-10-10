"""htslib — bgzip and tabix.

Exposes :func:`htslib_bgziptabix`, which bgzip-compresses a file (plain,
gzip, bzip2 or xz input, as upstream detects it) to
``<prefix>[.<ext>].gz`` and, optionally, tabix-indexes it.
"""

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

HTSLIB_IMAGE = "community.wave.seqera.io/library/htslib_xz:32f2772a564b3cd2"

HTSLIB_CPUS = 2
DEFAULT_RESOURCES = flyte.Resources(cpu=HTSLIB_CPUS, memory="2Gi")

# Upstream's compress path: an already-BGZF input is reused as is, other
# compressions are recompressed, and a plain file is compressed.
htslib_bgziptabix_cmd = shell.create(
    name="htslib_bgziptabix",
    image=HTSLIB_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"infile": File, "prefix": str, "out_ext": str, "make_index": bool, "args": str, "args2": str},
    defaults={"out_ext": "", "make_index": True, "args": "", "args2": ""},
    outputs={"results": Dir},
    script=rf"""
        IN=({{inputs.infile}})
        EXT={{inputs.out_ext}}
        ARGS={{inputs.args}}
        ARGS2={{inputs.args2}}
        OUT={{outputs.results}}/{{inputs.prefix}}${{EXT:+.$EXT}}.gz
        BGZIP="bgzip -c $ARGS -@ {HTSLIB_CPUS}"
        case "$(htsfile "${{IN[0]}}")" in
            *BGZF-compressed*) cp "${{IN[0]}}" "$OUT" ;;
            *gzip-compressed*) bgzip -d -c -@ {HTSLIB_CPUS} "${{IN[0]}}" | $BGZIP > "$OUT" ;;
            *bzip2-compressed*) bzcat "${{IN[0]}}" | $BGZIP > "$OUT" ;;
            *XZ-compressed*) xzcat "${{IN[0]}}" | $BGZIP > "$OUT" ;;
            *) $BGZIP "${{IN[0]}}" > "$OUT" ;;
        esac
        if [ {{inputs.make_index}} = true ]; then tabix -@ {HTSLIB_CPUS} $ARGS2 -f "$OUT"; fi
    """,
)


env = flyte.TaskEnvironment.from_task(
    "htslib",
    htslib_bgziptabix_cmd.as_task(),
)


async def htslib_bgziptabix(
    infile: File,
    prefix: str,
    out_ext: str = "",
    make_index: bool = True,
    args: str = "",
    args2: str = "",
) -> tuple[File, File | None]:
    """bgzip ``infile`` to ``<prefix>[.<out_ext>].gz``; return it and its ``.tbi`` (None without ``make_index``).

    ``args`` go to bgzip, ``args2`` to tabix (e.g. ``-p vcf``).
    """
    results = await htslib_bgziptabix_cmd(
        infile=infile, prefix=prefix, out_ext=out_ext, make_index=make_index, args=args, args2=args2
    )
    name = f"{prefix}.{out_ext}.gz" if out_ext else f"{prefix}.gz"
    compressed = await results.get_file(name)
    if compressed is None:
        raise FileNotFoundError(f"bgzip wrote no {name}")
    index = await results.get_file(f"{name}.tbi") if make_index else None
    if make_index and index is None:
        raise FileNotFoundError(f"tabix wrote no {name}.tbi")
    return compressed, index
