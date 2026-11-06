"""htslib — bgzip and tabix.

Exposes:

- :func:`htslib_bgziptabix` — bgzip-compress a file (plain, gzip, bzip2 or
  xz input, as upstream detects it) to ``<prefix>[.<ext>].gz`` and,
  optionally, tabix-index it.
- :func:`tabix_tabix` — tabix-index a bgzipped file (``<name>.tbi``, or
  ``<name>.csi`` with ``--csi``), or, given a regions file, extract those
  regions (with ``-R``) into ``<prefix><suffix>``, bgzipped. The suffix is the
  input name from its first dot (``.vcf.gz``).
- :func:`tabix_bgzip` — bgzip a file to ``<prefix>.<ext>.gz``, or, if its
  extension is ``gz``, ``bgz`` or ``bgzf``, decompress it to
  ``<prefix>.<ext>`` (``<ext>`` from the name without that extension). With
  ``-i`` in ``args`` it also writes the index ``<output>.gzi``.

:func:`tabix_tabix` and :func:`tabix_bgzip` follow upstream's tabix/tabix
and tabix/bgzip modules at their last working version (modules commit
2fb71fb): later upstream versions only fail with a deprecation notice that
points to htslib/bgziptabix. They run in those modules' htslib 1.21 image,
so ``env`` aggregates two environments.
"""

import os
import re

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

HTSLIB_IMAGE = "community.wave.seqera.io/library/htslib_xz:32f2772a564b3cd2"

# Upstream's tabix/tabix and tabix/bgzip image.
TABIX_IMAGE = "community.wave.seqera.io/library/htslib:1.21--ff8e28a189fbecaa"

HTSLIB_CPUS = 2
DEFAULT_RESOURCES = flyte.Resources(cpu=HTSLIB_CPUS, memory="2Gi")
# Upstream's process_single: one CPU.
TABIX_RESOURCES = flyte.Resources(cpu=1, memory="2Gi")

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
        {TASK_LIMITS}
        IN=({{inputs.infile}})
        EXT={{inputs.out_ext}}
        ARGS={{inputs.args}}
        ARGS2={{inputs.args2}}
        OUT={{outputs.results}}/{{inputs.prefix}}${{EXT:+.$EXT}}.gz
        BGZIP="bgzip -c $ARGS -@ $TASK_CPUS"
        case "$(htsfile "${{IN[0]}}")" in
            *BGZF-compressed*) cp "${{IN[0]}}" "$OUT" ;;
            *gzip-compressed*) bgzip -d -c -@ $TASK_CPUS "${{IN[0]}}" | $BGZIP > "$OUT" ;;
            *bzip2-compressed*) bzcat "${{IN[0]}}" | $BGZIP > "$OUT" ;;
            *XZ-compressed*) xzcat "${{IN[0]}}" | $BGZIP > "$OUT" ;;
            *) $BGZIP "${{IN[0]}}" > "$OUT" ;;
        esac
        if [ {{inputs.make_index}} = true ]; then tabix -@ $TASK_CPUS $ARGS2 -f "$OUT"; fi
    """,
)


# tabix writes its index next to the input, so it runs on a symlink in a
# working directory, beside the given index (which -R needs). Only an index
# made here is kept: with regions, nothing is indexed.
tabix_tabix_cmd = shell.create(
    name="tabix_tabix",
    image=TABIX_IMAGE,
    resources=TABIX_RESOURCES,
    inputs={"tab": File, "tai": File | None, "regions": File | None, "prefix": str, "args": str},
    defaults={"prefix": "", "args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        shopt -s nullglob
        TAB=({{inputs.tab}})
        TAI=({{inputs.tai}})
        REGIONS=({{inputs.regions}})
        ARGS={{inputs.args}}
        NAME=$(basename "${{TAB[0]}}")
        SUFFIX=""
        if [[ "$NAME" == *.* ]]; then SUFFIX=".${{NAME#*.}}"; fi
        W=$(mktemp -d)
        ln -s "${{TAB[0]}}" "$W/$NAME"
        for f in "${{TAI[@]}}"; do ln -s "$f" "$W/$(basename "$f")"; done
        cd "$W"
        if [ ${{#REGIONS[@]}} -gt 0 ]; then
            tabix -R "${{REGIONS[0]}}" --threads $TASK_CPUS $ARGS "$NAME" \
                | bgzip --threads $TASK_CPUS > {{outputs.results}}/{{inputs.prefix}}"$SUFFIX"
        else
            tabix --threads $TASK_CPUS $ARGS "$NAME"
            for f in "$NAME".tbi "$NAME".csi; do
                if [ -f "$f" ] && [ ! -L "$f" ]; then cp "$f" {{outputs.results}}/; fi
            done
        fi
    """,
)


# Runs in the output directory, so a relative `-I <name>.gzi` lands there.
tabix_bgzip_cmd = shell.create(
    name="tabix_bgzip",
    image=TABIX_IMAGE,
    resources=TABIX_RESOURCES,
    inputs={"infile": File, "output": str, "decompress": bool, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        IN=({{inputs.infile}})
        ARGS={{inputs.args}}
        D=""
        if [ {{inputs.decompress}} = true ]; then D=-d; fi
        cd {{outputs.results}}
        bgzip $D -c $ARGS -@$TASK_CPUS "${{IN[0]}}" > {{inputs.output}}
    """,
)


# One image per environment: tabix/tabix and tabix/bgzip use upstream's htslib 1.21 image.
htslib_bgziptabix_env = flyte.TaskEnvironment.from_task("htslib_bgziptabix", htslib_bgziptabix_cmd.as_task())
tabix_env = flyte.TaskEnvironment.from_task("tabix", tabix_tabix_cmd.as_task(), tabix_bgzip_cmd.as_task())
env = flyte.TaskEnvironment(name="htslib", depends_on=[htslib_bgziptabix_env, tabix_env])


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


def file_name(file: File) -> str:
    return file.name or os.path.basename(file.path.rstrip("/"))


async def tabix_tabix(
    tab: File,
    tai: File | None = None,
    regions: File | None = None,
    prefix: str = "",
    args: str = "",
) -> tuple[File | None, File | None]:
    """tabix ``tab``; return ``(index, extracted)``, one of them None.

    Without ``regions``: ``index`` is ``<tab name>.tbi`` (``.csi`` with
    ``--csi`` in ``args``, e.g. ``-p vcf --csi``). With ``regions`` (and
    ``tab``'s index as ``tai``): ``extracted`` holds those regions, bgzipped,
    as ``<prefix><suffix>``; ``prefix`` defaults to ``tab``'s name up to its
    first dot.
    """
    name = file_name(tab)
    stem, dot, rest = name.partition(".")
    prefix = prefix or stem
    results = await tabix_tabix_cmd(tab=tab, tai=tai, regions=regions, prefix=prefix, args=args)
    if regions is not None:
        out = f"{prefix}{dot}{rest}"
        extracted = await results.get_file(out)
        if extracted is None:
            raise FileNotFoundError(f"tabix wrote no {out}")
        return None, extracted
    for ext in ("tbi", "csi"):
        index = await results.get_file(f"{name}.{ext}")
        if index is not None:
            return index, None
    raise FileNotFoundError(f"tabix wrote no {name}.tbi or {name}.csi")


def bgzip_output_name(name: str, prefix: str) -> tuple[str, bool]:
    """Upstream's tabix/bgzip output name for input ``name``, and whether bgzip decompresses it."""
    base, dot, ext = name.rpartition(".")
    if not dot:
        base, ext = name, ""
    if ext in ("gz", "bgz", "bgzf"):
        return f"{prefix}.{base.split('.')[-1]}", True
    return f"{prefix}.{ext}.gz", False


async def tabix_bgzip(infile: File, prefix: str, args: str = "") -> tuple[File, File | None]:
    """bgzip ``infile`` (or decompress it, if it's ``.gz``/``.bgz``/``.bgzf``); return the output and its ``.gzi``.

    The ``.gzi`` index is written with ``-i`` in ``args`` (named
    ``<output>.gzi`` unless ``args`` names it with ``-I``); otherwise None.
    """
    output, decompress = bgzip_output_name(file_name(infile), prefix)
    split_args = re.split(" +|=", args)
    indexed = "-i" in split_args or "--index" in split_args
    if indexed and "-I" not in split_args and "--index-name" not in split_args:
        args = f"{args} -I {output}.gzi"
    results = await tabix_bgzip_cmd(infile=infile, output=output, decompress=decompress, args=args)
    out = await results.get_file(output)
    if out is None:
        raise FileNotFoundError(f"bgzip wrote no {output}")
    gzi = None
    if indexed:
        files = [f for f in await results.list_files() if file_name(f).endswith(".gzi")]
        gzi = files[0] if files else None
    return out, gzi
