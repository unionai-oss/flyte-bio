"""cat — file concatenation utilities.

Exposes:

- :func:`cat_fastq` — concatenate gzipped FASTQ shards into a single
  merged file, preserving input order.
- :func:`cat_cat` — concatenate any files into one, decompressing or
  compressing by the input and output extensions (see below).

Inputs are assumed already gzipped (.fastq.gz) and are joined with ``cat``
byte-for-byte without recompression (concatenated gzip members are a valid
gzip stream) — which keeps the output MD5 stable and snapshot-comparable.

Order matters: a sample's R1 and R2 shards must be merged in the same order
or read pairs stop lining up. A ``list[File]`` input can't carry that order
into the container — list items are staged under their original basenames
and a glob lists them alphabetically — so :func:`cat_fastq` folds the list
pairwise through :data:`cat_fastq_cmd`, whose two inputs are explicit.

The merged file is written as ``<prefix>.fastq.gz`` inside an output Dir
(as upstream names it ``<sample>.merged.fastq.gz``): tools such as
``fq`` decide whether to decompress from the extension, which a bare
output path wouldn't have.

:func:`cat_cat` follows the upstream module's last working version (later
upstream versions only fail with a deprecation notice). The output is named
``<filename>``, by default ``<sample>`` plus the first input's extension (with
a ``.gz``, the extension before it too: ``.fasta.gz``). The tool depends on
whether the first input and the output end in ``.gz``:

- both or neither: ``cat``;
- gzipped input, plain output: ``zcat``;
- plain input, gzipped output: ``cat | pigz``, on the task's CPUs.

As upstream, the output may not share a name with an input. The files are
joined in list order; since a ``list[File]`` is staged under its basenames
and globbed alphabetically, the wrapper passes the order as a list of names,
so the basenames must be distinct (upstream's staging requires that too).

:func:`concatenate` is upstream's replacement for cat/cat (its find/concatenate
module), in upstream's findutils + pigz image. Unlike :func:`cat_cat` it joins
the files in name order, as upstream does (``find | sort -z``), whatever the
list order. The output is named ``<filename>``, by default ``<sample>`` plus
the first input's extension (the last one, or the last two for a ``.gz``:
``.maf.gz``). The inputs must be all gzipped or all plain, and:

- gzipped or plain to the same: ``cat``;
- gzipped to plain output: ``pigz -cd``, on the task's CPUs;
- plain to a ``.gz`` output: ``cat``, then ``pigz`` in place, on the task's CPUs.
"""

import os
import re

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

# Seqera wave coreutils image.
CAT_IMAGE = "community.wave.seqera.io/library/coreutils_grep_gzip_lbzip2_pruned:838ba80435a629f8"

# Upstream's pigz container for the general concatenation.
PIGZ_IMAGE = "quay.io/biocontainers/pigz:2.8"

# Upstream's findutils + pigz container for find/concatenate.
FIND_PIGZ_IMAGE = "community.wave.seqera.io/library/findutils_pigz:c4dd5edc44402661"

# Single-threaded; fixed modest resources.
DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory="6Gi")
# Upstream's process_low gives pigz 2 CPUs; cat needs little memory.
CAT_CAT_RESOURCES = flyte.Resources(cpu=2, memory="4Gi")
CONCATENATE_RESOURCES = flyte.Resources(cpu=2, memory="4Gi")


cat_fastq_cmd = shell.create(
    name="cat_fastq",
    image=CAT_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"first": File, "second": File, "filename": str},
    outputs={"results": Dir},
    script=r"""
        cat {inputs.first} {inputs.second} > {outputs.results}/{inputs.filename}
    """,
)


# The names are joined with "/", which can't occur in a basename.
cat_cat_cmd = shell.create(
    name="cat_cat",
    image=PIGZ_IMAGE,
    resources=CAT_CAT_RESOURCES,
    inputs={"files": list[File], "names": str, "filename": str, "args": str, "args2": str},
    defaults={"args": "", "args2": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        F=({{inputs.files}})
        DIR=$(dirname "${{F[0]}}")
        IFS=/ read -ra NAMES <<< {{inputs.names}}
        IN=()
        for n in "${{NAMES[@]}}"; do IN+=("$DIR/$n"); done
        P={{inputs.filename}}
        ARGS={{inputs.args}}
        ARGS2={{inputs.args2}}
        OUT={{outputs.results}}/"$P"
        IN_ZIP=false; OUT_ZIP=false
        if [[ "${{NAMES[0]}}" == *.gz ]]; then IN_ZIP=true; fi
        if [[ "$P" == *.gz ]]; then OUT_ZIP=true; fi
        if [ $IN_ZIP = true ] && [ $OUT_ZIP = false ]; then
            zcat $ARGS "${{IN[@]}}" > "$OUT"
        elif [ $IN_ZIP = false ] && [ $OUT_ZIP = true ]; then
            cat $ARGS "${{IN[@]}}" | pigz -c -p $TASK_CPUS $ARGS2 > "$OUT"
        else
            cat $ARGS "${{IN[@]}}" > "$OUT"
        fi
    """,
)


# Upstream's loop, over the inputs linked into a `to_concatenate/` directory as upstream stages them.
# OUT_NAME is the file the loop appends to: for a plain-to-gzipped join, the output without `.gz`,
# which pigz then compresses in place.
concatenate_cmd = shell.create(
    name="concatenate",
    image=FIND_PIGZ_IMAGE,
    resources=CONCATENATE_RESOURCES,
    inputs={"files": list[File], "filename": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        F=({{inputs.files}})
        W=$(mktemp -d)
        mkdir "$W/to_concatenate"
        for f in "${{F[@]}}"; do ln -s "$f" "$W/to_concatenate/"; done
        P={{inputs.filename}}
        ARGS={{inputs.args}}
        IN_ZIP=false; OUT_ZIP=false
        if [[ "${{F[0]}}" == *.gz ]]; then IN_ZIP=true; fi
        if [[ "$P" == *.gz ]]; then OUT_ZIP=true; fi
        OUT_NAME="$P"
        if [ $IN_ZIP = false ] && [ $OUT_ZIP = true ]; then OUT_NAME="${{P%.gz}}"; fi
        OUT={{outputs.results}}/"$OUT_NAME"
        CMD1=cat
        if [ $IN_ZIP = true ] && [ $OUT_ZIP = false ]; then CMD1="pigz -cd -p $TASK_CPUS"; fi
        while IFS= read -r -d $'\0' file; do
            $CMD1 "$file" >> "$OUT"
        done < <( find "$W/to_concatenate/" -mindepth 1 -print0 | sort -z )
        if [ $IN_ZIP = false ] && [ $OUT_ZIP = true ]; then pigz -p $TASK_CPUS $ARGS "$OUT"; fi
        true
    """,
)


# One image per environment: cat_cat and concatenate run in upstream's images, so `env` aggregates three.
cat_fastq_env = flyte.TaskEnvironment.from_task("cat_fastq", cat_fastq_cmd.as_task())
cat_cat_env = flyte.TaskEnvironment.from_task("cat_cat", cat_cat_cmd.as_task())
concatenate_env = flyte.TaskEnvironment.from_task("concatenate", concatenate_cmd.as_task())
env = flyte.TaskEnvironment(name="cat", depends_on=[cat_fastq_env, cat_cat_env, concatenate_env])


async def cat_fastq(reads: list[File], prefix: str = "merged") -> File:
    """Concatenate ``reads`` in list order into ``<prefix>.fastq.gz``.

    A single shard is returned as-is (no copy), matching upstream, which
    only merges samples that have more than one run.
    """
    if not reads:
        raise ValueError("cat_fastq needs at least one file")
    filename = f"{prefix}.fastq.gz"
    merged = reads[0]
    for shard in reads[1:]:
        results = await cat_fastq_cmd(first=merged, second=shard, filename=filename)
        found = await results.get_file(filename)
        if found is None:
            raise FileNotFoundError(f"cat_fastq wrote no {filename}")
        merged = found
    return merged


def file_suffix(filename: str) -> str:
    """A file's extension as upstream takes it: for a ``.gz``, with the one before it (``.fasta.gz``)."""
    match = re.match(r"^.*?((\.\w{1,5})?(\.\w{1,5}\.gz$))", filename)
    if match:
        return match.group(1)
    dot = filename.rfind(".")
    return filename[dot:] if dot >= 0 else ""


async def cat_cat(
    files: list[File],
    sample: str = "sample",
    filename: str | None = None,
    args: str = "",
    args2: str = "",
) -> File:
    """Concatenate ``files`` in list order into ``filename``.

    ``filename`` defaults to ``<sample>`` plus the first file's extension.
    ``args`` go to ``cat``/``zcat`` and ``args2`` to ``pigz``.
    """
    if not files:
        raise ValueError("cat_cat needs at least one file")
    names = [f.name or os.path.basename(f.path) for f in files]
    if len(set(names)) != len(names):
        raise ValueError(f"cat_cat needs files with distinct names, got {names}")
    if filename is None:
        filename = f"{sample}{file_suffix(names[0])}"
    if filename.strip() in names:
        raise ValueError(
            f"The name of the input file can't be the same as for the output prefix (currently `{filename}`). "
            "Please choose a different one."
        )
    results = await cat_cat_cmd(files=files, names="/".join(names), filename=filename, args=args, args2=args2)
    found = await results.get_file(filename)
    if found is None:
        raise FileNotFoundError(f"cat_cat wrote no {filename}")
    return found


def trailing_extensions(filename: str) -> str:
    """A file's extension as find/concatenate takes it: the last one, or for a ``.gz`` the last two."""
    base = filename
    for _ in range(2 if filename.endswith(".gz") else 1):
        dot = base.rfind(".")
        if dot < 0:
            break
        base = base[:dot]
    return filename[len(base) :]


async def concatenate(
    files: list[File],
    sample: str = "sample",
    filename: str | None = None,
    args: str = "",
) -> File:
    """Concatenate ``files``, in name order, into ``filename``.

    ``filename`` defaults to ``<sample>`` plus the first file's extension; a
    ``.gz`` name gives a gzipped output. ``args`` go to ``pigz`` when it
    compresses plain inputs.
    """
    if not files:
        raise ValueError("concatenate needs at least one file")
    names = [f.name or os.path.basename(f.path) for f in files]
    if len(set(names)) != len(names):
        raise ValueError(f"concatenate needs files with distinct names, got {names}")
    zipped = [n.endswith(".gz") for n in names]
    if any(zipped) and not all(zipped):
        raise ValueError(
            "All files provided to this module must either be gzipped (and have the .gz extension) or unzipped "
            "(and not have the .gz extension). A mix of both is not allowed."
        )
    if filename is None:
        filename = f"{sample}{trailing_extensions(names[0])}"
    results = await concatenate_cmd(files=files, filename=filename, args=args)
    found = await results.get_file(filename)
    if found is None:
        raise FileNotFoundError(f"concatenate wrote no {filename}")
    return found
