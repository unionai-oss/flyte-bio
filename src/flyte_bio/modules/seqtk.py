"""seqtk — subsample and quality-trim FASTQ reads.

Exposes :func:`seqtk_sample`, which samples ``sample_size`` reads (a count,
or a fraction below 1) from each FASTQ file with ``seqtk sample``. As
upstream:

- each input ``<name>`` gives ``<prefix>_<name>``, with a ``.fq.gz`` suffix
  spelled ``.fastq.gz``, always gzipped (``gzip --no-name``);
- a seed is added (``-s100``) unless ``args`` sets one (`` -s<n>``), so
  paired files sampled with one seed stay paired.

and :func:`seqtk_trim`, which quality-trims each FASTQ file with ``seqtk
trimfq`` (its default: Phred error-probability trimming at 0.05; ``args``
passes ``-q``, ``-b``, ``-e`` and the rest). As upstream, each input
``<name>`` gives ``<prefix>_<name>``, gzipped (``gzip --no-name``).
"""

import re

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

SEQTK_IMAGE = "community.wave.seqera.io/library/seqtk:1.4--3fa479d3961a5b8d"
# Upstream's process_single.
DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory="2Gi")

# Upstream's test for a seed already in ``args``.
HAS_SEED = re.compile(r".* -s ?[0-9]+.*")

# A list[File] is staged under the files' names and globbed alphabetically, so
# ``names`` carries the caller's order (one name per line) into the script.
seqtk_sample_cmd = shell.create(
    name="seqtk_sample",
    image=SEQTK_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"reads": list[File], "names": str, "sample_size": str, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        set -o pipefail
        IN=({inputs.reads})
        NAMES={inputs.names}
        SIZE={inputs.sample_size}
        P={inputs.prefix}
        ARGS={inputs.args}
        OUT={outputs.results}
        set -f
        D=$(dirname "${IN[0]}")
        while IFS= read -r NAME; do
            NEW=$(printf '%s\n' "$NAME" | sed 's/\.f\(ast\)\?q\.gz$/.fastq.gz/')
            seqtk sample $ARGS "$D/$NAME" "$SIZE" | gzip --no-name > "$OUT/${P}_$NEW" || exit 1
        done <<< "$NAMES"
    """,
)

# Upstream's process_low asks for 2 CPUs and 12 GB; seqtk trimfq is
# single-threaded and streams, so 2Gi is plenty.
seqtk_trim_cmd = shell.create(
    name="seqtk_trim",
    image=SEQTK_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"reads": list[File], "names": str, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        set -o pipefail
        IN=({inputs.reads})
        NAMES={inputs.names}
        P={inputs.prefix}
        ARGS={inputs.args}
        OUT={outputs.results}
        D=$(dirname "${IN[0]}")
        while IFS= read -r NAME; do
            seqtk trimfq $ARGS "$D/$NAME" | gzip --no-name > "$OUT/${P}_$NAME" || exit 1
        done <<< "$NAMES"
    """,
)

env = flyte.TaskEnvironment.from_task("seqtk", seqtk_sample_cmd.as_task(), seqtk_trim_cmd.as_task())


def file_name(file: File) -> str:
    return file.path.rstrip("/").rsplit("/", 1)[-1]


def sampled_name(name: str) -> str:
    return re.sub(r"\.f(ast)?q\.gz$", ".fastq.gz", name)


async def seqtk_sample(
    reads: list[File], sample_size: int | float, prefix: str = "sample", args: str = ""
) -> list[File]:
    """``sample_size`` reads sampled from each of ``reads``: ``<prefix>_<name>``, in list order."""
    if not sample_size:
        raise ValueError("SEQTK/SAMPLE must have a sample_size value included")
    if not reads:
        raise ValueError("seqtk_sample needs at least one file")
    names = [file_name(f) for f in reads]
    if len(set(names)) != len(names):
        raise ValueError(f"seqtk_sample inputs need distinct names, got {names}")
    if not HAS_SEED.fullmatch(args):
        args += " -s100"
    results = await seqtk_sample_cmd(
        reads=reads, names="\n".join(names), sample_size=str(sample_size), prefix=prefix, args=args.strip()
    )
    outputs = []
    for name in names:
        out = f"{prefix}_{sampled_name(name)}"
        found = await results.get_file(out)
        if found is None:
            raise FileNotFoundError(f"seqtk sample wrote no {out}")
        outputs.append(found)
    return outputs


async def seqtk_trim(reads: list[File], prefix: str = "sample", args: str = "") -> list[File]:
    """Each of ``reads`` quality-trimmed with ``seqtk trimfq``: ``<prefix>_<name>``, gzipped, in list order."""
    if not reads:
        raise ValueError("seqtk_trim needs at least one file")
    names = [file_name(f) for f in reads]
    if len(set(names)) != len(names):
        raise ValueError(f"seqtk_trim inputs need distinct names, got {names}")
    results = await seqtk_trim_cmd(reads=reads, names="\n".join(names), prefix=prefix, args=args)
    outputs = []
    for name in names:
        found = await results.get_file(f"{prefix}_{name}")
        if found is None:
            raise FileNotFoundError(f"seqtk trimfq wrote no {prefix}_{name}")
        outputs.append(found)
    return outputs
