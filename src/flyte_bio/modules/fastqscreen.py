"""fastqscreen — screen reads against several genomes for contamination.

Exposes:

- :func:`fastqscreen_database` — a FastQ Screen database from Bowtie 2
  indexes: a Dir holding each index as ``dir1/``, ``dir2/``, ... and a
  ``fastq_screen.conf`` naming each genome, as upstream's buildfromindex
  module lays it out. The conf's paths are relative
  (``FastQ_Screen_Genomes/dir1/<index>``), so the database is used from a
  directory named ``FastQ_Screen_Genomes``, as upstream stages it.
- :func:`fastqscreen` — map a sample's reads against every genome in a
  database, writing ``<reads>_screen.{txt,png,html}`` and, with ``--tag``
  (and ``--filter``) in ``args``, ``<reads>.tagged.fastq.gz`` (and
  ``<reads>.tagged_filter.fastq.gz``), named after each FASTQ.

Shell tasks can't take a ``list[Dir]``, so :func:`fastqscreen_database`
copies the indexes through the caller's pod, as upstream's module copies
them into the database folder.
"""

import asyncio
import tempfile
from dataclasses import dataclass
from pathlib import Path

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

FASTQSCREEN_IMAGE = "community.wave.seqera.io/library/fastq-screen_perl-gdgraph:5c1786a5d5bc1309"

# Upstream's process_medium, capped at the 4 CPUs its tests run with.
DEFAULT_RESOURCES = flyte.Resources(cpu=4, memory="8Gi")

# The FASTQs are linked into in/ and the reports written to the output Dir
# (FastQ Screen's --outdir), named after each FASTQ.
fastqscreen_cmd = shell.create(
    name="fastqscreen",
    image=FASTQSCREEN_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"reads": list[File], "database": Dir, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        READS=({{inputs.reads}})
        ARGS={{inputs.args}}
        W=$(mktemp -d); mkdir "$W/in"; cd "$W"
        ln -s {{inputs.database}} FastQ_Screen_Genomes
        NAMES=()
        for r in "${{READS[@]}}"; do ln -s "$r" in/; NAMES+=("in/$(basename "$r")"); done
        fastq_screen --threads $TASK_CPUS \
            --conf FastQ_Screen_Genomes/fastq_screen.conf \
            --outdir {{outputs.results}} \
            "${{NAMES[@]}}" \
            $ARGS || exit 1
    """,
)

env = flyte.TaskEnvironment.from_task("fastqscreen", fastqscreen_cmd.as_task())


async def fastqscreen_database(genome_names: list[str], indexes: list[Dir]) -> Dir:
    """A FastQ Screen database of the Bowtie 2 ``indexes``, the i-th named ``genome_names[i]``.

    Each index is found, as upstream, by its ``<name>.rev.1.bt2`` file.
    """
    if len(genome_names) != len(indexes):
        raise ValueError(f"{len(genome_names)} genome names for {len(indexes)} indexes")
    root = Path(tempfile.mkdtemp(prefix="fastqscreen_")) / "FastQ_Screen_Genomes"
    root.mkdir()

    async def copy(i: int, index: Dir) -> str:
        folder = f"dir{i + 1}"
        await index.download(str(root / folder))
        found = sorted(p.name for p in (root / folder).glob("*.rev.1.bt2"))
        if not found:
            raise FileNotFoundError(f"index {i + 1} ({genome_names[i]}) has no *.rev.1.bt2 file")
        return f"{folder}/{found[0].removesuffix('.rev.1.bt2')}"

    paths = await asyncio.gather(*(copy(i, index) for i, index in enumerate(indexes)))
    conf = "\n\n".join(
        f"########## {name} \nDATABASE {name} FastQ_Screen_Genomes/{path}"
        for name, path in zip(genome_names, paths, strict=True)
    )
    (root / "fastq_screen.conf").write_text(conf + "\n")
    return await Dir.from_local(str(root))


@dataclass
class FastqScreenResult:
    txt: list[File]  # <reads>_screen.txt, one per FASTQ
    png: list[File]  # <reads>_screen.png
    html: list[File]  # <reads>_screen.html
    fastq: list[File]  # with --tag: <reads>.tagged.fastq.gz (and with --filter, <reads>.tagged_filter.fastq.gz)
    results: Dir  # everything the task wrote


async def fastqscreen(reads: list[File], database: Dir, args: str = "") -> FastqScreenResult:
    """Screen ``reads`` against every genome in ``database`` (from :func:`fastqscreen_database`)."""
    results = await fastqscreen_cmd(reads=reads, database=database, args=args)
    files = sorted([f async for f in results.walk()], key=lambda f: f.name)

    def ending(suffix: str) -> list[File]:
        return [f for f in files if f.name.endswith(suffix)]

    result = FastqScreenResult(
        txt=ending("_screen.txt"),
        png=ending("_screen.png"),
        html=ending("_screen.html"),
        fastq=ending(".fastq.gz"),
        results=results,
    )
    if not result.txt or not result.html:
        raise FileNotFoundError("fastq_screen wrote no _screen.txt / _screen.html report")
    return result
