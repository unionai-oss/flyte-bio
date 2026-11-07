"""homer — HOMER tag directories, peak calling, peak annotation and browser tracks.

Exposes:

- :func:`homer_maketagdirectory` — a tag directory, ``<prefix>_tagdir``, from
  aligned reads (BAM/SAM with ``-format sam``, BED with ``-format bed``)
  against a genome FASTA: one ``<chromosome>.tags.tsv`` per chromosome with
  tags, plus ``tagInfo.txt``, ``tagAutocorrelation.txt``,
  ``tagCountDistribution.txt`` and ``tagLengthDistribution.txt``.
- :func:`homer_findpeaks` — peaks from a tag directory, optionally with a
  mappability (``-uniqmap``) directory: ``<prefix>.peaks.txt``.
- :func:`homer_makeucscfile` — a bedGraph track from a tag directory:
  ``<prefix>.bedGraph.gz``.
- :func:`homer_pos2bed` — HOMER peaks as BED: ``<prefix>.bed``.
- :func:`homer_annotatepeaks` — annotate peaks against a genome FASTA and a
  GTF (no preinstalled HOMER genome is needed): ``<prefix>.annotatePeaks.txt``,
  and ``<prefix>.annStats.txt`` when ``args`` asks for it with ``-annStats``.

HOMER records the names of its inputs in its outputs (``tagInfo.txt``, the
peak file header, the annotation header), so every input is linked into the
working directory under its own name first: reads keep the caller's order,
and a tag or mappability directory is linked under its ``Dir`` name. The tag
directory :func:`homer_maketagdirectory` returns is named ``<prefix>_tagdir``,
as upstream.

A tag directory holds one file per chromosome that has tags. That suits
genomes of a few dozen chromosomes; a genome with hundreds of contigs
(scaffolds, alt haplotypes) makes a ``Dir`` of hundreds of files, which is
slow to stage, so align such samples to a primary assembly first.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

HOMER_IMAGE = (
    "community.wave.seqera.io/library/bioconductor-deseq2_bioconductor-edger_homer_samtools_pruned:a8f4c58755bb281b"
)
# Upstream's process_medium asks for more; HOMER on a test-sized genome needs little.
# Raise with <cmd>.as_task().override(resources=...) for whole genomes.
DEFAULT_RESOURCES = flyte.Resources(cpu=2, memory="8Gi")

# Shared by the commands that take a tag directory (and findPeaks' uniqmap):
# link a directory into the working directory under NAME, file by file, so
# HOMER can write beside its files even though inputs are read-only.
LINK_DIR = r"""
        link_dir() {
            local SRC=$1 NAME=$2 F
            mkdir -p "$NAME"
            for F in "$SRC"/* "$SRC"/.[!.]*; do
                if [ -e "$F" ]; then ln -s "$F" "$NAME/$(basename "$F")"; fi
            done
        }
"""

# A list[File] is staged under the files' names and globbed alphabetically, so
# the names ("/"-separated) carry the caller's order into the script.
homer_maketagdirectory_cmd = shell.create(
    name="homer_maketagdirectory",
    image=HOMER_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"reads": list[File], "read_names": str, "fasta": File, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        READS=({inputs.reads})
        FA=({inputs.fasta})
        READ_NAMES={inputs.read_names}
        P={inputs.prefix}
        ARGS={inputs.args}
        OUT={outputs.results}
        W=$(mktemp -d); cd "$W"

        IFS=/ read -ra NAMES <<< "$READ_NAMES"
        READS_DIR=$(dirname "${READS[0]}")
        for NAME in "${NAMES[@]}"; do ln -s "$READS_DIR/$NAME" "$NAME"; done
        FA_NAME=$(basename "${FA[0]}"); ln -s "${FA[0]}" "$FA_NAME"

        makeTagDirectory \
            "${P}_tagdir" \
            -genome "$FA_NAME" \
            $ARGS \
            "${NAMES[@]}" || exit 1
        mv "${P}_tagdir" "$OUT"/
    """,
)

homer_findpeaks_cmd = shell.create(
    name="homer_findpeaks",
    image=HOMER_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "tagdir": Dir,
        "tagdir_name": str,
        "uniqmap": Dir | None,
        "uniqmap_name": str,
        "prefix": str,
        "args": str,
    },
    defaults={"uniqmap_name": "", "args": ""},
    outputs={"results": Dir},
    script=LINK_DIR
    + r"""
        TAGDIR={inputs.tagdir}
        TAGDIR_NAME={inputs.tagdir_name}
        UNIQMAP={inputs.uniqmap}
        UNIQMAP_NAME={inputs.uniqmap_name}
        P={inputs.prefix}
        ARGS={inputs.args}
        OUT={outputs.results}
        W=$(mktemp -d); cd "$W"

        link_dir "$TAGDIR" "$TAGDIR_NAME"
        UNIQMAP_FLAG=()
        if [ -n "$UNIQMAP" ] && [ -d "$UNIQMAP" ]; then
            ln -s "$UNIQMAP" "$UNIQMAP_NAME"
            UNIQMAP_FLAG=(-uniqmap "$UNIQMAP_NAME")
        fi

        findPeaks \
            "$TAGDIR_NAME" \
            $ARGS \
            -o "$P.peaks.txt" \
            "${UNIQMAP_FLAG[@]}" || exit 1
        mv "$P.peaks.txt" "$OUT"/
    """,
)

homer_makeucscfile_cmd = shell.create(
    name="homer_makeucscfile",
    image=HOMER_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"tagdir": Dir, "tagdir_name": str, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=LINK_DIR
    + r"""
        TAGDIR={inputs.tagdir}
        TAGDIR_NAME={inputs.tagdir_name}
        P={inputs.prefix}
        ARGS={inputs.args}
        OUT={outputs.results}
        W=$(mktemp -d); cd "$W"

        link_dir "$TAGDIR" "$TAGDIR_NAME"
        makeUCSCfile \
            "$TAGDIR_NAME" \
            -o "$P.bedGraph" \
            $ARGS || exit 1
        mv "$P.bedGraph.gz" "$OUT"/
    """,
)

homer_pos2bed_cmd = shell.create(
    name="homer_pos2bed",
    image=HOMER_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"peaks": File, "prefix": str},
    outputs={"results": Dir},
    script=r"""
        PEAKS=({inputs.peaks})
        P={inputs.prefix}
        OUT={outputs.results}
        W=$(mktemp -d); cd "$W"

        NAME=$(basename "${PEAKS[0]}"); ln -s "${PEAKS[0]}" "$NAME"
        pos2bed.pl "$NAME" > "$P.bed" || exit 1
        mv "$P.bed" "$OUT"/
    """,
)

homer_annotatepeaks_cmd = shell.create(
    name="homer_annotatepeaks",
    image=HOMER_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"peaks": File, "fasta": File, "gtf": File, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        PEAKS=({{inputs.peaks}}); FA=({{inputs.fasta}}); GTF=({{inputs.gtf}})
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"

        PEAKS_NAME=$(basename "${{PEAKS[0]}}"); ln -s "${{PEAKS[0]}}" "$PEAKS_NAME"
        FA_NAME=$(basename "${{FA[0]}}"); ln -s "${{FA[0]}}" "$FA_NAME"
        GTF_NAME=$(basename "${{GTF[0]}}"); ln -s "${{GTF[0]}}" "$GTF_NAME"

        annotatePeaks.pl \
            "$PEAKS_NAME" \
            "$FA_NAME" \
            $ARGS \
            -gtf "$GTF_NAME" \
            -cpu $TASK_CPUS \
            > "$P.annotatePeaks.txt" || exit 1
        shopt -s nullglob
        mv ./*annotatePeaks.txt ./*annStats.txt "$OUT"/
    """,
)

env = flyte.TaskEnvironment.from_task(
    "homer",
    homer_maketagdirectory_cmd.as_task(),
    homer_findpeaks_cmd.as_task(),
    homer_makeucscfile_cmd.as_task(),
    homer_pos2bed_cmd.as_task(),
    homer_annotatepeaks_cmd.as_task(),
)


def file_name(file: File) -> str:
    return file.path.rstrip("/").rsplit("/", 1)[-1]


def dir_name(directory: Dir) -> str:
    """The name HOMER sees for ``directory``: its ``Dir`` name, else its path's last part."""
    return directory.name or directory.path.rstrip("/").rsplit("/", 1)[-1]


async def pick(results: Dir, name: str, tool: str) -> File:
    found = await results.get_file(name)
    if found is None:
        raise FileNotFoundError(f"HOMER {tool} wrote no {name}")
    return found


@dataclass
class HomerTagDirectory:
    tagdir: Dir  # <prefix>_tagdir: per-chromosome <chr>.tags.tsv and HOMER's tag statistics
    taginfo: File  # <prefix>_tagdir/tagInfo.txt


async def homer_maketagdirectory(
    reads: list[File], fasta: File, prefix: str = "sample", args: str = ""
) -> HomerTagDirectory:
    """A tag directory ``<prefix>_tagdir`` of ``reads`` against ``fasta``.

    ``args`` carries HOMER's options, such as ``-format sam`` for BAMs or
    ``-format bed``. The reads keep their names and the caller's order, which
    ``tagInfo.txt`` records.
    """
    if not reads:
        raise ValueError("homer_maketagdirectory needs reads")
    names = [file_name(f) for f in reads]
    if len(set(names)) != len(names):
        raise ValueError(f"homer_maketagdirectory needs reads with distinct names, got {names}")
    results = await homer_maketagdirectory_cmd(
        reads=reads, read_names="/".join(names), fasta=fasta, prefix=prefix, args=args
    )
    tagdir = Dir.from_existing_remote(results.path.rstrip("/") + f"/{prefix}_tagdir")
    return HomerTagDirectory(tagdir=tagdir, taginfo=await pick(tagdir, "tagInfo.txt", "makeTagDirectory"))


async def homer_findpeaks(tagdir: Dir, uniqmap: Dir | None = None, prefix: str = "sample", args: str = "") -> File:
    """Peaks from ``tagdir``: ``<prefix>.peaks.txt``.

    ``uniqmap`` is an optional directory of HOMER mappability files
    (``-uniqmap``). Both directories are passed to HOMER under their ``Dir``
    names, which the peak file's header records. Upstream names the output
    ``<sample>-<uniqmap name without its last extension>`` when a uniqmap is
    given; pass that as ``prefix`` to match.
    """
    results = await homer_findpeaks_cmd(
        tagdir=tagdir,
        tagdir_name=dir_name(tagdir),
        uniqmap=uniqmap,
        uniqmap_name=dir_name(uniqmap) if uniqmap is not None else "",
        prefix=prefix,
        args=args,
    )
    return await pick(results, f"{prefix}.peaks.txt", "findPeaks")


async def homer_makeucscfile(tagdir: Dir, prefix: str = "sample", args: str = "") -> File:
    """A bedGraph track of ``tagdir``: ``<prefix>.bedGraph.gz``."""
    results = await homer_makeucscfile_cmd(tagdir=tagdir, tagdir_name=dir_name(tagdir), prefix=prefix, args=args)
    return await pick(results, f"{prefix}.bedGraph.gz", "makeUCSCfile")


async def homer_pos2bed(peaks: File, prefix: str = "sample") -> File:
    """HOMER ``peaks`` (e.g. from :func:`homer_findpeaks`) as ``<prefix>.bed``."""
    results = await homer_pos2bed_cmd(peaks=peaks, prefix=prefix)
    return await pick(results, f"{prefix}.bed", "pos2bed.pl")


@dataclass
class HomerAnnotatePeaksResult:
    txt: File  # <prefix>.annotatePeaks.txt
    stats: File | None  # the -annStats table, when args asks for one
    results: Dir  # everything annotatePeaks.pl wrote


async def homer_annotatepeaks(
    peaks: File, fasta: File, gtf: File, prefix: str = "sample", args: str = ""
) -> HomerAnnotatePeaksResult:
    """Annotate ``peaks`` against ``fasta`` and ``gtf``: ``<prefix>.annotatePeaks.txt``.

    Upstream's pipelines add ``-annStats <prefix>.annStats.txt`` (and options
    such as ``-gid``) through ``args``.
    """
    results = await homer_annotatepeaks_cmd(peaks=peaks, fasta=fasta, gtf=gtf, prefix=prefix, args=args)
    files = {file_name(f): f for f in await results.list_files()}
    stats = [f for name, f in sorted(files.items()) if name.endswith("annStats.txt")]
    return HomerAnnotatePeaksResult(
        txt=await pick(results, f"{prefix}.annotatePeaks.txt", "annotatePeaks.pl"),
        stats=stats[0] if stats else None,
        results=results,
    )
