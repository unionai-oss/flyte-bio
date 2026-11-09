"""gapseq — metabolic pathway prediction and genome-scale metabolic model reconstruction for bacteria.

Exposes:

- :func:`gapseq_requestdb` — download gapseq's reference sequences for a
  taxon (``gapseq update-sequences -t <taxon>``) as ``gapseq_db.tar.gz``, an
  archive of the ``gapseq_db`` folder (thousands of FASTA files, too many for
  a Dir to stage reliably). :func:`gapseq_find` and :func:`gapseq_doall`
  extract it and pass it as ``-D``; without it they use the sequences the
  container ships.
- :func:`gapseq_find` — predict reactions and pathways (``gapseq find -p all
  -u all -b 200 -t Bacteria``): ``<name>-all-Reactions.tbl``,
  ``<name>-all-Pathways.tbl``, plus any ``.fna`` and ``.log`` files, where
  ``<name>`` is the FASTA's name without its extension.
- :func:`gapseq_findtransport` — predict transporters (``gapseq
  find-transport``): ``<name>-Transporter.tbl``, plus any ``.fna`` and
  ``.log`` files.
- :func:`gapseq_draft` — build a draft model from those tables (``gapseq
  draft``): ``<prefix>-draft.RDS`` and its SBML ``<prefix>-draft.xml``.
- :func:`gapseq_medium` — predict a growth medium for a draft model
  (``gapseq medium``): ``<prefix>-medium.csv``.
- :func:`gapseq_fill` — gap-fill a draft model on a medium (``gapseq fill``),
  by default gapseq's ``ALLmed.csv``: ``<prefix>-filled.RDS`` and
  ``<prefix>-filled.xml``. As upstream, the draft is staged under its own
  name, which must end in ``-draft.RDS``.
- :func:`gapseq_doall` — the whole reconstruction in one step (``gapseq
  doall -t Bacteria``), optionally on a given medium: the models
  (``<name>-draft.RDS``, ``<name>.RDS`` and their ``.xml``), the ``.tbl``
  tables and any ``.fna`` and ``.log`` files.

``find``, ``find-transport`` and ``doall`` get the task's CPUs as threads
(``-K``). gapseq solves its linear programs with GLPK (open source) through
cobrar, so no solver licence is needed. gapseq is AGPL-3.0.
"""

import os
from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

GAPSEQ_IMAGE = "community.wave.seqera.io/library/gapseq:2.1.0--c32b876ebb5e5f5b"
# Upstream's process_single, process_medium and process_high, scaled down.
SINGLE_RESOURCES = flyte.Resources(cpu=1, memory="4Gi", disk="20Gi")
MEDIUM_RESOURCES = flyte.Resources(cpu=6, memory="16Gi", disk="20Gi")
HIGH_RESOURCES = flyte.Resources(cpu=8, memory="24Gi", disk="20Gi")

# Extracts the optional database archive into $W/gapseq_db and sets DB_ARG.
EXTRACT_DB = r"""
        DB_ARG=()
        if [ ${#DB[@]} -gt 0 ]; then
            mkdir -p "$W/db"
            tar -xzf "${DB[0]}" -C "$W/db"
            DB_ARG=(-D "$W/db/gapseq_db")
        fi
"""

gapseq_requestdb_cmd = shell.create(
    name="gapseq_requestdb",
    image=GAPSEQ_IMAGE,
    resources=SINGLE_RESOURCES,
    inputs={"taxon": str},
    outputs={"results": Dir},
    script=r"""
        TAXON={inputs.taxon}
        OUT={outputs.results}
        W=$(mktemp -d); cd "$W"
        export HOME="$W"
        # gapseq resolves -D with readlink -f, so it writes here, not under $HOME.
        gapseq update-sequences -t "$TAXON" -D gapseq_db || exit 1
        tar -czf "$OUT/gapseq_db.tar.gz" gapseq_db
    """,
)

gapseq_find_cmd = shell.create(
    name="gapseq_find",
    image=GAPSEQ_IMAGE,
    resources=MEDIUM_RESOURCES,
    inputs={"fasta": File, "db": File | None, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        FA=({{inputs.fasta}})
        shopt -s nullglob; DB=({{inputs.db}}); shopt -u nullglob
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        export HOME="$W"
        {EXTRACT_DB}
        ln -s "${{FA[0]}}" "$(basename "${{FA[0]}}")"
        gapseq find -p all -u all -b 200 -t Bacteria -K $TASK_CPUS "${{DB_ARG[@]}}" $ARGS \
            "$(basename "${{FA[0]}}")" || exit 1
        shopt -s nullglob
        mv ./*-all-Reactions.tbl ./*-all-Pathways.tbl ./*.fna ./*.log "$OUT"/
    """,
)

gapseq_findtransport_cmd = shell.create(
    name="gapseq_findtransport",
    image=GAPSEQ_IMAGE,
    resources=MEDIUM_RESOURCES,
    inputs={"fasta": File, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        FA=({{inputs.fasta}})
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        export HOME="$W"
        export LANG=C LC_ALL=C TZ=UTC R_LIBS_USER=
        ln -s "${{FA[0]}}" "$(basename "${{FA[0]}}")"
        gapseq find-transport -K $TASK_CPUS $ARGS "$(basename "${{FA[0]}}")" || exit 1
        shopt -s nullglob
        mv ./*.tbl ./*.fna ./*.log "$OUT"/
    """,
)

gapseq_draft_cmd = shell.create(
    name="gapseq_draft",
    image=GAPSEQ_IMAGE,
    resources=MEDIUM_RESOURCES,
    inputs={"reactions": File, "transporters": File | None, "pathways": File, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        RXN=({inputs.reactions}); PWY=({inputs.pathways})
        shopt -s nullglob; TRANSPORT=({inputs.transporters}); shopt -u nullglob
        P={inputs.prefix}
        ARGS={inputs.args}
        OUT={outputs.results}
        W=$(mktemp -d); cd "$W"
        export HOME="$W"
        TRANSPORT_ARG=()
        if [ ${#TRANSPORT[@]} -gt 0 ]; then TRANSPORT_ARG=(-t "${TRANSPORT[0]}"); fi
        gapseq draft -r "${RXN[0]}" "${TRANSPORT_ARG[@]}" -p "${PWY[0]}" -n "$P" $ARGS || exit 1
        shopt -s nullglob
        mv ./*.RDS ./*.xml ./*.log "$OUT"/
    """,
)

gapseq_medium_cmd = shell.create(
    name="gapseq_medium",
    image=GAPSEQ_IMAGE,
    resources=MEDIUM_RESOURCES,
    inputs={"draft": File, "pathways": File, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        DRAFT=({inputs.draft}); PWY=({inputs.pathways})
        P={inputs.prefix}
        ARGS={inputs.args}
        OUT={outputs.results}
        W=$(mktemp -d); cd "$W"
        export HOME="$W"
        gapseq medium --model "${DRAFT[0]}" --pathway.pred "${PWY[0]}" $ARGS --output.file "$P-medium.csv" || exit 1
        shopt -s nullglob
        mv ./*.csv ./*.log "$OUT"/
    """,
)

gapseq_fill_cmd = shell.create(
    name="gapseq_fill",
    image=GAPSEQ_IMAGE,
    resources=HIGH_RESOURCES,
    inputs={"draft": File, "medium": File | None, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        DRAFT=({inputs.draft})
        shopt -s nullglob; MEDIUM=({inputs.medium}); shopt -u nullglob
        P={inputs.prefix}
        ARGS={inputs.args}
        OUT={outputs.results}
        W=$(mktemp -d); cd "$W"
        export HOME="$W"
        GAPSEQ_ROOT=$(dirname "$(readlink -f "$(which gapseq)")")
        MEDIUM_FILE="$GAPSEQ_ROOT/dat/media/ALLmed.csv"
        if [ ${#MEDIUM[@]} -gt 0 ]; then MEDIUM_FILE="${MEDIUM[0]}"; fi
        ln -s "${DRAFT[0]}" "$(basename "${DRAFT[0]}")"
        DRAFT_MODEL=$(ls ./*-draft.RDS)
        gapseq fill -m "$DRAFT_MODEL" -n "$MEDIUM_FILE" $ARGS || exit 1
        FILLED_MODEL=$(ls ./*.RDS | grep -Ev '(-draft)\.RDS' | head -n1)
        mv "$FILLED_MODEL" "$OUT/$P-filled.RDS"
        FILLED_XML=$(ls ./*.xml | grep -Ev '(-draft)\.xml' | head -n1)
        mv "$FILLED_XML" "$OUT/$P-filled.xml"
        shopt -s nullglob
        LOGS=(./*.log)
        if [ ${#LOGS[@]} -gt 0 ]; then mv "${LOGS[@]}" "$OUT"/; fi
    """,
)

gapseq_doall_cmd = shell.create(
    name="gapseq_doall",
    image=GAPSEQ_IMAGE,
    resources=HIGH_RESOURCES,
    inputs={"fasta": File, "medium": File | None, "db": File | None, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        FA=({{inputs.fasta}})
        shopt -s nullglob; MEDIUM=({{inputs.medium}}); DB=({{inputs.db}}); shopt -u nullglob
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        export HOME="$W"
        {EXTRACT_DB}
        MEDIUM_ARG=()
        if [ ${{#MEDIUM[@]}} -gt 0 ]; then MEDIUM_ARG=(-m "${{MEDIUM[0]}}"); fi
        ln -s "${{FA[0]}}" "$(basename "${{FA[0]}}")"
        gapseq doall -t Bacteria "${{MEDIUM_ARG[@]}}" -K $TASK_CPUS "${{DB_ARG[@]}}" $ARGS \
            "$(basename "${{FA[0]}}")" || exit 1
        shopt -s nullglob
        mv ./*.RDS ./*.xml ./*.tbl ./*.fna ./*.log "$OUT"/
    """,
)

env = flyte.TaskEnvironment.from_task(
    "gapseq",
    gapseq_requestdb_cmd.as_task(),
    gapseq_find_cmd.as_task(),
    gapseq_findtransport_cmd.as_task(),
    gapseq_draft_cmd.as_task(),
    gapseq_medium_cmd.as_task(),
    gapseq_fill_cmd.as_task(),
    gapseq_doall_cmd.as_task(),
)


def file_name(file: File) -> str:
    return file.name or os.path.basename(file.path.rstrip("/"))


async def files_ending(results: Dir, suffix: str) -> list[File]:
    """The files in ``results`` whose names end in ``suffix``, sorted by name."""
    files = [f for f in await results.list_files() if file_name(f).endswith(suffix)]
    return sorted(files, key=file_name)


async def one_ending(results: Dir, suffix: str, step: str) -> File:
    found = await files_ending(results, suffix)
    if len(found) != 1:
        raise FileNotFoundError(f"gapseq {step} wrote {len(found)} files ending in {suffix}, expected 1")
    return found[0]


async def gapseq_requestdb(taxon: str = "Bacteria") -> File:
    """gapseq's reference sequences for ``taxon``: ``gapseq_db.tar.gz``."""
    results = await gapseq_requestdb_cmd(taxon=taxon)
    return await one_ending(results, "gapseq_db.tar.gz", "update-sequences")


@dataclass
class GapseqFindResult:
    reactions: File  # <name>-all-Reactions.tbl
    pathways: File  # <name>-all-Pathways.tbl
    fna: list[File]
    log: list[File]
    results: Dir


async def gapseq_find(fasta: File, db: File | None = None, args: str = "") -> GapseqFindResult:
    """Reactions and pathways predicted for the genome or proteome ``fasta``."""
    results = await gapseq_find_cmd(fasta=fasta, db=db, args=args)
    return GapseqFindResult(
        reactions=await one_ending(results, "-all-Reactions.tbl", "find"),
        pathways=await one_ending(results, "-all-Pathways.tbl", "find"),
        fna=await files_ending(results, ".fna"),
        log=await files_ending(results, ".log"),
        results=results,
    )


@dataclass
class GapseqFindTransportResult:
    tbl: File  # <name>-Transporter.tbl
    fna: list[File]
    log: list[File]
    results: Dir


async def gapseq_findtransport(fasta: File, args: str = "") -> GapseqFindTransportResult:
    """Transporters predicted for the genome or proteome ``fasta``."""
    results = await gapseq_findtransport_cmd(fasta=fasta, args=args)
    return GapseqFindTransportResult(
        tbl=await one_ending(results, ".tbl", "find-transport"),
        fna=await files_ending(results, ".fna"),
        log=await files_ending(results, ".log"),
        results=results,
    )


@dataclass
class GapseqModelResult:
    model: File  # .RDS
    xml: File  # the SBML model
    log: list[File]
    results: Dir


async def gapseq_draft(
    reactions: File,
    pathways: File,
    transporters: File | None = None,
    prefix: str = "sample",
    args: str = "",
) -> GapseqModelResult:
    """A draft model from ``find``'s reactions and pathways and ``find-transport``'s transporters."""
    results = await gapseq_draft_cmd(
        reactions=reactions, transporters=transporters, pathways=pathways, prefix=prefix, args=args
    )
    return GapseqModelResult(
        model=await one_ending(results, f"{prefix}-draft.RDS", "draft"),
        xml=await one_ending(results, f"{prefix}-draft.xml", "draft"),
        log=await files_ending(results, ".log"),
        results=results,
    )


async def gapseq_medium(draft: File, pathways: File, prefix: str = "sample", args: str = "") -> File:
    """The growth medium predicted for ``draft``: ``<prefix>-medium.csv``."""
    results = await gapseq_medium_cmd(draft=draft, pathways=pathways, prefix=prefix, args=args)
    return await one_ending(results, f"{prefix}-medium.csv", "medium")


async def gapseq_fill(
    draft: File, medium: File | None = None, prefix: str = "sample", args: str = ""
) -> GapseqModelResult:
    """``draft`` (named ``*-draft.RDS``) gap-filled on ``medium`` (default: gapseq's ``ALLmed.csv``)."""
    results = await gapseq_fill_cmd(draft=draft, medium=medium, prefix=prefix, args=args)
    return GapseqModelResult(
        model=await one_ending(results, f"{prefix}-filled.RDS", "fill"),
        xml=await one_ending(results, f"{prefix}-filled.xml", "fill"),
        log=await files_ending(results, ".log"),
        results=results,
    )


@dataclass
class GapseqDoallResult:
    model: list[File]  # the .RDS models: draft and gap-filled
    xml: list[File]  # their SBML
    tbl: list[File]  # the find and find-transport tables
    fna: list[File]
    log: list[File]
    results: Dir


async def gapseq_doall(
    fasta: File, medium: File | None = None, db: File | None = None, args: str = ""
) -> GapseqDoallResult:
    """The whole gapseq reconstruction of ``fasta``, optionally on ``medium``."""
    results = await gapseq_doall_cmd(fasta=fasta, medium=medium, db=db, args=args)
    return GapseqDoallResult(
        model=await files_ending(results, ".RDS"),
        xml=await files_ending(results, ".xml"),
        tbl=await files_ending(results, ".tbl"),
        fna=await files_ending(results, ".fna"),
        log=await files_ending(results, ".log"),
        results=results,
    )
