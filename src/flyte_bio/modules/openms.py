"""openms — OpenMS tools for mass-spectrometry proteomics and immunopeptidomics.

Exposes:

- :func:`openms_decoydatabase` — a FASTA with decoy proteins added
  (``DecoyDatabase``): ``<prefix>.fasta``, by default ``<fasta stem>_decoy.fasta``.
- :func:`openms_featurefinderidentification` — quantify the features of an
  mzML that peptide identifications point at (``FeatureFinderIdentification``):
  ``<prefix>.<out_type>`` (``featureXML``).
- :func:`openms_fileconverter` — convert between MS file formats
  (``FileConverter``): ``<prefix>.<out_type>``.
- :func:`openms_filefilter` — filter an mzML, featureXML or consensusXML
  (``FileFilter``): ``<prefix>.<input extension>``.
- :func:`openms_idconflictresolver` — keep one identification per feature
  (``IDConflictResolver``): ``<prefix>.<input extension>``.
- :func:`openms_idfilter` — filter peptide and protein identifications
  (``IDFilter``), optionally against a filter file: ``<prefix>.<input extension>``.
- :func:`openms_idmerger` — merge idXML files (``IDMerger``): ``<prefix>.idXML``.
- :func:`openms_idripper` — split a merged idXML back into one idXML per
  originating file (``IDRipper``).
- :func:`openms_idscoreswitcher` — switch the main score of identifications
  (``IDScoreSwitcher``): ``<prefix>.idXML``.
- :func:`openms_mapaligneridentification` — retention-time transformations
  aligning runs by their identifications (``MapAlignerIdentification``): one
  ``<input stem>.trafoXML`` per input.
- :func:`openms_maprttransformer` — apply a retention-time transformation
  (``MapRTTransformer``): ``<prefix>.<input extension>``.
- :func:`openms_mztabexporter` — export identifications or features as mzTab
  (``MzTabExporter``): ``<prefix>.mzTab``.
- :func:`openms_peakpickerhires` — centroid profile spectra
  (``PeakPickerHiRes``): ``<prefix>.mzML``.
- :func:`openms_peptideindexer` — map peptides to their proteins in a FASTA,
  marking targets and decoys (``PeptideIndexer``): ``<prefix>.idXML``.
- :func:`openms_psmfeatureextractor` — add search-engine features for
  rescoring (``PSMFeatureExtractor``): ``<prefix>.idXML``.
- :func:`openms_textexporter` — export identifications or features as text
  (``TextExporter``): ``<prefix>.tsv``.

As upstream, each tool runs with the task's CPUs as ``-threads``, on inputs
linked into its working directory under their own names (OpenMS records the
names it is given in its output). ``args`` is parsed with shell quoting, so a
value with spaces can be quoted as upstream writes it
(``-enzyme 'unspecific cleavage'``). Upstream's default prefixes come from the
sample name (``meta.id``, plus ``_resolved``, ``_aligned`` or ``_indexed`` for
IDConflictResolver, MapRTTransformer and PeptideIndexer); here ``prefix`` is
the whole output name before its extension. A prefix that would overwrite the
input is refused, as upstream.
"""

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

OPENMS_IMAGE = "community.wave.seqera.io/library/openms:3.6.0--88df48c0a426b7f7"
# Upstream labels most tools process_single, PeakPickerHiRes process_low and
# FeatureFinderIdentification process_medium.
SINGLE_RESOURCES = flyte.Resources(cpu=1, memory="4Gi")
LOW_RESOURCES = flyte.Resources(cpu=2, memory="4Gi")
MEDIUM_RESOURCES = flyte.Resources(cpu=4, memory="8Gi")

# Run in a scratch directory with a writable HOME (OpenMS keeps settings
# there); parse ``args`` with shell quoting into the array ARGS. ``stage``
# links files into the working directory under their own names.
OPENMS_PREAMBLE = (
    TASK_LIMITS
    + r"""
        OUT={outputs.results}
        W=$(mktemp -d); cd "$W"
        export HOME="$W"
        ARGS_TEXT={inputs.args}
        set -f; eval "ARGS=($ARGS_TEXT)"; set +f
        stage() { local f; for f in "$@"; do ln -s "$f" "$W/$(basename "$f")"; done; }
"""
)


def openms_cmd(name: str, inputs: dict, defaults: dict, body: str, resources: flyte.Resources = SINGLE_RESOURCES):
    """A shell task running ``body`` (after :data:`OPENMS_PREAMBLE`) in the OpenMS image, writing into ``results``."""
    return shell.create(
        name=name,
        image=OPENMS_IMAGE,
        resources=resources,
        inputs={**inputs, "args": str},
        defaults={**defaults, "args": ""},
        outputs={"results": Dir},
        script=OPENMS_PREAMBLE + body,
    )


openms_decoydatabase_cmd = openms_cmd(
    "openms_decoydatabase",
    {"fasta": File, "prefix": str},
    {},
    r"""
        FA=({inputs.fasta}); stage "${FA[0]}"
        P={inputs.prefix}
        DecoyDatabase -in "$(basename "${FA[0]}")" -out "$P.fasta" -threads $TASK_CPUS "${ARGS[@]}" || exit 1
        mv "$P.fasta" "$OUT"/
    """,
)

openms_featurefinderidentification_cmd = openms_cmd(
    "openms_featurefinderidentification",
    {"mzml": File, "id_int": File, "id_ext": File | None, "out_type": str, "prefix": str},
    {"out_type": "featureXML"},
    r"""
        MZML=({inputs.mzml}); ID_INT=({inputs.id_int})
        shopt -s nullglob; ID_EXT=({inputs.id_ext}); shopt -u nullglob
        stage "${MZML[0]}" "${ID_INT[0]}" "${ID_EXT[@]}"
        EXT=()
        if [ ${#ID_EXT[@]} -gt 0 ]; then EXT=(-id_ext "$(basename "${ID_EXT[0]}")"); fi
        OUT_NAME={inputs.prefix}.{inputs.out_type}
        FeatureFinderIdentification -in "$(basename "${MZML[0]}")" -id "$(basename "${ID_INT[0]}")" "${EXT[@]}" \
            -out "$OUT_NAME" -threads $TASK_CPUS "${ARGS[@]}" || exit 1
        mv "$OUT_NAME" "$OUT"/
    """,
    MEDIUM_RESOURCES,
)

openms_fileconverter_cmd = openms_cmd(
    "openms_fileconverter",
    {"input_file": File, "out_type": str, "prefix": str},
    {},
    r"""
        IN=({inputs.input_file}); stage "${IN[0]}"
        OUT_NAME={inputs.prefix}.{inputs.out_type}
        FileConverter -in "$(basename "${IN[0]}")" -out "$OUT_NAME" -threads $TASK_CPUS "${ARGS[@]}" || exit 1
        mv "$OUT_NAME" "$OUT"/
    """,
)

# Tools that read one file and write <prefix>.<the input's extension>.
SAME_TYPE_BODY = r"""
        IN=({inputs.input_file}); stage "${IN[0]}"
        NAME=$(basename "${IN[0]}")
        OUT_NAME={inputs.prefix}."${NAME##*.}"
        TOOL -in "$NAME" -out "$OUT_NAME" -threads $TASK_CPUS "${ARGS[@]}" || exit 1
        mv "$OUT_NAME" "$OUT"/
"""

openms_filefilter_cmd = openms_cmd(
    "openms_filefilter", {"input_file": File, "prefix": str}, {}, SAME_TYPE_BODY.replace("TOOL", "FileFilter")
)

openms_idconflictresolver_cmd = openms_cmd(
    "openms_idconflictresolver",
    {"input_file": File, "prefix": str},
    {},
    SAME_TYPE_BODY.replace("TOOL", "IDConflictResolver"),
)

openms_idfilter_cmd = openms_cmd(
    "openms_idfilter",
    {"input_file": File, "filter_file": File | None, "filter_criterion": str, "prefix": str},
    {"filter_criterion": "-whitelist:peptides"},
    r"""
        IN=({inputs.input_file})
        shopt -s nullglob; FILTER_FILE=({inputs.filter_file}); shopt -u nullglob
        stage "${IN[0]}" "${FILTER_FILE[@]}"
        NAME=$(basename "${IN[0]}")
        OUT_NAME={inputs.prefix}."${NAME##*.}"
        FILTER=()
        if [ ${#FILTER_FILE[@]} -gt 0 ]; then FILTER=({inputs.filter_criterion} "$(basename "${FILTER_FILE[0]}")"); fi
        IDFilter -in "$NAME" -out "$OUT_NAME" -threads $TASK_CPUS "${FILTER[@]}" "${ARGS[@]}" || exit 1
        mv "$OUT_NAME" "$OUT"/
    """,
)

openms_idmerger_cmd = openms_cmd(
    "openms_idmerger",
    {"idxmls": list[File], "prefix": str},
    {},
    r"""
        IN=({inputs.idxmls}); stage "${IN[@]}"
        NAMES=(); for f in "${IN[@]}"; do NAMES+=("$(basename "$f")"); done
        P={inputs.prefix}
        IDMerger -in "${NAMES[@]}" -out "$P.idXML" -threads $TASK_CPUS "${ARGS[@]}" || exit 1
        mv "$P.idXML" "$OUT"/
    """,
)

openms_idripper_cmd = openms_cmd(
    "openms_idripper",
    {"merged_idxml": File},
    {},
    r"""
        IN=({inputs.merged_idxml}); stage "${IN[0]}"
        NAME=$(basename "${IN[0]}")
        IDRipper -in "$NAME" -out . -threads $TASK_CPUS "${ARGS[@]}" || exit 1
        for f in ./*.idXML; do
            if [ "$(basename "$f")" != "$NAME" ] || [ ! -L "$f" ]; then mv "$f" "$OUT"/; fi
        done
    """,
)

openms_idscoreswitcher_cmd = openms_cmd(
    "openms_idscoreswitcher",
    {"idxml": File, "prefix": str},
    {},
    r"""
        IN=({inputs.idxml}); stage "${IN[0]}"
        P={inputs.prefix}
        IDScoreSwitcher -in "$(basename "${IN[0]}")" -out "$P.idXML" -threads $TASK_CPUS "${ARGS[@]}" || exit 1
        mv "$P.idXML" "$OUT"/
    """,
)

openms_mapaligneridentification_cmd = openms_cmd(
    "openms_mapaligneridentification",
    {"id_files": list[File]},
    {},
    r"""
        IN=({inputs.id_files}); stage "${IN[@]}"
        NAMES=(); TRAFOS=()
        for f in "${IN[@]}"; do
            NAME=$(basename "$f"); NAMES+=("$NAME"); TRAFOS+=("${NAME%.*}.trafoXML")
        done
        MapAlignerIdentification -in "${NAMES[@]}" -trafo_out "${TRAFOS[@]}" -threads $TASK_CPUS "${ARGS[@]}" || exit 1
        mv "${TRAFOS[@]}" "$OUT"/
    """,
)

openms_maprttransformer_cmd = openms_cmd(
    "openms_maprttransformer",
    {"input_file": File, "trafoxml": File, "prefix": str},
    {},
    r"""
        IN=({inputs.input_file}); TRAFO=({inputs.trafoxml}); stage "${IN[0]}" "${TRAFO[0]}"
        NAME=$(basename "${IN[0]}")
        OUT_NAME={inputs.prefix}."${NAME##*.}"
        MapRTTransformer -in "$NAME" -trafo_in "$(basename "${TRAFO[0]}")" -out "$OUT_NAME" \
            -threads $TASK_CPUS "${ARGS[@]}" || exit 1
        mv "$OUT_NAME" "$OUT"/
    """,
)

openms_mztabexporter_cmd = openms_cmd(
    "openms_mztabexporter",
    {"input_file": File, "prefix": str},
    {},
    r"""
        IN=({inputs.input_file}); stage "${IN[0]}"
        P={inputs.prefix}
        MzTabExporter -in "$(basename "${IN[0]}")" -out "$P.mzTab" -threads $TASK_CPUS "${ARGS[@]}" || exit 1
        mv "$P.mzTab" "$OUT"/
    """,
)

openms_peakpickerhires_cmd = openms_cmd(
    "openms_peakpickerhires",
    {"mzml": File, "prefix": str},
    {},
    r"""
        IN=({inputs.mzml}); stage "${IN[0]}"
        P={inputs.prefix}
        PeakPickerHiRes -in "$(basename "${IN[0]}")" -out "$P.mzML" -threads $TASK_CPUS "${ARGS[@]}" || exit 1
        mv "$P.mzML" "$OUT"/
    """,
    LOW_RESOURCES,
)

openms_peptideindexer_cmd = openms_cmd(
    "openms_peptideindexer",
    {"idxml": File, "fasta": File | None, "prefix": str},
    {},
    r"""
        IN=({inputs.idxml})
        shopt -s nullglob; FA=({inputs.fasta}); shopt -u nullglob
        stage "${IN[0]}" "${FA[@]}"
        FASTA=()
        if [ ${#FA[@]} -gt 0 ]; then FASTA=(-fasta "$(basename "${FA[0]}")"); fi
        P={inputs.prefix}
        PeptideIndexer -in "$(basename "${IN[0]}")" "${FASTA[@]}" -out "$P.idXML" \
            -threads $TASK_CPUS "${ARGS[@]}" || exit 1
        mv "$P.idXML" "$OUT"/
    """,
)

openms_psmfeatureextractor_cmd = openms_cmd(
    "openms_psmfeatureextractor",
    {"id_file": File, "prefix": str},
    {},
    r"""
        IN=({inputs.id_file}); stage "${IN[0]}"
        P={inputs.prefix}
        PSMFeatureExtractor -in "$(basename "${IN[0]}")" -out "$P.idXML" -threads $TASK_CPUS "${ARGS[@]}" || exit 1
        mv "$P.idXML" "$OUT"/
    """,
)

openms_textexporter_cmd = openms_cmd(
    "openms_textexporter",
    {"input_file": File, "prefix": str},
    {},
    r"""
        IN=({inputs.input_file}); stage "${IN[0]}"
        P={inputs.prefix}
        TextExporter -in "$(basename "${IN[0]}")" -out "$P.tsv" -threads $TASK_CPUS "${ARGS[@]}" || exit 1
        mv "$P.tsv" "$OUT"/
    """,
)

env = flyte.TaskEnvironment.from_task(
    "openms",
    openms_decoydatabase_cmd.as_task(),
    openms_featurefinderidentification_cmd.as_task(),
    openms_fileconverter_cmd.as_task(),
    openms_filefilter_cmd.as_task(),
    openms_idconflictresolver_cmd.as_task(),
    openms_idfilter_cmd.as_task(),
    openms_idmerger_cmd.as_task(),
    openms_idripper_cmd.as_task(),
    openms_idscoreswitcher_cmd.as_task(),
    openms_mapaligneridentification_cmd.as_task(),
    openms_maprttransformer_cmd.as_task(),
    openms_mztabexporter_cmd.as_task(),
    openms_peakpickerhires_cmd.as_task(),
    openms_peptideindexer_cmd.as_task(),
    openms_psmfeatureextractor_cmd.as_task(),
    openms_textexporter_cmd.as_task(),
)


def file_name(file: File) -> str:
    return file.path.rstrip("/").rsplit("/", 1)[-1]


def file_stem(file: File) -> str:
    return file_name(file).rsplit(".", 1)[0]


def file_extension(file: File) -> str:
    return file_name(file).rsplit(".", 1)[-1]


def check_prefix(file: File, out_name: str) -> None:
    if file_name(file) == out_name:
        raise ValueError(f"input and output are both named {out_name}: choose another prefix")


async def openms_pick(results: Dir, name: str) -> File:
    """The file ``name`` from an OpenMS tool's ``results``."""
    found = await results.get_file(name)
    if found is None:
        raise FileNotFoundError(f"OpenMS wrote no {name}")
    return found


async def openms_decoydatabase(fasta: File, prefix: str | None = None, args: str = "") -> File:
    """``fasta`` with decoys added: ``<prefix>.fasta``, by default ``<fasta stem>_decoy.fasta``."""
    prefix = prefix or f"{file_stem(fasta)}_decoy"
    results = await openms_decoydatabase_cmd(fasta=fasta, prefix=prefix, args=args)
    return await openms_pick(results, f"{prefix}.fasta")


async def openms_featurefinderidentification(
    mzml: File,
    id_int: File,
    id_ext: File | None = None,
    out_type: str = "featureXML",
    prefix: str = "sample",
    args: str = "",
) -> File:
    """Features of ``mzml`` at the identifications ``id_int`` (and external ``id_ext``): ``<prefix>.<out_type>``."""
    results = await openms_featurefinderidentification_cmd(
        mzml=mzml, id_int=id_int, id_ext=id_ext, out_type=out_type, prefix=prefix, args=args
    )
    return await openms_pick(results, f"{prefix}.{out_type}")


async def openms_fileconverter(input_file: File, out_type: str, prefix: str = "sample", args: str = "") -> File:
    """``input_file`` converted to ``out_type``: ``<prefix>.<out_type>``."""
    results = await openms_fileconverter_cmd(input_file=input_file, out_type=out_type, prefix=prefix, args=args)
    return await openms_pick(results, f"{prefix}.{out_type}")


async def openms_filefilter(input_file: File, prefix: str = "sample", args: str = "") -> File:
    """``input_file`` (mzML, featureXML or consensusXML) filtered: ``<prefix>.<its extension>``."""
    out_name = f"{prefix}.{file_extension(input_file)}"
    check_prefix(input_file, out_name)
    results = await openms_filefilter_cmd(input_file=input_file, prefix=prefix, args=args)
    return await openms_pick(results, out_name)


async def openms_idconflictresolver(input_file: File, prefix: str = "sample_resolved", args: str = "") -> File:
    """``input_file`` (featureXML or consensusXML) with one identification per feature."""
    out_name = f"{prefix}.{file_extension(input_file)}"
    check_prefix(input_file, out_name)
    results = await openms_idconflictresolver_cmd(input_file=input_file, prefix=prefix, args=args)
    return await openms_pick(results, out_name)


async def openms_idfilter(
    input_file: File,
    filter_file: File | None = None,
    filter_criterion: str = "-whitelist:peptides",
    prefix: str = "sample",
    args: str = "",
) -> File:
    """``input_file`` (idXML or consensusXML) filtered, optionally by ``filter_criterion filter_file``."""
    out_name = f"{prefix}.{file_extension(input_file)}"
    check_prefix(input_file, out_name)
    results = await openms_idfilter_cmd(
        input_file=input_file, filter_file=filter_file, filter_criterion=filter_criterion, prefix=prefix, args=args
    )
    return await openms_pick(results, out_name)


async def openms_idmerger(idxmls: list[File], prefix: str = "sample", args: str = "") -> File:
    """``idxmls`` merged into ``<prefix>.idXML``."""
    results = await openms_idmerger_cmd(idxmls=idxmls, prefix=prefix, args=args)
    return await openms_pick(results, f"{prefix}.idXML")


async def openms_idripper(merged_idxml: File, args: str = "") -> list[File]:
    """``merged_idxml`` split into one idXML per originating file, sorted by name."""
    results = await openms_idripper_cmd(merged_idxml=merged_idxml, args=args)
    return sorted(await results.list_files(), key=file_name)


async def openms_idscoreswitcher(idxml: File, prefix: str = "sample", args: str = "") -> File:
    """``idxml`` with its main score switched (as ``args`` say): ``<prefix>.idXML``."""
    check_prefix(idxml, f"{prefix}.idXML")
    results = await openms_idscoreswitcher_cmd(idxml=idxml, prefix=prefix, args=args)
    return await openms_pick(results, f"{prefix}.idXML")


async def openms_mapaligneridentification(id_files: list[File], args: str = "") -> list[File]:
    """Retention-time transformations aligning ``id_files``: ``<stem>.trafoXML`` for each, in order."""
    results = await openms_mapaligneridentification_cmd(id_files=id_files, args=args)
    return [await openms_pick(results, f"{file_stem(f)}.trafoXML") for f in id_files]


async def openms_maprttransformer(
    input_file: File, trafoxml: File, prefix: str = "sample_aligned", args: str = ""
) -> File:
    """``input_file`` with the retention-time transformation ``trafoxml`` applied: ``<prefix>.<its extension>``."""
    out_name = f"{prefix}.{file_extension(input_file)}"
    check_prefix(input_file, out_name)
    results = await openms_maprttransformer_cmd(input_file=input_file, trafoxml=trafoxml, prefix=prefix, args=args)
    return await openms_pick(results, out_name)


async def openms_mztabexporter(input_file: File, prefix: str = "sample", args: str = "") -> File:
    """``input_file`` exported as ``<prefix>.mzTab``."""
    results = await openms_mztabexporter_cmd(input_file=input_file, prefix=prefix, args=args)
    return await openms_pick(results, f"{prefix}.mzTab")


async def openms_peakpickerhires(mzml: File, prefix: str = "sample", args: str = "") -> File:
    """``mzml``'s spectra centroided: ``<prefix>.mzML``."""
    results = await openms_peakpickerhires_cmd(mzml=mzml, prefix=prefix, args=args)
    return await openms_pick(results, f"{prefix}.mzML")


async def openms_peptideindexer(
    idxml: File, fasta: File | None = None, prefix: str = "sample_indexed", args: str = ""
) -> File:
    """``idxml``'s peptides mapped to their proteins in ``fasta``: ``<prefix>.idXML``."""
    results = await openms_peptideindexer_cmd(idxml=idxml, fasta=fasta, prefix=prefix, args=args)
    return await openms_pick(results, f"{prefix}.idXML")


async def openms_psmfeatureextractor(id_file: File, prefix: str = "sample", args: str = "") -> File:
    """``id_file`` with rescoring features added: ``<prefix>.idXML``."""
    check_prefix(id_file, f"{prefix}.idXML")
    results = await openms_psmfeatureextractor_cmd(id_file=id_file, prefix=prefix, args=args)
    return await openms_pick(results, f"{prefix}.idXML")


async def openms_textexporter(input_file: File, prefix: str = "sample", args: str = "") -> File:
    """``input_file`` exported as text: ``<prefix>.tsv``."""
    results = await openms_textexporter_cmd(input_file=input_file, prefix=prefix, args=args)
    return await openms_pick(results, f"{prefix}.tsv")
