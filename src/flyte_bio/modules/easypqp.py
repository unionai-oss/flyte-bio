"""easypqp — spectral libraries from peptide identifications.

Exposes:

- :func:`easypqp_convert` — pair identifications (pepXML, idXML, …) with
  their spectra (mzML, mzXML, MGF): ``<spectra stem>.psmpkl`` (PSMs) and
  ``<spectra stem>.peakpkl`` (peaks), optionally with a custom Unimod XML.
- :func:`easypqp_library` — a spectral library, ``<prefix>.tsv``, from
  converted PSM and peak pickles, with per-run peaks
  ``<run>_run_peaks.tsv`` when written.

As upstream, the output names of ``convert`` follow the spectra file, since
``library`` pairs PSM and peak pickles by that name; matplotlib's and the
cache's directories are made writable in the working directory.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

EASYPQP_IMAGE = "quay.io/biocontainers/easypqp:0.1.60--pyhdfd78af_1"
# Upstream labels both tools process_single.
DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory="6Gi")

WRITABLE_CACHES = r"""
        W=$(mktemp -d); cd "$W"
        export HOME="$W" MPLCONFIGDIR="$W/.matplotlib" XDG_CACHE_HOME="$W/.cache"
        mkdir -p "$MPLCONFIGDIR" "$XDG_CACHE_HOME"
"""

easypqp_convert_cmd = shell.create(
    name="easypqp_convert",
    image=EASYPQP_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"id_file": File, "spectra": File, "unimod": File | None, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        ID=({inputs.id_file}); SPECTRA=({inputs.spectra})
        shopt -s nullglob; UNIMOD=({inputs.unimod}); shopt -u nullglob
        ARGS={inputs.args}
        OUT={outputs.results}
"""
    + WRITABLE_CACHES
    + r"""
        for f in "${ID[0]}" "${SPECTRA[0]}" "${UNIMOD[@]}"; do ln -s "$f" .; done
        EXTRA=()
        if [ ${#UNIMOD[@]} -gt 0 ]; then EXTRA=(--unimod "$(basename "${UNIMOD[0]}")"); fi
        easypqp convert --pepxml "$(basename "${ID[0]}")" --spectra "$(basename "${SPECTRA[0]}")" \
            "${EXTRA[@]}" $ARGS || exit 1
        mv ./*.psmpkl ./*.peakpkl "$OUT"/
    """,
)

easypqp_library_cmd = shell.create(
    name="easypqp_library",
    image=EASYPQP_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"psmpkls": list[File], "peakpkls": list[File], "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        PSM=({inputs.psmpkls}); PEAK=({inputs.peakpkls})
        P={inputs.prefix}
        ARGS={inputs.args}
        OUT={outputs.results}
"""
    + WRITABLE_CACHES
    + r"""
        NAMES=()
        for f in "${PSM[@]}" "${PEAK[@]}"; do ln -s "$f" .; NAMES+=("$(basename "$f")"); done
        easypqp library --out "$P.tsv" $ARGS "${NAMES[@]}" || exit 1
        shopt -s nullglob
        mv "$P.tsv" ./*_run_peaks.tsv "$OUT"/
    """,
)

env = flyte.TaskEnvironment.from_task("easypqp", easypqp_convert_cmd.as_task(), easypqp_library_cmd.as_task())


async def easypqp_pick(results: Dir, name: str) -> File:
    found = await results.get_file(name)
    if found is None:
        raise FileNotFoundError(f"easypqp wrote no {name}")
    return found


@dataclass
class EasypqpConvertResult:
    psmpkl: File  # <spectra stem>.psmpkl
    peakpkl: File  # <spectra stem>.peakpkl


async def easypqp_convert(
    id_file: File, spectra: File, unimod: File | None = None, args: str = ""
) -> EasypqpConvertResult:
    """PSM and peak pickles of ``id_file`` matched to ``spectra``, named after the spectra file."""
    results = await easypqp_convert_cmd(id_file=id_file, spectra=spectra, unimod=unimod, args=args)
    stem = spectra.path.rstrip("/").rsplit("/", 1)[-1].rsplit(".", 1)[0]
    return EasypqpConvertResult(
        psmpkl=await easypqp_pick(results, f"{stem}.psmpkl"),
        peakpkl=await easypqp_pick(results, f"{stem}.peakpkl"),
    )


@dataclass
class EasypqpLibraryResult:
    tsv: File  # <prefix>.tsv, the library
    run_peaks: list[File]  # <run>_run_peaks.tsv, when written


async def easypqp_library(
    psmpkls: list[File], peakpkls: list[File], prefix: str = "sample", args: str = ""
) -> EasypqpLibraryResult:
    """A spectral library, ``<prefix>.tsv``, from ``convert``'s PSM and peak pickles."""
    results = await easypqp_library_cmd(psmpkls=psmpkls, peakpkls=peakpkls, prefix=prefix, args=args)
    files = await results.list_files()
    run_peaks = sorted((f for f in files if f.path.endswith("_run_peaks.tsv")), key=lambda f: f.path)
    return EasypqpLibraryResult(tsv=await easypqp_pick(results, f"{prefix}.tsv"), run_peaks=run_peaks)
