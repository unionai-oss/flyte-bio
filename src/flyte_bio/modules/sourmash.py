"""sourmash — MinHash sketches of sequence files, and searches among them.

Exposes:

- :func:`sourmash_sketch`, which sketches one or more sequence files
  (FASTA/FASTQ, optionally gzipped) into ``<prefix>.sig``. Each file is staged
  as ``library/<name>``, as upstream, so signatures record that file name. With
  ``merge_sigs`` all files go into one signature named ``<prefix>``; otherwise
  there's one signature per file, all in ``<prefix>.sig``. ``args`` sets what
  to sketch; the default is upstream's ``dna --param-string
  scaled=1000,k=21,k=31,k=51,abund`` (its quotes dropped, as ``args`` isn't
  shell-parsed).
- :func:`sourmash_index`, which indexes signatures at one k-mer size into a
  Sequence Bloom Tree, ``<prefix>.sbt.zip``.
- :func:`sourmash_gather`, which finds the shortest list of a database's
  signatures (signature files or indexes) covering a query signature, written
  as ``<prefix>.csv.gz``; there's no CSV when nothing matches. Optionally it
  also saves the query's unmatched hashes (``<prefix>_unassigned.sig.zip``),
  the matching signatures (``<prefix>_matches.sig.zip``) and the prefetch
  matches (``<prefix>_prefetch.sig.zip``, ``<prefix>_prefetch.csv.gz``).
  ``args`` are extra options, such as ``--ksize 31`` or ``--threshold-bp``.

- :func:`sourmash_compare`, which compares signatures pairwise (similarity,
  or containment with ``--containment`` in ``args``) and writes the matrix as
  NumPy (``<prefix>_comp.npy`` with its labels,
  ``<prefix>_comp.npy.labels.txt``) and/or CSV (``<prefix>_comp.csv``). The
  signatures come as files (passed in name order, as upstream) and/or a
  ``file_list`` of signature paths (``--from-file``). It uses the task's CPUs
  (``--processes``).

As upstream, index, gather and compare inputs are staged in the working directory
under their own names, which gather's CSV records.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

SOURMASH_IMAGE = "quay.io/biocontainers/sourmash:4.9.4--hdfd78af_0"
DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory="4Gi")
DEFAULT_ARGS = "dna --param-string scaled=1000,k=21,k=31,k=51,abund"

sourmash_sketch_cmd = shell.create(
    name="sourmash_sketch",
    image=SOURMASH_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"library": list[File], "merge_sigs": bool, "prefix": str, "args": str},
    defaults={"merge_sigs": False, "args": DEFAULT_ARGS},
    outputs={"results": Dir},
    script=r"""
        LIBRARY=({inputs.library})
        P={inputs.prefix}
        ARGS={inputs.args}
        OUT={outputs.results}
        W=$(mktemp -d); cd "$W"
        mkdir library
        for F in "${LIBRARY[@]}"; do ln -s "$F" library/; done
        find -L library/ -type f > library.txt
        MERGE=()
        if [ {inputs.merge_sigs} = true ]; then MERGE=(--merge "$P"); fi
        sourmash sketch \
            $ARGS \
            "${MERGE[@]}" \
            --output "$P.sig" \
            --from-file library.txt
        mv "$P.sig" "$OUT"/
    """,
)

sourmash_index_cmd = shell.create(
    name="sourmash_index",
    image=SOURMASH_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"signatures": list[File], "ksize": int, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        SIGNATURES=({inputs.signatures})
        P={inputs.prefix}
        ARGS={inputs.args}
        OUT={outputs.results}
        W=$(mktemp -d); cd "$W"
        NAMES=()
        for F in "${SIGNATURES[@]}"; do ln -s "$F" .; NAMES+=("$(basename "$F")"); done
        sourmash index \
            --ksize {inputs.ksize} \
            $ARGS \
            "$P.sbt.zip" \
            "${NAMES[@]}" || exit 1
        mv "$P.sbt.zip" "$OUT"/
    """,
)

sourmash_gather_cmd = shell.create(
    name="sourmash_gather",
    image=SOURMASH_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "signature": File,
        "database": list[File],
        "prefix": str,
        "save_unassigned": bool,
        "save_matches_sig": bool,
        "save_prefetch": bool,
        "save_prefetch_csv": bool,
        "args": str,
    },
    defaults={
        "save_unassigned": False,
        "save_matches_sig": False,
        "save_prefetch": False,
        "save_prefetch_csv": False,
        "args": "",
    },
    outputs={"results": Dir},
    script=r"""
        SIGNATURE=({inputs.signature})
        DATABASE=({inputs.database})
        P={inputs.prefix}
        ARGS={inputs.args}
        OUT={outputs.results}
        W=$(mktemp -d); cd "$W"
        ln -s "${SIGNATURE[0]}" .
        QUERY=$(basename "${SIGNATURE[0]}")
        NAMES=()
        for F in "${DATABASE[@]}"; do ln -s "$F" .; NAMES+=("$(basename "$F")"); done
        SAVE=()
        if [ {inputs.save_unassigned} = true ]; then SAVE+=(--output-unassigned "${P}_unassigned.sig.zip"); fi
        if [ {inputs.save_matches_sig} = true ]; then SAVE+=(--save-matches "${P}_matches.sig.zip"); fi
        if [ {inputs.save_prefetch} = true ]; then SAVE+=(--save-prefetch "${P}_prefetch.sig.zip"); fi
        if [ {inputs.save_prefetch_csv} = true ]; then SAVE+=(--save-prefetch-csv "${P}_prefetch.csv.gz"); fi
        sourmash gather \
            $ARGS \
            --output "$P.csv.gz" \
            "${SAVE[@]}" \
            "$QUERY" \
            "${NAMES[@]}" || exit 1
        shopt -s nullglob
        mv ./*.csv.gz ./*.sig.zip "$OUT"/
    """,
)

sourmash_compare_cmd = shell.create(
    name="sourmash_compare",
    image=SOURMASH_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "signatures": list[File] | None,
        "file_list": File | None,
        "save_numpy_matrix": bool,
        "save_csv": bool,
        "prefix": str,
        "args": str,
    },
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        shopt -s nullglob; SIGNATURES=({{inputs.signatures}}); FILE_LIST=({{inputs.file_list}}); shopt -u nullglob
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        NAMES=()
        for F in "${{SIGNATURES[@]}}"; do ln -s "$F" .; NAMES+=("$(basename "$F")"); done
        mapfile -t NAMES < <(if [ ${{#NAMES[@]}} -gt 0 ]; then printf '%s\n' "${{NAMES[@]}}" | LC_ALL=C sort; fi)
        SAVE=()
        if [ {{inputs.save_numpy_matrix}} = true ]; then SAVE+=(--output "${{P}}_comp.npy"); fi
        if [ {{inputs.save_csv}} = true ]; then SAVE+=(--csv "${{P}}_comp.csv"); fi
        FROM_FILE=()
        if [ ${{#FILE_LIST[@]}} -gt 0 ]; then FROM_FILE=(--from-file "${{FILE_LIST[0]}}"); fi
        sourmash compare \
            $ARGS \
            --processes $TASK_CPUS \
            "${{SAVE[@]}}" \
            "${{FROM_FILE[@]}}" \
            "${{NAMES[@]}}" || exit 1
        shopt -s nullglob
        mv ./*_comp.npy ./*_comp.npy.labels.txt ./*_comp.csv "$OUT"/
    """,
)

env = flyte.TaskEnvironment.from_task(
    "sourmash",
    sourmash_sketch_cmd.as_task(),
    sourmash_index_cmd.as_task(),
    sourmash_gather_cmd.as_task(),
    sourmash_compare_cmd.as_task(),
)


async def sourmash_sketch(
    library: list[File], prefix: str = "sample", merge_sigs: bool = False, args: str = DEFAULT_ARGS
) -> File:
    """The sourmash signatures ``<prefix>.sig`` of ``library``'s files, merged into one with ``merge_sigs``."""
    results = await sourmash_sketch_cmd(library=library, merge_sigs=merge_sigs, prefix=prefix, args=args)
    signatures = await results.get_file(f"{prefix}.sig")
    if signatures is None:
        raise FileNotFoundError(f"sourmash wrote no {prefix}.sig")
    return signatures


async def sourmash_index(signatures: list[File], ksize: int, prefix: str = "sample", args: str = "") -> File:
    """A Sequence Bloom Tree ``<prefix>.sbt.zip`` of ``signatures`` at k-mer size ``ksize``."""
    results = await sourmash_index_cmd(signatures=signatures, ksize=ksize, prefix=prefix, args=args)
    index = await results.get_file(f"{prefix}.sbt.zip")
    if index is None:
        raise FileNotFoundError(f"sourmash wrote no {prefix}.sbt.zip")
    return index


@dataclass
class SourmashGatherResult:
    result: File | None  # <prefix>.csv.gz; None when nothing matched
    unassigned: File | None  # <prefix>_unassigned.sig.zip, with save_unassigned
    matches: File | None  # <prefix>_matches.sig.zip, with save_matches_sig
    prefetch: File | None  # <prefix>_prefetch.sig.zip, with save_prefetch
    prefetch_csv: File | None  # <prefix>_prefetch.csv.gz, with save_prefetch_csv
    results: Dir  # everything gather wrote


async def sourmash_gather(
    signature: File,
    database: list[File],
    prefix: str = "sample",
    save_unassigned: bool = False,
    save_matches_sig: bool = False,
    save_prefetch: bool = False,
    save_prefetch_csv: bool = False,
    args: str = "",
) -> SourmashGatherResult:
    """The shortest list of ``database``'s signatures that covers the query ``signature``."""
    results = await sourmash_gather_cmd(
        signature=signature,
        database=database,
        prefix=prefix,
        save_unassigned=save_unassigned,
        save_matches_sig=save_matches_sig,
        save_prefetch=save_prefetch,
        save_prefetch_csv=save_prefetch_csv,
        args=args,
    )
    return SourmashGatherResult(
        result=await results.get_file(f"{prefix}.csv.gz"),
        unassigned=await results.get_file(f"{prefix}_unassigned.sig.zip"),
        matches=await results.get_file(f"{prefix}_matches.sig.zip"),
        prefetch=await results.get_file(f"{prefix}_prefetch.sig.zip"),
        prefetch_csv=await results.get_file(f"{prefix}_prefetch.csv.gz"),
        results=results,
    )


@dataclass
class SourmashCompareResult:
    matrix: File | None  # <prefix>_comp.npy, with save_numpy_matrix
    labels: File | None  # <prefix>_comp.npy.labels.txt, with save_numpy_matrix
    csv: File | None  # <prefix>_comp.csv, with save_csv


async def sourmash_compare(
    signatures: list[File] | None = None,
    file_list: File | None = None,
    save_numpy_matrix: bool = True,
    save_csv: bool = False,
    prefix: str = "sample",
    args: str = "",
) -> SourmashCompareResult:
    """The pairwise comparison matrix of ``signatures`` and/or the signatures listed in ``file_list``."""
    if not save_numpy_matrix and not save_csv:
        raise ValueError("Supply either save_numpy_matrix, save_csv, or both or no output will be created")
    if not signatures and file_list is None:
        raise ValueError("Supply either signatures, file_list, or both")
    names = [f.name for f in signatures or []]
    if len(set(names)) != len(names):
        raise ValueError(f"sourmash_compare needs distinct signature names, got {names}")
    results = await sourmash_compare_cmd(
        signatures=signatures or None,
        file_list=file_list,
        save_numpy_matrix=save_numpy_matrix,
        save_csv=save_csv,
        prefix=prefix,
        args=args,
    )
    return SourmashCompareResult(
        matrix=await results.get_file(f"{prefix}_comp.npy"),
        labels=await results.get_file(f"{prefix}_comp.npy.labels.txt"),
        csv=await results.get_file(f"{prefix}_comp.csv"),
    )
