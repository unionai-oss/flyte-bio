"""transdecoder — coding regions in transcripts with TransDecoder.

Exposes:

- :func:`transdecoder_longorf`, ``TransDecoder.LongOrfs`` on a transcript
  FASTA (may be gzipped): the long open reading frames, in its
  ``<fasta name without .gz>.transdecoder_dir`` (``longest_orfs.pep``,
  ``longest_orfs.gff3``, ``longest_orfs.cds`` and ``base_freqs.dat``, which
  the result also picks out).
- :func:`transdecoder_predict`, ``TransDecoder.Predict`` on the same FASTA
  and that directory: the likely coding regions, as
  ``<fasta name without .gz>.transdecoder.{pep,gff3,cds,bed}``.

As upstream, Predict works on a private directory of links to LongOrfs'
files, so it never changes LongOrfs' output. ``args`` passes each tool's own
options (``-m`` minimum protein length, ``-S`` strand-specific,
``--single_best_only`` and so on).
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

TRANSDECODER_IMAGE = "quay.io/biocontainers/transdecoder:5.7.1--pl5321hdfd78af_0"
# Upstream's process_medium; both tools are single-threaded unless args asks
# for more (Predict's --cpu with homology searches).
DEFAULT_RESOURCES = flyte.Resources(cpu=2, memory="8Gi")

transdecoder_longorf_cmd = shell.create(
    name="transdecoder_longorf",
    image=TRANSDECODER_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"fasta": File, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        FA=({inputs.fasta})
        ARGS={inputs.args}
        OUT={outputs.results}
        W=$(mktemp -d); cd "$W"
        NAME=$(basename "${FA[0]}")
        ln -s "${FA[0]}" "$NAME"
        TransDecoder.LongOrfs \
            $ARGS \
            -O longorfs \
            -t "$NAME" || exit 1
        mv "longorfs/${NAME%.gz}.transdecoder_dir" "$OUT/"
    """,
)

transdecoder_predict_cmd = shell.create(
    name="transdecoder_predict",
    image=TRANSDECODER_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"fasta": File, "folder": Dir, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        FA=({inputs.fasta})
        FOLDER={inputs.folder}
        ARGS={inputs.args}
        OUT={outputs.results}
        W=$(mktemp -d); cd "$W"
        NAME=$(basename "${FA[0]}")
        ln -s "${FA[0]}" "$NAME"
        if [ -z "$(ls -A "$FOLDER")" ]; then echo "TransDecoder.LongOrfs' directory is empty: $FOLDER" >&2; exit 1; fi
        # Predict reads and writes <name>.transdecoder_dir beside its output.
        mkdir "${NAME%.gz}.transdecoder_dir"
        find "$FOLDER" -mindepth 1 -maxdepth 1 -exec ln -s {} "${NAME%.gz}.transdecoder_dir"/ \;
        TransDecoder.Predict \
            $ARGS \
            -O . \
            -t "$NAME" || exit 1
        mv "${NAME%.gz}".transdecoder.pep "${NAME%.gz}".transdecoder.gff3 \
            "${NAME%.gz}".transdecoder.cds "${NAME%.gz}".transdecoder.bed "$OUT/"
    """,
)

env = flyte.TaskEnvironment.from_task(
    "transdecoder", transdecoder_longorf_cmd.as_task(), transdecoder_predict_cmd.as_task()
)


def stem(fasta: File) -> str:
    name = fasta.name
    return name.removesuffix(".gz")


async def picked(results: Dir, name: str, tool: str) -> File:
    found = await results.get_file(name)
    if found is None:
        raise FileNotFoundError(f"TransDecoder.{tool} wrote no {name}")
    return found


@dataclass
class TransdecoderLongorfResult:
    pep: File  # longest_orfs.pep
    gff3: File  # longest_orfs.gff3
    cds: File  # longest_orfs.cds
    dat: File  # base_freqs.dat
    folder: Dir  # <fasta name>.transdecoder_dir, for transdecoder_predict


async def transdecoder_longorf(fasta: File, args: str = "") -> TransdecoderLongorfResult:
    """``TransDecoder.LongOrfs`` on ``fasta``: its ``<fasta name>.transdecoder_dir`` and the files in it."""
    results = await transdecoder_longorf_cmd(fasta=fasta, args=args)
    folder_name = f"{stem(fasta)}.transdecoder_dir"
    folder = Dir.from_existing_remote(f"{results.path.rstrip('/')}/{folder_name}")
    return TransdecoderLongorfResult(
        pep=await picked(folder, "longest_orfs.pep", "LongOrfs"),
        gff3=await picked(folder, "longest_orfs.gff3", "LongOrfs"),
        cds=await picked(folder, "longest_orfs.cds", "LongOrfs"),
        dat=await picked(folder, "base_freqs.dat", "LongOrfs"),
        folder=folder,
    )


@dataclass
class TransdecoderPredictResult:
    pep: File  # <fasta name>.transdecoder.pep
    gff3: File  # <fasta name>.transdecoder.gff3
    cds: File  # <fasta name>.transdecoder.cds
    bed: File  # <fasta name>.transdecoder.bed


async def transdecoder_predict(fasta: File, folder: Dir, args: str = "") -> TransdecoderPredictResult:
    """``TransDecoder.Predict`` on ``fasta`` with ``folder`` (from :func:`transdecoder_longorf`)."""
    results = await transdecoder_predict_cmd(fasta=fasta, folder=folder, args=args)
    base = f"{stem(fasta)}.transdecoder"
    return TransdecoderPredictResult(
        pep=await picked(results, f"{base}.pep", "Predict"),
        gff3=await picked(results, f"{base}.gff3", "Predict"),
        cds=await picked(results, f"{base}.cds", "Predict"),
        bed=await picked(results, f"{base}.bed", "Predict"),
    )
