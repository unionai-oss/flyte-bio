"""artic — the ARTIC network's tools for nanopore amplicon sequencing of viruses.

Exposes:

- :func:`artic_guppyplex` — gathers a directory of basecalled FASTQs into one
  length-filtered FASTQ, ``<prefix>.fastq.gz`` (``args`` carries guppyplex's
  options, e.g. ``--min-length 400 --max-length 700``).
- :func:`artic_minion` — the ARTIC consensus pipeline on one sample's reads:
  alignment to the reference, primer trimming against the scheme's BED,
  variant calling with Clair3 and a consensus sequence. It returns every
  ``<prefix>.*`` file artic writes and picks out the sorted, trimmed and
  primer-trimmed BAMs (with their ``.bai``), the consensus
  ``<prefix>.consensus.fasta`` and the passing variants
  ``<prefix>.pass.vcf.gz`` (with its ``.tbi``).

As upstream, minion's Clair3 model is ``model`` (e.g. ``r941_prom_hac_g360+g422``)
read from ``model_dir``; without a ``model_dir`` it is one of the models
bundled in the image (artic's ``bin`` directory's sibling ``models``).
``HDF5_PLUGIN_PATH`` is ``hdf5_plugin_path`` when given, else upstream's
default. minion's ``--threads`` follow the task's CPUs.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

ARTIC_IMAGE = "community.wave.seqera.io/library/artic:1.6.2--d4956cdc155b8612"
# Upstream: process_high.
GUPPYPLEX_RESOURCES = flyte.Resources(cpu=2, memory="4Gi")
MINION_RESOURCES = flyte.Resources(cpu=4, memory="16Gi")

artic_guppyplex_cmd = shell.create(
    name="artic_guppyplex",
    image=ARTIC_IMAGE,
    resources=GUPPYPLEX_RESOURCES,
    inputs={"fastq_dir": Dir, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        FASTQ_DIR={{inputs.fastq_dir}}
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        artic \
            guppyplex \
            $ARGS \
            --directory "$FASTQ_DIR" \
            --output "$P.fastq" || exit 1
        pigz -p $TASK_CPUS "$P.fastq" || exit 1
        mv "$P.fastq.gz" "$OUT"/
    """,
)

artic_minion_cmd = shell.create(
    name="artic_minion",
    image=ARTIC_IMAGE,
    resources=MINION_RESOURCES,
    inputs={
        "fastq": File,
        "model_dir": Dir | None,
        "model": str,
        "fasta": File,
        "bed": File,
        "hdf5_plugin_path": Dir | None,
        "prefix": str,
        "args": str,
    },
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        FQ=({{inputs.fastq}}); FA=({{inputs.fasta}}); BED=({{inputs.bed}})
        MODEL={{inputs.model}}
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        # An unset optional Dir is staged as a file, not a directory.
        MODEL_DIR="$(which artic | sed 's/artic/models/')"
        GIVEN={{inputs.model_dir}}
        if [ -d "$GIVEN" ]; then MODEL_DIR="$GIVEN"; fi
        HDF5_PLUGIN_PATH=/usr/local/lib/python3.6/site-packages/ont_fast5_api/vbz_plugin
        GIVEN={{inputs.hdf5_plugin_path}}
        if [ -d "$GIVEN" ]; then HDF5_PLUGIN_PATH="$GIVEN"; fi
        export HDF5_PLUGIN_PATH
        # Inputs are staged by name in their own directory, so <prefix>.* is only artic's outputs.
        W=$(mktemp -d); cd "$W"; mkdir inputs
        for F in "${{FQ[0]}}" "${{FA[0]}}" "${{BED[0]}}"; do ln -s "$F" inputs/"$(basename "$F")"; done
        artic \
            minion \
            $ARGS \
            --threads $TASK_CPUS \
            --read-file inputs/"$(basename "${{FQ[0]}}")" \
            --bed inputs/"$(basename "${{BED[0]}}")" \
            --ref inputs/"$(basename "${{FA[0]}}")" \
            --model-dir "$MODEL_DIR" \
            --model "$MODEL" \
            "$P" || exit 1
        mv "$P".* "$OUT"/
        shopt -s nullglob
        for F in ./*.json; do mv "$F" "$OUT"/; done
    """,
)

env = flyte.TaskEnvironment.from_task("artic", artic_guppyplex_cmd.as_task(), artic_minion_cmd.as_task())


async def pick_file(results: Dir, name: str, tool: str) -> File:
    found = await results.get_file(name)
    if found is None:
        raise FileNotFoundError(f"{tool} wrote no {name}")
    return found


async def artic_guppyplex(fastq_dir: Dir, prefix: str = "sample", args: str = "") -> File:
    """The reads in ``fastq_dir``'s FASTQs, filtered and gathered into ``<prefix>.fastq.gz``."""
    results = await artic_guppyplex_cmd(fastq_dir=fastq_dir, prefix=prefix, args=args)
    return await pick_file(results, f"{prefix}.fastq.gz", "artic guppyplex")


@dataclass
class ArticMinionResult:
    bam: File  # <prefix>.sorted.bam
    bai: File
    bam_trimmed: File  # <prefix>.trimmed.rg.sorted.bam
    bai_trimmed: File
    bam_primertrimmed: File  # <prefix>.primertrimmed.rg.sorted.bam
    bai_primertrimmed: File
    fasta: File  # <prefix>.consensus.fasta
    vcf: File  # <prefix>.pass.vcf.gz
    tbi: File
    json: list[File]  # any *.json artic wrote
    results: Dir  # every <prefix>.* file, and the JSONs


async def artic_minion(
    fastq: File,
    fasta: File,
    bed: File,
    model: str,
    model_dir: Dir | None = None,
    hdf5_plugin_path: Dir | None = None,
    prefix: str = "sample",
    args: str = "",
) -> ArticMinionResult:
    """The ARTIC pipeline on ``fastq`` against the scheme's reference ``fasta`` and primer ``bed``."""
    results = await artic_minion_cmd(
        fastq=fastq,
        model_dir=model_dir,
        model=model,
        fasta=fasta,
        bed=bed,
        hdf5_plugin_path=hdf5_plugin_path,
        prefix=prefix,
        args=args,
    )

    async def pick(suffix: str) -> File:
        return await pick_file(results, f"{prefix}.{suffix}", "artic minion")

    files = await results.list_files()
    return ArticMinionResult(
        bam=await pick("sorted.bam"),
        bai=await pick("sorted.bam.bai"),
        bam_trimmed=await pick("trimmed.rg.sorted.bam"),
        bai_trimmed=await pick("trimmed.rg.sorted.bam.bai"),
        bam_primertrimmed=await pick("primertrimmed.rg.sorted.bam"),
        bai_primertrimmed=await pick("primertrimmed.rg.sorted.bam.bai"),
        fasta=await pick("consensus.fasta"),
        vcf=await pick("pass.vcf.gz"),
        tbi=await pick("pass.vcf.gz.tbi"),
        json=sorted((f for f in files if f.name.endswith(".json")), key=lambda f: f.name),
        results=results,
    )
