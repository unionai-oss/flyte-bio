"""kallisto — pseudoalignment index and transcript quantification.

Exposes:

- :func:`kallisto_index`, which builds a kallisto index from a transcript
  FASTA (optionally gzipped). As upstream, the index is a single file named
  ``kallisto``.
- :func:`kallisto_quant`, which quantifies one sample's reads against that
  index. Single-end reads need the fragment length and its standard deviation
  (``--single --fragment-length --sd``), as upstream. Strandedness
  (``--fr-stranded``/``--rf-stranded``) rides in ``args``. Outputs keep
  upstream's names: the results directory ``<prefix>`` (``abundance.h5``,
  ``abundance.tsv``, ``kallisto_quant.log``, ``run_info.json``), and copies of
  its log and run info as ``<prefix>.log`` and ``<prefix>.run_info.json``.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

KALLISTO_IMAGE = "community.wave.seqera.io/library/kallisto:0.52.0--31c771060d82d25c"
# Upstream labels it process_medium; a human transcriptome needs about 16 GiB.
INDEX_RESOURCES = flyte.Resources(cpu=1, memory="16Gi")
# Upstream labels quant process_high; it needs the index in memory and little more.
QUANT_RESOURCES = flyte.Resources(cpu=4, memory="16Gi")

# Upstream passes no thread count, so kallisto indexes with its default (1);
# `-t` in `args` speeds it up.
kallisto_index_cmd = shell.create(
    name="kallisto_index",
    image=KALLISTO_IMAGE,
    resources=INDEX_RESOURCES,
    inputs={"fasta": File, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        FA=({inputs.fasta})
        ARGS={inputs.args}
        OUT={outputs.results}
        W=$(mktemp -d); cd "$W"
        kallisto index $ARGS -i kallisto "${FA[0]}"
        mv kallisto "$OUT"/kallisto
    """,
)

# The reads and index are linked into a working directory under their own
# names, so kallisto's log names them as upstream's does.
kallisto_quant_cmd = shell.create(
    name="kallisto_quant",
    image=KALLISTO_IMAGE,
    resources=QUANT_RESOURCES,
    inputs={
        "reads_1": File,
        "reads_2": File | None,
        "index": File,
        "gtf": File | None,
        "chromosomes": File | None,
        "fragment_length": int,
        "fragment_length_sd": int,
        "prefix": str,
        "args": str,
    },
    defaults={"fragment_length": 0, "fragment_length_sd": 0, "args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        R1=({{inputs.reads_1}}); IDX=({{inputs.index}})
        shopt -s nullglob
        R2=({{inputs.reads_2}}); GTF=({{inputs.gtf}}); CHROMS=({{inputs.chromosomes}})
        shopt -u nullglob
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        ln -s "${{IDX[0]}}" "$(basename "${{IDX[0]}}")"
        READS=("$(basename "${{R1[0]}}")"); ln -s "${{R1[0]}}" "${{READS[0]}}"
        EXTRA=()
        if [ ${{#R2[@]}} -gt 0 ]; then
            READS+=("$(basename "${{R2[0]}}")"); ln -s "${{R2[0]}}" "${{READS[1]}}"
        else
            EXTRA+=(--single --fragment-length={{inputs.fragment_length}} --sd={{inputs.fragment_length_sd}})
        fi
        if [ ${{#GTF[@]}} -gt 0 ]; then EXTRA+=(--gtf "${{GTF[0]}}"); fi
        if [ ${{#CHROMS[@]}} -gt 0 ]; then EXTRA+=(--chromosomes "${{CHROMS[0]}}"); fi
        mkdir -p "$P"
        kallisto quant \
            --threads $TASK_CPUS \
            --index "$(basename "${{IDX[0]}}")" \
            "${{EXTRA[@]}}" \
            $ARGS \
            -o "$P" \
            "${{READS[@]}}" 2> "$P"/kallisto_quant.log \
            || {{ cat "$P"/kallisto_quant.log >&2; exit 1; }}
        cat "$P"/kallisto_quant.log >&2
        cp "$P"/kallisto_quant.log "$P".log
        cp "$P"/run_info.json "$P".run_info.json
        mv "$P" "$P".log "$P".run_info.json "$OUT"/
    """,
)

env = flyte.TaskEnvironment.from_task("kallisto", kallisto_index_cmd.as_task(), kallisto_quant_cmd.as_task())


async def kallisto_index(fasta: File, args: str = "") -> File:
    """The kallisto index of ``fasta``: one file named ``kallisto``."""
    results = await kallisto_index_cmd(fasta=fasta, args=args)
    index = await results.get_file("kallisto")
    if index is None:
        raise FileNotFoundError("kallisto index wrote no kallisto")
    return index


@dataclass
class KallistoQuantResult:
    results: Dir  # <prefix>: abundance.h5, abundance.tsv, kallisto_quant.log, run_info.json
    json_info: File  # <prefix>.run_info.json
    log: File  # <prefix>.log


async def kallisto_quant(
    reads_1: File,
    index: File,
    reads_2: File | None = None,
    gtf: File | None = None,
    chromosomes: File | None = None,
    fragment_length: int | None = None,
    fragment_length_sd: int | None = None,
    prefix: str = "sample",
    args: str = "",
) -> KallistoQuantResult:
    """Quantify one sample's reads (``reads_2`` None for single-end) against a kallisto ``index``."""
    if reads_2 is None and (fragment_length is None or fragment_length_sd is None):
        raise ValueError("kallisto quant needs fragment_length and fragment_length_sd for single-end reads")
    results = await kallisto_quant_cmd(
        reads_1=reads_1,
        reads_2=reads_2,
        index=index,
        gtf=gtf,
        chromosomes=chromosomes,
        fragment_length=fragment_length or 0,
        fragment_length_sd=fragment_length_sd or 0,
        prefix=prefix,
        args=args,
    )

    async def pick(name: str) -> File:
        found = await results.get_file(name)
        if found is None:
            raise FileNotFoundError(f"kallisto quant wrote no {name}")
        return found

    return KallistoQuantResult(
        results=Dir.from_existing_remote(results.path.rstrip("/") + "/" + prefix),
        json_info=await pick(f"{prefix}.run_info.json"),
        log=await pick(f"{prefix}.log"),
    )
