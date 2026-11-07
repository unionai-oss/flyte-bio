"""porechop — adapter trimming of Oxford Nanopore reads with Porechop or Porechop_ABI.

Exposes two functions:

- :func:`porechop` runs Porechop on one sample's reads, writing
  ``<prefix>.fastq.gz`` (default prefix ``sample``) and Porechop's report
  ``<prefix>.log`` (its stdout);
- :func:`porechop_abi` runs Porechop_ABI, which can also infer adapters ab
  initio (``--ab_initio`` in ``args``) or take a custom adapter file
  (``--custom_adapters``), writing ``<prefix>.fastq.gz`` (default prefix
  ``sample.porechop_abi``) and its report ``<prefix>.log``.

The two tools use different images, so each has an environment of its own
under the aggregate ``env``. Their ``-t``/``--threads`` is the task's CPUs.
"""

import os
from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

PORECHOP_IMAGE = "community.wave.seqera.io/library/porechop_pigz:d1655e5b5bad786c"
PORECHOP_ABI_IMAGE = "quay.io/biocontainers/porechop_abi:0.5.0post1--py310h275bdba_0"

# Upstream's label asks for 6 CPUs and 36 GB, far more than Porechop needs for
# a sample's reads; the 4 CPUs match upstream's tests.
DEFAULT_RESOURCES = flyte.Resources(cpu=4, memory="8Gi")

porechop_cmd = shell.create(
    name="porechop_porechop",
    image=PORECHOP_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"reads": File, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        READS=({{inputs.reads}})
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        porechop \
            -i "${{READS[0]}}" \
            -t $TASK_CPUS \
            $ARGS \
            -o "$P.fastq.gz" \
            > "$P.log"
        mv "$P.fastq.gz" "$P.log" "$OUT"/
    """,
)

porechop_abi_cmd = shell.create(
    name="porechop_abi",
    image=PORECHOP_ABI_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"reads": File, "custom_adapters": File | None, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        READS=({{inputs.reads}})
        shopt -s nullglob; ADAPTERS=({{inputs.custom_adapters}}); shopt -u nullglob
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        ADAPTER_ARGS=()
        if [ ${{#ADAPTERS[@]}} -gt 0 ]; then ADAPTER_ARGS=(--custom_adapters "${{ADAPTERS[0]}}"); fi
        W=$(mktemp -d); cd "$W"
        set -o pipefail
        porechop_abi \
            --input "${{READS[0]}}" \
            "${{ADAPTER_ARGS[@]}}" \
            --threads $TASK_CPUS \
            $ARGS \
            --output "$P.fastq.gz" \
            | tee "$P.log"
        mv "$P.fastq.gz" "$P.log" "$OUT"/
    """,
)

porechop_env = flyte.TaskEnvironment.from_task("porechop_porechop", porechop_cmd.as_task())
porechop_abi_env = flyte.TaskEnvironment.from_task("porechop_abi", porechop_abi_cmd.as_task())
env = flyte.TaskEnvironment(name="porechop", depends_on=[porechop_env, porechop_abi_env])


@dataclass
class PorechopResult:
    reads: File  # <prefix>.fastq.gz, the trimmed reads
    log: File  # <prefix>.log, the trimming report
    results: Dir  # both


async def read_porechop_results(results: Dir, prefix: str, tool: str) -> PorechopResult:
    async def pick(name: str) -> File:
        found = await results.get_file(name)
        if found is None:
            raise FileNotFoundError(f"{tool} wrote no {name}")
        return found

    return PorechopResult(reads=await pick(f"{prefix}.fastq.gz"), log=await pick(f"{prefix}.log"), results=results)


async def porechop(reads: File, prefix: str = "sample", args: str = "") -> PorechopResult:
    """Trim adapters from ``reads`` with Porechop."""
    results = await porechop_cmd(reads=reads, prefix=prefix, args=args)
    return await read_porechop_results(results, prefix, "porechop")


async def porechop_abi(
    reads: File,
    custom_adapters: File | None = None,
    prefix: str = "sample.porechop_abi",
    args: str = "",
) -> PorechopResult:
    """Trim adapters from ``reads`` with Porechop_ABI, optionally with ``custom_adapters``."""
    if (reads.name or os.path.basename(reads.path)) == f"{prefix}.fastq.gz":
        raise ValueError('Input and output names are the same, use "prefix" to disambiguate!')
    results = await porechop_abi_cmd(reads=reads, custom_adapters=custom_adapters, prefix=prefix, args=args)
    return await read_porechop_results(results, prefix, "porechop_abi")
