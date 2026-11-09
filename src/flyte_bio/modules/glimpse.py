"""glimpse — GLIMPSE (v1) genotype imputation and phasing of low-coverage data.

Exposes the three steps of GLIMPSE 1.1.1, which run in this order:

- :func:`glimpse_chunk` — split a region of a reference panel into imputation
  chunks, ``<prefix>.txt`` (one line per chunk: its input region with buffer
  and its output region).
- :func:`glimpse_phase` — impute and phase one chunk of a sample's genotype
  likelihoods against the reference panel (and optionally a genetic map),
  into ``<prefix>.<suffix>`` (upstream's default prefix is
  ``<id>_<input region>``, with ``:`` turned into ``_``).
- :func:`glimpse_ligate` — ligate the phased chunks (each with its index) into
  ``<prefix>.<suffix>``. As upstream, the chunks are listed in name order.

Inputs are linked into a working directory under their own names, indexes
beside their VCFs. ``--thread`` follows the task's CPUs.
"""

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

GLIMPSE_IMAGE = "quay.io/biocontainers/glimpse-bio:1.1.1--hce55b13_1"
# Upstream's process_medium for chunk and phase, at the 4 CPUs its tests run
# with; ligate is process_low.
MEDIUM_RESOURCES = flyte.Resources(cpu=4, memory="8Gi")
LOW_RESOURCES = flyte.Resources(cpu=2, memory="4Gi")

# `link FILE` links FILE into the working directory under its own name and
# prints that name; `link_index VCF INDEX` links INDEX beside VCF.
LINK_HELPERS = r"""
        link() { ln -s "$1" "$(basename "$1")"; basename "$1"; }
        link_index() { case "$2" in *.csi) ln -s "$2" "$1.csi" ;; *) ln -s "$2" "$1.tbi" ;; esac; }
        W=$(mktemp -d); cd "$W"
"""

glimpse_chunk_cmd = shell.create(
    name="glimpse_chunk",
    image=GLIMPSE_IMAGE,
    resources=MEDIUM_RESOURCES,
    inputs={"vcf": File, "index": File, "region": str, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        V=({{inputs.vcf}}); I=({{inputs.index}})
        REGION={{inputs.region}}; P={{inputs.prefix}}
        ARGS={{inputs.args}}
        {LINK_HELPERS}
        NAME=$(link "${{V[0]}}"); link_index "$NAME" "${{I[0]}}"
        GLIMPSE_chunk $ARGS --input "$NAME" --region "$REGION" --thread $TASK_CPUS --output "$P.txt" || exit 1
        mv "$P.txt" {{outputs.results}}/
    """,
)

glimpse_phase_cmd = shell.create(
    name="glimpse_phase",
    image=GLIMPSE_IMAGE,
    resources=MEDIUM_RESOURCES,
    inputs={
        "vcf": File,
        "index": File,
        "samples_file": File | None,
        "input_region": str,
        "output_region": str,
        "reference": File,
        "reference_index": File,
        "genetic_map": File | None,
        "prefix": str,
        "suffix": str,
        "args": str,
    },
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        V=({{inputs.vcf}}); I=({{inputs.index}}); REF=({{inputs.reference}}); REF_I=({{inputs.reference_index}})
        shopt -s nullglob; S=({{inputs.samples_file}}); M=({{inputs.genetic_map}}); shopt -u nullglob
        INPUT_REGION={{inputs.input_region}}; OUTPUT_REGION={{inputs.output_region}}
        P={{inputs.prefix}}; SUFFIX={{inputs.suffix}}
        ARGS={{inputs.args}}
        {LINK_HELPERS}
        NAME=$(link "${{V[0]}}"); link_index "$NAME" "${{I[0]}}"
        REF_NAME=$(link "${{REF[0]}}"); link_index "$REF_NAME" "${{REF_I[0]}}"
        EXTRA=()
        if [ ${{#M[@]}} -gt 0 ]; then EXTRA+=(--map "$(link "${{M[0]}}")"); fi
        if [ ${{#S[@]}} -gt 0 ]; then EXTRA+=(--samples-file "$(link "${{S[0]}}")"); fi
        GLIMPSE_phase $ARGS --input "$NAME" --reference "$REF_NAME" "${{EXTRA[@]}}" \
            --input-region "$INPUT_REGION" --output-region "$OUTPUT_REGION" \
            --thread $TASK_CPUS --output "$P.$SUFFIX" || exit 1
        mv "$P.$SUFFIX" {{outputs.results}}/
    """,
)

# As upstream, the chunks are listed (all_files.txt) in name order.
glimpse_ligate_cmd = shell.create(
    name="glimpse_ligate",
    image=GLIMPSE_IMAGE,
    resources=LOW_RESOURCES,
    inputs={"vcfs": list[File], "indexes": list[File], "prefix": str, "suffix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        shopt -s nullglob; VCFS=({{inputs.vcfs}}); IDX=({{inputs.indexes}}); shopt -u nullglob
        P={{inputs.prefix}}; SUFFIX={{inputs.suffix}}
        ARGS={{inputs.args}}
        {LINK_HELPERS}
        for F in "${{VCFS[@]}}" "${{IDX[@]}}"; do link "$F" > /dev/null; done
        for F in "${{VCFS[@]}}"; do basename "$F"; done | LC_ALL=C sort > all_files.txt
        GLIMPSE_ligate $ARGS --input all_files.txt --thread $TASK_CPUS --output "$P.$SUFFIX" || exit 1
        mv "$P.$SUFFIX" {{outputs.results}}/
    """,
)

env = flyte.TaskEnvironment.from_task(
    "glimpse", glimpse_chunk_cmd.as_task(), glimpse_phase_cmd.as_task(), glimpse_ligate_cmd.as_task()
)


async def pick(results: Dir, name: str, tool: str) -> File:
    found = await results.get_file(name)
    if found is None:
        raise FileNotFoundError(f"{tool} wrote no {name}")
    return found


async def glimpse_chunk(vcf: File, index: File, region: str, prefix: str, args: str = "") -> File:
    """``<prefix>.txt``: the imputation chunks of ``region`` (e.g. ``chr22``) of the panel ``vcf``.

    Chunk sizes go in ``args``, e.g. ``--window-size 2000000 --buffer-size 200000``.
    """
    results = await glimpse_chunk_cmd(vcf=vcf, index=index, region=region, prefix=prefix, args=args)
    return await pick(results, f"{prefix}.txt", "GLIMPSE_chunk")


async def glimpse_phase(
    vcf: File,
    index: File,
    reference: File,
    reference_index: File,
    input_region: str,
    output_region: str,
    prefix: str,
    samples_file: File | None = None,
    genetic_map: File | None = None,
    suffix: str = "vcf.gz",
    args: str = "",
) -> File:
    """Impute and phase ``vcf``'s genotype likelihoods in one chunk: ``<prefix>.<suffix>``.

    ``input_region`` (the chunk with its buffer) and ``output_region`` come
    from :func:`glimpse_chunk`. ``reference`` is the phased panel;
    ``samples_file`` gives each sample's ploidy (``--samples-file``) and
    ``genetic_map`` is passed as ``--map``. Upstream's default prefix is
    ``<id>_<input_region>`` with ``:`` replaced by ``_``.
    """
    results = await glimpse_phase_cmd(
        vcf=vcf,
        index=index,
        samples_file=samples_file,
        input_region=input_region,
        output_region=output_region,
        reference=reference,
        reference_index=reference_index,
        genetic_map=genetic_map,
        prefix=prefix,
        suffix=suffix,
        args=args,
    )
    return await pick(results, f"{prefix}.{suffix}", "GLIMPSE_phase")


async def glimpse_ligate(
    vcfs: list[File], indexes: list[File], prefix: str, suffix: str = "vcf.gz", args: str = ""
) -> File:
    """Ligate phased chunks ``vcfs`` (with their ``indexes``, named after them) into ``<prefix>.<suffix>``."""
    names = [f.name for f in vcfs]
    if len(set(names)) != len(names):
        raise ValueError(f"GLIMPSE_ligate inputs need distinct names, got {names}")
    results = await glimpse_ligate_cmd(vcfs=vcfs, indexes=indexes, prefix=prefix, suffix=suffix, args=args)
    return await pick(results, f"{prefix}.{suffix}", "GLIMPSE_ligate")
