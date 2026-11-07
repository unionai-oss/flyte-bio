"""mageck — MAGeCK, the CRISPR screen analysis toolkit.

Exposes:

- :func:`mageck_count` — count sgRNA reads per sample from FASTQs (each
  FASTQ a sample, labelled by ``sample_label``, a comma-separated list) or
  re-read an existing count table: ``<prefix>.count.txt`` and
  ``<prefix>.count_normalized.txt``, plus MAGeCK's ``<prefix>.countsummary.txt``
  and ``<prefix>.log``.
- :func:`mageck_mle` — gene essentiality by maximum likelihood over a count
  table and a design matrix: ``<prefix>.gene_summary.txt`` and
  ``<prefix>.sgrna_summary.txt``.
- :func:`mageck_test` — robust rank aggregation between treatment and
  control samples (``-t``/``-c`` or ``--day0-label`` in ``args``):
  ``*.gene_summary.txt`` and ``*.sgrna_summary.txt`` (one pair per
  comparison with ``--day0-label``) and the plotting script ``<prefix>.R``.

As upstream, :func:`mageck_count` passes FASTQs (inputs ending in
``.fastq.gz``) with ``--fastq`` and anything else as a count table with
``-k``. ``mageck mle`` runs on the task's CPUs.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

MAGECK_IMAGE = "quay.io/biocontainers/mageck:0.5.9.5--py39h1f90b4d_3"
DEFAULT_RESOURCES = flyte.Resources(cpu=4, memory="8Gi")

# A list[File] is staged under the files' names and globbed alphabetically, so
# ``names`` carries the caller's order (one name per line) into the script.
mageck_count_cmd = shell.create(
    name="mageck_count",
    image=MAGECK_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "inputs": list[File],
        "names": str,
        "library": File,
        "sample_label": str,
        "prefix": str,
        "args": str,
    },
    defaults={"sample_label": "", "args": ""},
    outputs={"results": Dir},
    script=r"""
        IN=({inputs.inputs})
        LIB=({inputs.library})
        NAMES={inputs.names}
        LABEL={inputs.sample_label}
        P={inputs.prefix}
        ARGS={inputs.args}
        OUT={outputs.results}
        D=$(dirname "${IN[0]}")
        W=$(mktemp -d); cd "$W"

        FILES=()
        while IFS= read -r NAME; do FILES+=("$D/$NAME"); done <<< "$NAMES"
        LAST=${FILES[-1]}
        if [[ $LAST == *.fastq.gz ]]; then INPUT=(--fastq "${FILES[@]}"); else INPUT=(-k "${FILES[@]}"); fi
        LABELS=()
        if [[ $LAST == *.fastq.gz || $LAST == *.fq.gz ]] && [ -n "$LABEL" ]; then LABELS=(--sample-label "$LABEL"); fi
        mageck count $ARGS -l "${LIB[0]}" -n "$P" "${LABELS[@]}" "${INPUT[@]}" || exit 1
        shopt -s nullglob
        mv ./*count.txt ./*.count_normalized.txt ./*.countsummary.txt ./*.log "$OUT"/
    """,
)

mageck_mle_cmd = shell.create(
    name="mageck_mle",
    image=MAGECK_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"count_table": File, "design_matrix": File, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        COUNTS=({{inputs.count_table}})
        DESIGN=({{inputs.design_matrix}})
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        mageck mle $ARGS --threads $TASK_CPUS -k "${{COUNTS[0]}}" -d "${{DESIGN[0]}}" -n "$P" || exit 1
        shopt -s nullglob
        mv ./*.gene_summary.txt ./*.sgrna_summary.txt ./*.log "$OUT"/
    """,
)

mageck_test_cmd = shell.create(
    name="mageck_test",
    image=MAGECK_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"count_table": File, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        COUNTS=({inputs.count_table})
        P={inputs.prefix}
        ARGS={inputs.args}
        OUT={outputs.results}
        W=$(mktemp -d); cd "$W"
        mageck test $ARGS -k "${COUNTS[0]}" -n "$P" || exit 1
        shopt -s nullglob
        mv ./*.gene_summary.txt ./*.sgrna_summary.txt ./*.R ./*.log ./*.pdf ./*.Rnw "$OUT"/
    """,
)

env = flyte.TaskEnvironment.from_task(
    "mageck", mageck_count_cmd.as_task(), mageck_mle_cmd.as_task(), mageck_test_cmd.as_task()
)


def file_name(file: File) -> str:
    return file.path.rstrip("/").rsplit("/", 1)[-1]


async def pick(results: Dir, name: str, tool: str) -> File:
    found = await results.get_file(name)
    if found is None:
        raise FileNotFoundError(f"{tool} wrote no {name}")
    return found


async def ending(results: Dir, suffix: str) -> list[File]:
    files = [f for f in await results.list_files() if file_name(f).endswith(suffix)]
    return sorted(files, key=file_name)


@dataclass
class MageckCountResult:
    count: File  # <prefix>.count.txt
    norm: File  # <prefix>.count_normalized.txt
    results: Dir  # also <prefix>.countsummary.txt and <prefix>.log


async def mageck_count(
    inputs: list[File], library: File, sample_label: str = "", prefix: str = "sample", args: str = ""
) -> MageckCountResult:
    """Count sgRNAs in ``inputs`` (FASTQs, one per sample, in order; or a count table) against ``library``.

    ``sample_label`` (e.g. ``"day0,day7"``) names the FASTQ samples, as
    upstream's ``meta.id``.
    """
    if not inputs:
        raise ValueError("mageck_count needs at least one input")
    names = [file_name(f) for f in inputs]
    if len(set(names)) != len(names):
        raise ValueError(f"mageck_count inputs need distinct names, got {names}")
    results = await mageck_count_cmd(
        inputs=inputs, names="\n".join(names), library=library, sample_label=sample_label, prefix=prefix, args=args
    )
    return MageckCountResult(
        count=await pick(results, f"{prefix}.count.txt", "mageck count"),
        norm=await pick(results, f"{prefix}.count_normalized.txt", "mageck count"),
        results=results,
    )


@dataclass
class MageckSummaryResult:
    gene_summary: list[File]  # *.gene_summary.txt, by name
    sgrna_summary: list[File]  # *.sgrna_summary.txt, by name
    r_script: File | None  # <prefix>.R (mageck test)
    results: Dir


async def mageck_mle(
    count_table: File, design_matrix: File, prefix: str = "sample", args: str = ""
) -> MageckSummaryResult:
    """``mageck mle`` over ``count_table`` with ``design_matrix``."""
    results = await mageck_mle_cmd(count_table=count_table, design_matrix=design_matrix, prefix=prefix, args=args)
    return MageckSummaryResult(
        gene_summary=[await pick(results, f"{prefix}.gene_summary.txt", "mageck mle")],
        sgrna_summary=[await pick(results, f"{prefix}.sgrna_summary.txt", "mageck mle")],
        r_script=None,
        results=results,
    )


async def mageck_test(count_table: File, prefix: str = "sample", args: str = "") -> MageckSummaryResult:
    """``mageck test`` over ``count_table``; ``args`` names the comparison (``-t``/``-c`` or ``--day0-label``)."""
    results = await mageck_test_cmd(count_table=count_table, prefix=prefix, args=args)
    gene_summary = await ending(results, ".gene_summary.txt")
    sgrna_summary = await ending(results, ".sgrna_summary.txt")
    if not gene_summary or not sgrna_summary:
        raise FileNotFoundError("mageck test wrote no gene or sgRNA summary")
    return MageckSummaryResult(
        gene_summary=gene_summary,
        sgrna_summary=sgrna_summary,
        r_script=await results.get_file(f"{prefix}.R"),
        results=results,
    )
