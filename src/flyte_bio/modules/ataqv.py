"""ataqv — ATAC-seq quality control metrics and their interactive report.

Exposes:

- :func:`ataqv_ataqv` — ATAC-seq metrics of one BAM (fragment lengths, TSS
  enrichment, reads in peaks, mitochondrial and duplicate fractions):
  ``<prefix>.ataqv.json``, and with ``--log-problematic-reads`` in ``args``,
  one ``<read group>.problems`` per read group.
- :func:`ataqv_mkarv` — the interactive HTML report of one or more ataqv
  metrics files, a Dir with ``index.html`` and its ``css``, ``fonts``, ``js``
  and ``data`` folders.

``organism`` (``human``, ``mouse``, ``fly``, ``worm``, ``rat`` or ``yeast``)
picks ataqv's autosomes, unless ``autosomal_reference_file`` lists them.
``peak_file``, ``tss_file``, ``excluded_regions_file`` and
``mitochondrial_reference_name`` add ``--peak-file``, ``--tss-file``,
``--excluded-region-file`` and ``--mitochondrial-reference-name``, as
upstream. TSS enrichment needs the BAM's index (``bai``). Inputs are linked
into the working directory under their own names, and both tools use the
task's CPUs.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

ATAQV_IMAGE = "quay.io/biocontainers/ataqv:1.3.1--py310ha155cf9_1"
DEFAULT_RESOURCES = flyte.Resources(cpu=4, memory="8Gi")

ataqv_ataqv_cmd = shell.create(
    name="ataqv_ataqv",
    image=ATAQV_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "bam": File,
        "bai": File | None,
        "peak_file": File | None,
        "organism": str,
        "mitochondrial_reference_name": str,
        "tss_file": File | None,
        "excluded_regions_file": File | None,
        "autosomal_reference_file": File | None,
        "prefix": str,
        "args": str,
    },
    defaults={"mitochondrial_reference_name": "", "args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        BAM=({{inputs.bam}})
        shopt -s nullglob
        BAI=({{inputs.bai}}); PEAK=({{inputs.peak_file}}); TSS=({{inputs.tss_file}})
        EXCLUDED=({{inputs.excluded_regions_file}}); AUTOSOMES=({{inputs.autosomal_reference_file}})
        shopt -u nullglob
        ORGANISM={{inputs.organism}}
        MITO={{inputs.mitochondrial_reference_name}}
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"

        link() {{ ln -s "$1" "$(basename "$1")"; }}
        link "${{BAM[0]}}"
        OPTS=()
        if [ ${{#BAI[@]}} -gt 0 ]; then link "${{BAI[0]}}"; fi
        if [ -n "$MITO" ]; then OPTS+=(--mitochondrial-reference-name "$MITO"); fi
        if [ ${{#PEAK[@]}} -gt 0 ]; then link "${{PEAK[0]}}"; OPTS+=(--peak-file "$(basename "${{PEAK[0]}}")"); fi
        if [ ${{#TSS[@]}} -gt 0 ]; then link "${{TSS[0]}}"; OPTS+=(--tss-file "$(basename "${{TSS[0]}}")"); fi
        if [ ${{#EXCLUDED[@]}} -gt 0 ]; then
            link "${{EXCLUDED[0]}}"; OPTS+=(--excluded-region-file "$(basename "${{EXCLUDED[0]}}")")
        fi
        if [ ${{#AUTOSOMES[@]}} -gt 0 ]; then
            link "${{AUTOSOMES[0]}}"; OPTS+=(--autosomal-reference-file "$(basename "${{AUTOSOMES[0]}}")")
        fi

        ataqv $ARGS "${{OPTS[@]}}" --metrics-file "$P.ataqv.json" --threads $TASK_CPUS --name "$P" \
            "$ORGANISM" "$(basename "${{BAM[0]}}")"
        shopt -s nullglob
        mv "$P.ataqv.json" ./*.problems "$OUT"/
    """,
)

# As upstream, the metrics files are linked into jsons/ under their own names.
ataqv_mkarv_cmd = shell.create(
    name="ataqv_mkarv",
    image=ATAQV_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"jsons": list[File], "args": str},
    defaults={"args": ""},
    outputs={"html": Dir},
    script=rf"""
        {TASK_LIMITS}
        JSONS=({{inputs.jsons}})
        ARGS={{inputs.args}}
        W=$(mktemp -d); cd "$W"
        mkdir jsons
        for JSON in "${{JSONS[@]}}"; do ln -s "$JSON" jsons/; done
        mkarv $ARGS --concurrency $TASK_CPUS --force ./html/ jsons/*
        cp -RL html/. {{outputs.html}}/
    """,
)

env = flyte.TaskEnvironment.from_task("ataqv", ataqv_ataqv_cmd.as_task(), ataqv_mkarv_cmd.as_task())


@dataclass
class AtaqvResult:
    json: File  # <prefix>.ataqv.json
    problems: list[File]  # <read group>.problems, with --log-problematic-reads
    results: Dir  # everything the task wrote


async def ataqv_ataqv(
    bam: File,
    organism: str,
    prefix: str = "sample",
    bai: File | None = None,
    peak_file: File | None = None,
    mitochondrial_reference_name: str = "",
    tss_file: File | None = None,
    excluded_regions_file: File | None = None,
    autosomal_reference_file: File | None = None,
    args: str = "",
) -> AtaqvResult:
    """ATAC-seq quality metrics of ``bam``."""
    results = await ataqv_ataqv_cmd(
        bam=bam,
        bai=bai,
        peak_file=peak_file,
        organism=organism,
        mitochondrial_reference_name=mitochondrial_reference_name,
        tss_file=tss_file,
        excluded_regions_file=excluded_regions_file,
        autosomal_reference_file=autosomal_reference_file,
        prefix=prefix,
        args=args,
    )
    json = await results.get_file(f"{prefix}.ataqv.json")
    if json is None:
        raise FileNotFoundError(f"ataqv wrote no {prefix}.ataqv.json")
    problems = sorted((f for f in await results.list_files() if f.name.endswith(".problems")), key=lambda f: f.name)
    return AtaqvResult(json=json, problems=problems, results=results)


async def ataqv_mkarv(jsons: list[File], args: str = "") -> Dir:
    """The interactive HTML report of ataqv metrics files ``jsons``."""
    if not jsons:
        raise ValueError("ataqv_mkarv needs metrics files")
    return await ataqv_mkarv_cmd(jsons=jsons, args=args)
