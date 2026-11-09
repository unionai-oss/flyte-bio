"""fusionreport — merge, rank and report gene fusions from several callers with fusion-report.

Exposes:

- :func:`fusionreport_download`, which downloads fusion-report's databases
  (``fusion_report download``): Mitelman, FusionGDB2 and COSMIC, as SQLite
  files in a Dir. COSMIC needs an account: pass ``--cosmic_usr <user>
  --cosmic_passwd <password>`` (and ``--qiagen`` for a commercial QIAGEN
  login) in ``args``, as upstream, or leave COSMIC out with ``--no-cosmic``
  (``--no-mitelman`` and ``--no-fusiongdb2`` leave out the others). Note that
  ``args`` is a task input, recorded with the run. ``fusionreport_detect``
  needs the FusionGDB2 database, whose server (``compbio.uth.edu``) was
  refusing scripted downloads (HTTP 429) on 2026-11-09.
- :func:`fusionreport_detect`, which merges the fusions found by Arriba,
  STAR-Fusion and FusionCatcher for one sample against those databases
  (``fusion_report run``, with ``--allow-multiple-gene-symbols`` and
  ``--tool-cutoff``, as upstream), writing ``<prefix>.fusionreport.tsv``,
  ``<prefix>.fusionreport_filtered.tsv``, the HTML report
  ``<prefix>_fusionreport_index.html``, MultiQC's ``fusion_genes_mqc.json``
  and, with ``--export csv`` / ``--export json`` in ``args``,
  ``<prefix>.fusions.csv`` / ``<prefix>.fusions.json``.

fusion-report writes one HTML page per fusion beside the index; a sample can
have hundreds, which Flyte stages poorly, so the pages (``*_*.html``, the
index among them, as upstream emits them) come back as one archive,
``<prefix>_fusionreport_html.tar.gz``.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

FUSIONREPORT_IMAGE = (
    "community.wave.seqera.io/library/fusion-report_beautifulsoup4_click_colorlog_pruned:15d1184d4eac76b8"
)
DEFAULT_RESOURCES = flyte.Resources(cpu=2, memory="8Gi")

# Fails the task with the end of fusion-report's log $1, which the shell task
# otherwise drops on failure: on stderr, and as the container's termination
# message, which the action's error shows.
REPORT_FAILURE = r"""
        report_failure() {
            tail -c 2000 "$1" >&2
            tail -c 2000 "$1" > /dev/termination-log 2> /dev/null
            exit 1
        }
""".strip()

fusionreport_download_cmd = shell.create(
    name="fusionreport_download",
    image=FUSIONREPORT_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"args": str},
    defaults={"args": ""},
    outputs={"db": Dir},
    script=rf"""
        {TASK_LIMITS}
        {REPORT_FAILURE}
        ARGS={{inputs.args}}
        W=$(mktemp -d); cd "$W"
        if ! fusion_report download $ARGS {{outputs.db}} > download.log 2>&1; then
            report_failure download.log
        fi
        cat download.log >&2
    """,
)

# As upstream: only the tools given are passed; the outputs are renamed to the prefix.
fusionreport_detect_cmd = shell.create(
    name="fusionreport_detect",
    image=FUSIONREPORT_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "arriba_fusions": File | None,
        "starfusion_fusions": File | None,
        "fusioncatcher_fusions": File | None,
        "db": Dir,
        "tools_cutoff": int,
        "prefix": str,
        "args": str,
    },
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        {REPORT_FAILURE}
        shopt -s nullglob
        ARRIBA=({{inputs.arriba_fusions}}); STARFUSION=({{inputs.starfusion_fusions}})
        FUSIONCATCHER=({{inputs.fusioncatcher_fusions}})
        shopt -u nullglob
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        TOOLS=()
        if [ ${{#ARRIBA[@]}} -gt 0 ]; then TOOLS+=(--arriba "${{ARRIBA[0]}}"); fi
        if [ ${{#STARFUSION[@]}} -gt 0 ]; then TOOLS+=(--starfusion "${{STARFUSION[0]}}"); fi
        if [ ${{#FUSIONCATCHER[@]}} -gt 0 ]; then TOOLS+=(--fusioncatcher "${{FUSIONCATCHER[0]}}"); fi
        if ! fusion_report run "$P" . \
            {{inputs.db}} "${{TOOLS[@]}}" \
            --allow-multiple-gene-symbols \
            --tool-cutoff {{inputs.tools_cutoff}} \
            $ARGS > run.log 2>&1; then
            report_failure run.log
        fi
        cat run.log >&2
        mv fusion_list.tsv "$OUT/$P.fusionreport.tsv"
        mv fusion_list_filtered.tsv "$OUT/$P.fusionreport_filtered.tsv"
        mv index.html "${{P}}_fusionreport_index.html"
        [ ! -f fusions.csv ] || mv fusions.csv "$OUT/$P.fusions.csv"
        [ ! -f fusions.json ] || mv fusions.json "$OUT/$P.fusions.json"
        shopt -s nullglob
        PAGES=(./*_*.html)
        tar -czf "$OUT/${{P}}_fusionreport_html.tar.gz" "${{PAGES[@]}}" || exit 1
        cp "${{P}}_fusionreport_index.html" "$OUT"/
        for f in ./*.json ./*.csv; do mv "$f" "$OUT"/; done
    """,
)

env = flyte.TaskEnvironment.from_task(
    "fusionreport",
    fusionreport_download_cmd.as_task(),
    fusionreport_detect_cmd.as_task(),
)


async def fusionreport_download(args: str = "") -> Dir:
    """fusion-report's databases (``*.db``), as a Dir; ``--no-cosmic`` or COSMIC credentials go in ``args``."""
    return await fusionreport_download_cmd(args=args)


@dataclass
class FusionReportResult:
    fusion_list: File  # <prefix>.fusionreport.tsv
    fusion_list_filtered: File  # <prefix>.fusionreport_filtered.tsv
    report: File  # <prefix>_fusionreport_index.html
    html: File  # <prefix>_fusionreport_html.tar.gz: the index and one page per fusion
    csv: File | None  # <prefix>.fusions.csv, with --export csv
    json: list[File]  # fusion_genes_mqc.json, and <prefix>.fusions.json with --export json
    results: Dir


async def fusionreport_detect(
    db: Dir,
    prefix: str,
    tools_cutoff: int = 1,
    arriba_fusions: File | None = None,
    starfusion_fusions: File | None = None,
    fusioncatcher_fusions: File | None = None,
    args: str = "",
) -> FusionReportResult:
    """Merge one sample's fusion calls against the fusion-report databases ``db``.

    A fusion is kept in the filtered list when at least ``tools_cutoff`` tools found it.
    """
    if arriba_fusions is None and starfusion_fusions is None and fusioncatcher_fusions is None:
        raise ValueError("fusionreport_detect needs at least one tool's fusions")
    results = await fusionreport_detect_cmd(
        arriba_fusions=arriba_fusions,
        starfusion_fusions=starfusion_fusions,
        fusioncatcher_fusions=fusioncatcher_fusions,
        db=db,
        tools_cutoff=tools_cutoff,
        prefix=prefix,
        args=args,
    )

    async def pick(name: str) -> File:
        found = await results.get_file(name)
        if found is None:
            raise FileNotFoundError(f"fusion_report wrote no {name}")
        return found

    files = sorted([f async for f in results.walk()], key=lambda f: f.name)
    return FusionReportResult(
        fusion_list=await pick(f"{prefix}.fusionreport.tsv"),
        fusion_list_filtered=await pick(f"{prefix}.fusionreport_filtered.tsv"),
        report=await pick(f"{prefix}_fusionreport_index.html"),
        html=await pick(f"{prefix}_fusionreport_html.tar.gz"),
        csv=await results.get_file(f"{prefix}.fusions.csv"),
        json=[f for f in files if f.name.endswith(".json")],
        results=results,
    )
