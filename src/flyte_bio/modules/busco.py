"""busco — assembly, proteome and transcriptome completeness with BUSCO.

Exposes two functions:

- :func:`busco_download` fetches a BUSCO lineage dataset (``bacteria_odb10``,
  ``eukaryota_odb10``, ...) as ``<lineage>.tar.gz``, an archive of BUSCO's
  ``--download_path`` (``lineages/<lineage>``, ``file_versions.tsv``, ...).
  It is a task of its own, so each lineage is downloaded once and cached.
  It's an archive rather than a Dir because lineages hold hundreds of files
  (``eukaryota_odb10`` has over 500), more than a task's input staging copes
  with.
- :func:`busco` checks one or more FASTA files (each may be gzipped) against a
  lineage in ``genome``, ``proteins`` or ``transcriptome`` mode. ``lineage``
  may also be ``auto``, ``auto_prok`` or ``auto_euk`` for BUSCO's
  ``--auto-lineage*``. With ``lineages`` (a :func:`busco_download` archive,
  or any ``.tar.gz`` of a ``--download_path``) BUSCO reads datasets from
  there; without, it downloads what it needs itself and the downloads are
  returned too.

Outputs keep upstream's names, under the prefix ``<prefix>`` (default
``sample-<lineage>``):

- ``<prefix>-busco.batch_summary.txt``: one line per input (always), and
  ``<prefix>-busco.batch_summary.failed.txt`` listing failed runs, if any;
- ``short_summary.*.{txt,json}``: one pair per input with genes found;
- ``<prefix>-busco.log``: BUSCO's log;
- ``<prefix>-busco/``: BUSCO's output directory, with per input
  ``<input>/run_<lineage>/{full_table.tsv,missing_busco_list.tsv,busco_sequences/}``
  and, in transcriptome mode, ``single_copy_proteins.faa`` and
  ``<input>/translated_proteins/``;
- ``busco_downloads/``: the lineages BUSCO downloaded, without ``lineages``.

``clean_intermediates`` removes BUSCO's auto-lineage, miniprot, hmmer and
BBTools intermediates and temporary files, as upstream. Pass ``--tar`` in
``args`` to have BUSCO tar its sequence directories. BBTools gets a quarter of
the task's memory, as upstream (at least 120 MiB).
"""

from dataclasses import dataclass
from fnmatch import fnmatchcase

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

BUSCO_IMAGE = "community.wave.seqera.io/library/busco:6.1.0--c502c35c2ab05540"
# Upstream labels BUSCO process_medium; large eukaryote genomes may need more memory.
DEFAULT_RESOURCES = flyte.Resources(cpu=4, memory="16Gi")
DOWNLOAD_RESOURCES = flyte.Resources(cpu=1, memory="2Gi")
MODES = ("genome", "proteins", "transcriptome")
AUTO_LINEAGES = {"auto": "--auto-lineage", "auto_prok": "--auto-lineage-prok", "auto_euk": "--auto-lineage-euk"}

busco_download_cmd = shell.create(
    name="busco_download",
    image=BUSCO_IMAGE,
    resources=DOWNLOAD_RESOURCES,
    inputs={"lineage": str},
    outputs={"results": Dir},
    script=r"""
        LINEAGE={inputs.lineage}
        OUT={outputs.results}
        W=$(mktemp -d); cd "$W"
        busco --download_path "$W/downloads" --download "$LINEAGE"
        if [ ! -d "$W/downloads/lineages/$LINEAGE" ]; then
            echo "busco downloaded no lineages/$LINEAGE" >&2
            exit 1
        fi
        tar -czf "$OUT/$LINEAGE.tar.gz" -C "$W/downloads" .
    """,
)

busco_cmd = shell.create(
    name="busco",
    image=BUSCO_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "fasta": list[File],
        "mode": str,
        "lineage_option": str,
        "lineages": File | None,
        "config_file": File | None,
        "clean_intermediates": bool,
        "prefix": str,
        "args": str,
    },
    defaults={"clean_intermediates": False, "args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        FASTA=({{inputs.fasta}})
        shopt -s nullglob; CONFIG=({{inputs.config_file}}); LINEAGES=({{inputs.lineages}}); shopt -u nullglob
        MODE={{inputs.mode}}
        LINEAGE_OPTION={{inputs.lineage_option}}
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"

        BBTOOLS_MB=$(( TASK_MEMORY_MB / 4 ))
        if [ $BBTOOLS_MB -ge 1024 ]; then
            export BUSCO_BBTOOLS_MEMORY=$(( BBTOOLS_MB / 1024 ))g
        elif [ $BBTOOLS_MB -gt 120 ]; then
            export BUSCO_BBTOOLS_MEMORY=${{BBTOOLS_MB}}m
        else
            export BUSCO_BBTOOLS_MEMORY=120m
        fi

        # Make Augustus's config available and writable, as upstream.
        ENV_AUGUSTUS=/opt/conda/etc/conda/activate.d/augustus.sh
        set +u
        if [ -z "${{AUGUSTUS_CONFIG_PATH}}" ] && [ -f "${{ENV_AUGUSTUS}}" ]; then
            source "${{ENV_AUGUSTUS}}"
        fi
        if [ ! -w "${{AUGUSTUS_CONFIG_PATH}}" ]; then
            AUG_CONF_DIR=$(mktemp -d)
            cp -r "$AUGUSTUS_CONFIG_PATH"/* "$AUG_CONF_DIR"
            export AUGUSTUS_CONFIG_PATH=$AUG_CONF_DIR
            echo "New AUGUSTUS_CONFIG_PATH=${{AUGUSTUS_CONFIG_PATH}}"
        fi

        # BUSCO reads a directory of inputs, uncompressed.
        INPUT_SEQS=$(mktemp -d)
        for F in "${{FASTA[@]}}"; do
            NAME=$(basename "$F")
            if [ "${{NAME##*.}}" = gz ]; then
                gzip -cdf "$F" > "$INPUT_SEQS/${{NAME%.gz}}"
            else
                ln -s "$F" "$INPUT_SEQS/$NAME"
            fi
        done

        EXTRA=()
        # Extracted to a writable directory: BUSCO checks its datasets for
        # updates and may rewrite them.
        if [ ${{#LINEAGES[@]}} -gt 0 ]; then
            DOWNLOADS=$(mktemp -d)
            tar -xzf "${{LINEAGES[0]}}" -C "$DOWNLOADS"
            EXTRA+=(--download_path "$DOWNLOADS")
        fi
        if [ ${{#CONFIG[@]}} -gt 0 ]; then EXTRA+=(--config "${{CONFIG[0]}}"); fi

        busco \
            --cpu $TASK_CPUS \
            --in "$INPUT_SEQS" \
            --out "$P-busco" \
            --mode "$MODE" \
            $LINEAGE_OPTION \
            "${{EXTRA[@]}}" \
            $ARGS

        rm -rf "$INPUT_SEQS"
        if [ {{inputs.clean_intermediates}} = true ]; then
            rm -fr ./*-busco/*/auto_lineage ./*-busco/*/*/{{miniprot,hmmer,.bbtools}}_output \
                ./*-busco/*/prodigal_output/predicted_genes/tmp/ ./*-busco/*/tmp/
        fi
        # Remove symlinks the cleanup broke.
        find . -xtype l -delete

        grep -v 'Run failed; check logs' "$P-busco/batch_summary.txt" > "$P-busco.batch_summary.txt"
        mv "$P-busco"/*/short_summary.*.{{json,txt}} . \
            || echo "Short summaries were not available: No genes were found."
        mv "$P-busco/logs/busco.log" "$P-busco.log"
        if grep -q 'Run failed; check logs' "$P-busco/batch_summary.txt"; then
            grep 'Run failed; check logs' "$P-busco/batch_summary.txt" > "$P-busco.batch_summary.failed.txt"
        fi
        mv ./* "$OUT"/
    """,
)

env = flyte.TaskEnvironment.from_task("busco", busco_download_cmd.as_task(), busco_cmd.as_task())


@dataclass
class BuscoResult:
    batch_summary: File  # <prefix>-busco.batch_summary.txt
    batch_summary_failed: File | None  # <prefix>-busco.batch_summary.failed.txt, if any run failed
    short_summaries_txt: list[File]  # short_summary.*.txt, one per input with genes found
    short_summaries_json: list[File]  # short_summary.*.json
    log: File | None  # <prefix>-busco.log
    full_table: list[File]  # <prefix>-busco/<input>/run_<lineage>/full_table.tsv, sorted by path
    missing_busco_list: list[File]  # .../run_<lineage>/missing_busco_list.tsv
    single_copy_proteins: list[File]  # .../run_<lineage>/single_copy_proteins.faa (transcriptome mode)
    seq_dir: list[Dir]  # .../run_<lineage>/busco_sequences
    translated_dir: list[Dir]  # <prefix>-busco/<input>/translated_proteins (transcriptome mode)
    single_copy_faa: list[File]  # .../busco_sequences/single_copy_busco_sequences/*.faa (without --tar)
    single_copy_fna: list[File]  # .../busco_sequences/single_copy_busco_sequences/*.fna (without --tar)
    downloaded_lineages: list[Dir]  # busco_downloads/lineages/*, when BUSCO downloaded them itself
    busco_dir: Dir  # <prefix>-busco
    results: Dir  # everything above


def path_matches(path: str, pattern: str) -> bool:
    """Whether ``path`` matches ``pattern`` segment by segment (``*`` doesn't cross ``/``)."""
    parts, wanted = path.split("/"), pattern.split("/")
    return len(parts) == len(wanted) and all(fnmatchcase(p, w) for p, w in zip(parts, wanted, strict=True))


async def busco_download(lineage: str) -> File:
    """Download the BUSCO dataset ``lineage``: ``<lineage>.tar.gz``, to pass to :func:`busco` as ``lineages``."""
    results = await busco_download_cmd(lineage=lineage)
    archive = await results.get_file(f"{lineage}.tar.gz")
    if archive is None:
        raise FileNotFoundError(f"busco_download wrote no {lineage}.tar.gz")
    return archive


async def busco(
    fasta: list[File],
    mode: str,
    lineage: str,
    lineages: File | None = None,
    config_file: File | None = None,
    clean_intermediates: bool = False,
    prefix: str = "",
    args: str = "",
) -> BuscoResult:
    """Run BUSCO on ``fasta`` against ``lineage`` (or ``auto*``); ``prefix`` defaults to ``sample-<lineage>``."""
    if mode not in MODES:
        raise ValueError("Mode must be one of 'genome', 'proteins', or 'transcriptome'.")
    prefix = prefix or f"sample-{lineage}"
    lineage_option = AUTO_LINEAGES.get(lineage, f"--lineage_dataset {lineage}")
    results = await busco_cmd(
        fasta=fasta,
        mode=mode,
        lineage_option=lineage_option,
        lineages=lineages,
        config_file=config_file,
        clean_intermediates=clean_intermediates,
        prefix=prefix,
        args=args,
    )

    base = results.path.rstrip("/")
    files = sorted([f.path[len(base) :].lstrip("/") async for f in results.walk()])
    out_dir = f"{prefix}-busco"

    def matching(pattern: str) -> list[File]:
        return [File.from_existing_remote(f"{base}/{p}") for p in files if path_matches(p, pattern)]

    def dirs_holding(pattern: str) -> list[Dir]:
        depth = len(pattern.split("/"))
        parents = {"/".join(p.split("/")[:depth]) for p in files if len(p.split("/")) > depth}
        return [Dir.from_existing_remote(f"{base}/{d}") for d in sorted(parents) if path_matches(d, pattern)]

    def single(name: str) -> File | None:
        return File.from_existing_remote(f"{base}/{name}") if name in files else None

    batch_summary = single(f"{prefix}-busco.batch_summary.txt")
    if batch_summary is None:
        raise FileNotFoundError(f"busco wrote no {prefix}-busco.batch_summary.txt")
    sequences = f"{out_dir}/*/run_*/busco_sequences/single_copy_busco_sequences"
    return BuscoResult(
        batch_summary=batch_summary,
        batch_summary_failed=single(f"{prefix}-busco.batch_summary.failed.txt"),
        short_summaries_txt=matching("short_summary.*.txt"),
        short_summaries_json=matching("short_summary.*.json"),
        log=single(f"{prefix}-busco.log"),
        full_table=matching(f"{out_dir}/*/run_*/full_table.tsv"),
        missing_busco_list=matching(f"{out_dir}/*/run_*/missing_busco_list.tsv"),
        single_copy_proteins=matching(f"{out_dir}/*/run_*/single_copy_proteins.faa"),
        seq_dir=dirs_holding(f"{out_dir}/*/run_*/busco_sequences"),
        translated_dir=dirs_holding(f"{out_dir}/*/translated_proteins"),
        single_copy_faa=matching(f"{sequences}/*.faa"),
        single_copy_fna=matching(f"{sequences}/*.fna"),
        downloaded_lineages=dirs_holding("busco_downloads/lineages/*"),
        busco_dir=Dir.from_existing_remote(f"{base}/{out_dir}"),
        results=results,
    )
