"""ncbigenomedownload — download genomes from NCBI with ncbi-genome-download.

Exposes :func:`ncbigenomedownload`, which downloads the assemblies of one or
more taxonomic ``groups`` (``bacteria``, ``viral``, ``all``, … or a
comma-separated list) from NCBI's RefSeq or GenBank at run time. Select
assemblies with a file of ``accessions`` (``-A``) or ``taxids`` (``-t``), or
with ncbi-genome-download's own options in ``args`` (``-A <accession>``,
``--formats genbank,fasta``, ``--section genbank``, …). The task is cached on
its inputs, so a download runs once; pass ``-N`` (no cache) for the tool's
own cache to stay out of the way.

Files are written flat (``--flat-output``) under NCBI's names,
``<assembly accession>_<assembly name>_<kind>``. The result groups them as
upstream's output globs do (``fna`` matches every ``*_genomic.fna.gz``, so it
includes ``cds`` and ``rna_fna``). ``--parallel`` follows the task's CPUs.
"""

from dataclasses import dataclass
from fnmatch import fnmatchcase

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

NCBIGENOMEDOWNLOAD_IMAGE = "quay.io/biocontainers/ncbi-genome-download:0.3.3--pyh7cba7a3_0"
# Upstream labels it process_low; the downloads are network-bound.
DEFAULT_RESOURCES = flyte.Resources(cpu=2, memory="4Gi")

ncbigenomedownload_cmd = shell.create(
    name="ncbigenomedownload",
    image=NCBIGENOMEDOWNLOAD_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"accessions": File | None, "taxids": File | None, "groups": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        shopt -s nullglob; ACCESSIONS=({{inputs.accessions}}); TAXIDS=({{inputs.taxids}}); shopt -u nullglob
        GROUP_LIST={{inputs.groups}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        SELECT=()
        if [ ${{#ACCESSIONS[@]}} -gt 0 ]; then SELECT+=(-A "${{ACCESSIONS[0]}}"); fi
        if [ ${{#TAXIDS[@]}} -gt 0 ]; then SELECT+=(-t "${{TAXIDS[0]}}"); fi

        ncbi-genome-download \
            $ARGS \
            "${{SELECT[@]}}" \
            --output-folder ./ \
            --flat-output \
            --parallel $TASK_CPUS \
            $GROUP_LIST || exit 1
        shopt -s nullglob
        FILES=(./*)
        if [ ${{#FILES[@]}} -gt 0 ]; then mv "${{FILES[@]}}" "$OUT"/; fi
    """,
)

env = flyte.TaskEnvironment.from_task("ncbigenomedownload", ncbigenomedownload_cmd.as_task())

# Upstream's output globs, by output name.
OUTPUT_GLOBS = {
    "gbk": "*_genomic.gbff.gz",
    "fna": "*_genomic.fna.gz",
    "rm": "*_rm.out.gz",
    "features": "*_feature_table.txt.gz",
    "gff": "*_genomic.gff.gz",
    "faa": "*_protein.faa.gz",
    "gpff": "*_protein.gpff.gz",
    "wgs_gbk": "*_wgsmaster.gbff.gz",
    "cds": "*_cds_from_genomic.fna.gz",
    "rna": "*_rna.fna.gz",
    "rna_fna": "*_rna_from_genomic.fna.gz",
    "report": "*_assembly_report.txt",
    "stats": "*_assembly_stats.txt",
}


@dataclass
class NcbiGenomeDownloadResult:
    gbk: list[File]  # *_genomic.gbff.gz
    fna: list[File]  # *_genomic.fna.gz (with the CDS and RNA from genomic)
    rm: list[File]  # *_rm.out.gz
    features: list[File]  # *_feature_table.txt.gz
    gff: list[File]  # *_genomic.gff.gz
    faa: list[File]  # *_protein.faa.gz
    gpff: list[File]  # *_protein.gpff.gz
    wgs_gbk: list[File]  # *_wgsmaster.gbff.gz
    cds: list[File]  # *_cds_from_genomic.fna.gz
    rna: list[File]  # *_rna.fna.gz
    rna_fna: list[File]  # *_rna_from_genomic.fna.gz
    report: list[File]  # *_assembly_report.txt
    stats: list[File]  # *_assembly_stats.txt
    results: Dir  # everything downloaded


async def ncbigenomedownload(
    groups: str,
    accessions: File | None = None,
    taxids: File | None = None,
    args: str = "",
) -> NcbiGenomeDownloadResult:
    """Download the NCBI assemblies of ``groups`` that ``accessions``, ``taxids`` and ``args`` select."""
    results = await ncbigenomedownload_cmd(accessions=accessions, taxids=taxids, groups=groups, args=args)
    files = sorted(await results.list_files(), key=lambda f: f.name)
    matched = {key: [f for f in files if fnmatchcase(f.name, pattern)] for key, pattern in OUTPUT_GLOBS.items()}
    return NcbiGenomeDownloadResult(**matched, results=results)
