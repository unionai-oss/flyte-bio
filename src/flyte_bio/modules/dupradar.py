"""dupradar — duplication rate vs expression QC for RNA-seq BAMs.

Fits dupRadar's duplication-rate model per gene and writes the duplication
matrix, the intercept/slope summary, three PDF plots and two MultiQC custom
content tables (``<prefix>_dup_intercept_mqc.txt`` and
``<prefix>_duprateExpDensCurve_mqc.txt``). The R helper is vendored and
staged into the stock dupRadar container as a File input.

As upstream, the helper runs in the output directory with the sample id as
its prefix: it uses the prefix both for file names and as the sample name
inside the MultiQC tables, so it must not be a path.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.scripts import path

DUPRADAR_IMAGE = "community.wave.seqera.io/library/bioconductor-dupradar:1.38.0--831da16eb40a64ab"

DUPRADAR_THREADS = 4
DEFAULT_RESOURCES = flyte.Resources(cpu=DUPRADAR_THREADS, memory="8Gi")

dupradar_cmd = shell.create(
    name="dupradar",
    image=DUPRADAR_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "script": File,
        "bam": File,
        "gtf": File,
        "prefix": str,
        "strandedness": str,
        "single_end": bool,
        "args": str,
    },
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        BAM=({{inputs.bam}}); GTF=({{inputs.gtf}}); SCRIPT=({{inputs.script}})
        ARGS={{inputs.args}}
        cd {{outputs.results}}
        Rscript "${{SCRIPT[0]}}" \
            --bam "${{BAM[0]}}" \
            --gtf "${{GTF[0]}}" \
            --prefix {{inputs.prefix}} \
            --threads {DUPRADAR_THREADS} \
            --strandedness {{inputs.strandedness}} \
            --single_end {{inputs.single_end}} \
            $ARGS
    """,
)


env = flyte.TaskEnvironment.from_task(
    "dupradar",
    dupradar_cmd.as_task(),
)


@dataclass
class DupradarResult:
    dup_matrix: File  # <prefix>_dupMatrix.txt
    intercept_slope: File
    intercept_mqc: File  # <prefix>_dup_intercept_mqc.txt
    curve_mqc: File  # <prefix>_duprateExpDensCurve_mqc.txt
    results: Dir  # everything, incl. the PDF plots and R session log


async def dupradar(
    bam: File, gtf: File, prefix: str, strandedness: str, single_end: bool, args: str = ""
) -> DupradarResult:
    """Run dupRadar on ``bam``; ``strandedness`` is forward / reverse / anything else (unstranded).

    ``args`` takes the helper's options, e.g. ``--feature_type CDS`` for a GTF
    without exon features (the default feature type is ``exon``).
    """
    script = await File.from_local(str(path("dupradar.r")))
    results = await dupradar_cmd(
        script=script, bam=bam, gtf=gtf, prefix=prefix, strandedness=strandedness, single_end=single_end, args=args
    )

    async def pick(suffix: str) -> File:
        found = await results.get_file(f"{prefix}{suffix}")
        if found is None:
            raise FileNotFoundError(f"dupradar wrote no {prefix}{suffix}")
        return found

    return DupradarResult(
        dup_matrix=await pick("_dupMatrix.txt"),
        intercept_slope=await pick("_intercept_slope.txt"),
        intercept_mqc=await pick("_dup_intercept_mqc.txt"),
        curve_mqc=await pick("_duprateExpDensCurve_mqc.txt"),
        results=results,
    )
