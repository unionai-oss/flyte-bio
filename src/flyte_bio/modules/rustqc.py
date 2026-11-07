"""rustqc — RNA-seq alignment QC in one pass (a reimplementation of the rnaseq QC tools).

Exposes :func:`rustqc`, which runs ``rustqc rna`` on a coordinate-sorted BAM
(with its index) and a GTF and returns its output directory ``<prefix>``, laid
out as upstream by the tool each part reimplements:

- ``dupradar/``: duplication-rate matrix and plots;
- ``featurecounts/``: ``<prefix>.featureCounts.tsv`` (+ ``.summary``), the
  biotype counts and their MultiQC tables;
- ``preseq/``: ``<prefix>.lc_extrap.txt``;
- ``samtools/``: ``<prefix>.{flagstat,idxstats,stats}``;
- ``rseqc/``: one subdirectory per RSeQC tool (bam_stat, infer_experiment,
  inner_distance, junction_annotation, junction_saturation,
  read_distribution, read_duplication, tin);
- ``qualimap/``: ``rnaseq_qc_results.txt``, the HTML report and its raw data
  and images.

Paired-end BAMs add ``--paired``, as upstream for samples not marked single-end.
"""

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

RUSTQC_IMAGE = "community.wave.seqera.io/library/rustqc:0.2.1--00df1502b490e005"
# Upstream labels it process_high.
DEFAULT_RESOURCES = flyte.Resources(cpu=4, memory="16Gi")

# The BAM and its index are linked side by side under their own names, so the
# reports name them as upstream's do.
rustqc_cmd = shell.create(
    name="rustqc",
    image=RUSTQC_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"bam": File, "bai": File, "gtf": File, "paired": bool, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        B=({{inputs.bam}}); I=({{inputs.bai}}); GTF=({{inputs.gtf}})
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        BAM=$(basename "${{B[0]}}"); ln -s "${{B[0]}}" "$BAM"; ln -s "${{I[0]}}" "$(basename "${{I[0]}}")"
        ln -s "${{GTF[0]}}" "$(basename "${{GTF[0]}}")"
        PAIRED=()
        if [ {{inputs.paired}} = true ]; then PAIRED=(--paired); fi
        rustqc rna \
            "$BAM" \
            --gtf "$(basename "${{GTF[0]}}")" \
            "${{PAIRED[@]}}" \
            --threads $TASK_CPUS \
            --outdir "$P" \
            --sample-name "$P" \
            $ARGS \
            || exit 1
        mv "$P" "$OUT"/
    """,
)

env = flyte.TaskEnvironment.from_task("rustqc", rustqc_cmd.as_task())


async def rustqc(bam: File, bai: File, gtf: File, paired: bool = True, prefix: str = "sample", args: str = "") -> Dir:
    """RNA-seq QC of a sorted, indexed ``bam`` against ``gtf``: the output directory ``<prefix>``."""
    results = await rustqc_cmd(bam=bam, bai=bai, gtf=gtf, paired=paired, prefix=prefix, args=args)
    return Dir.from_existing_remote(results.path.rstrip("/") + "/" + prefix)
