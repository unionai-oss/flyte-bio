"""Run the variant_calling pipeline (BWA-MEM + GATK HaplotypeCaller) on a samplesheet.

The pipeline functions run inside your own task, which fans out to the
flyte-bio module tasks. So the task's image needs ``flyte-bio`` installed,
and its environment must ``depends_on`` the flyte-bio module environments
(which deploys every tool's task and image alongside it).

Run it on the upstream test profile's data (one sample, two lanes, a 40 kb
slice of chr22)::

    DATA=https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/genome
    curl -sO https://raw.githubusercontent.com/nf-core/sarek/3.10.0/tests/csv/3.0/fastq_single.csv
    flyte run examples/variant_calling.py variant_calling_example \\
        --samplesheet fastq_single.csv \\
        --fasta $DATA/genome.fasta \\
        --dbsnp $DATA/vcf/dbsnp_146.hg38.vcf.gz \\
        --known_indels $DATA/vcf/mills_and_1000G.indels.vcf.gz \\
        --intervals $DATA/genome.interval_list

(The test reads come from the whole genome, so on that tiny reference they
yield no calls after base-quality recalibration; real data does.) For your
own data, give FASTQ paths the cluster can read (``s3://…``, ``gs://…``) and
add ``--outdir s3://my-bucket/variant-calling`` to publish the results there.
"""

import flyte
from flyte.io import Dir, File

from flyte_bio import env as bio_env
from flyte_bio.pipelines.variant_calling import variant_calling

# flyte-bio isn't on PyPI yet: install it from GitHub. To run your local
# checkout instead, use `.with_uv_project("pyproject.toml",
# project_install_mode="install_project")` from the repo root.
FLYTE_BIO = "flyte-bio @ https://github.com/unionai-oss/flyte-bio/archive/refs/heads/main.zip"

env = flyte.TaskEnvironment(
    name="variant_calling_example",
    image=flyte.Image.from_debian_base().with_pip_packages(FLYTE_BIO),
    resources=flyte.Resources(cpu=2, memory="4Gi"),
    depends_on=[bio_env],
)


@env.task
async def variant_calling_example(
    samplesheet: File,
    fasta: File,
    dbsnp: File,
    known_indels: File,
    intervals: File | None = None,
    outdir: str | None = None,
) -> Dir:
    """Call germline variants for every sample; return the results in upstream's --outdir layout."""
    result = await variant_calling(
        samplesheet,
        fasta=fasta,
        dbsnp=dbsnp,
        known_indels=[known_indels],
        intervals=intervals,
        outdir=outdir,
        publish_results=True,
        # Any other upstream option is a keyword argument, e.g. wes=True for exome data,
        # skip_baserecalibrator=True, skip_haplotypecaller_filter=True or nucleotides_per_second=...
    )
    assert result.outdir is not None
    return result.outdir
