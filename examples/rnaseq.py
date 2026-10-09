"""Run the rnaseq pipeline (STAR + salmon) on a samplesheet.

The pipeline functions run inside your own task, which fans out to the
flyte-bio module tasks. So the task's image needs ``flyte-bio`` installed,
and its environment must ``depends_on`` the flyte-bio module environments
(which deploys every tool's task and image alongside it).

Run it on rnaseq's test profile data (5 small yeast samples)::

    DATA=https://raw.githubusercontent.com/nf-core/test-datasets/626c8fab639062eade4b10747e919341cbf9b41a
    curl -sO $DATA/samplesheet/v3.10/samplesheet_test.csv
    flyte run examples/rnaseq.py rnaseq_example \\
        --samplesheet samplesheet_test.csv \\
        --fasta $DATA/reference/genome.fasta \\
        --gtf $DATA/reference/genes_with_empty_tid.gtf.gz \\
        --transcript_fasta $DATA/reference/transcriptome.fasta

The samplesheet is passed as a local file (``flyte run`` uploads it) because
the task reads it itself, and reading ``https://`` files inside a task isn't
reliable; the FASTQ and reference URLs only go to the tool tasks, which
fetch ``https://`` fine. For your own data, use a samplesheet whose FASTQ
paths are URIs the cluster can read (``s3://…``, ``gs://…``).
"""

import flyte
from flyte.io import Dir, File

from flyte_bio import env as bio_env
from flyte_bio.pipelines.rnaseq import rnaseq

# flyte-bio isn't on PyPI yet: install it from GitHub. To run your local
# checkout instead, use `.with_uv_project("pyproject.toml",
# project_install_mode="install_project")` from the repo root.
FLYTE_BIO = "flyte-bio @ https://github.com/unionai-oss/flyte-bio/archive/refs/heads/main.zip"

env = flyte.TaskEnvironment(
    name="rnaseq_example",
    image=flyte.Image.from_debian_base().with_pip_packages(FLYTE_BIO),
    resources=flyte.Resources(cpu=2, memory="4Gi"),
    depends_on=[bio_env],
)


@env.task
async def rnaseq_example(
    samplesheet: File,
    fasta: File,
    gtf: File,
    transcript_fasta: File | None = None,
    outdir: str | None = None,
) -> Dir:
    """Quantify every sample; return the results in upstream's --outdir layout.

    With ``outdir`` (e.g. ``s3://my-bucket/rnaseq``) the results land there;
    without it, in Flyte's storage. Either way the returned Dir is that tree.
    """
    result = await rnaseq(
        samplesheet,
        fasta=fasta,
        gtf=gtf,
        transcript_fasta=transcript_fasta,
        outdir=outdir,
        publish_results=True,
        # Any other upstream option is a keyword argument, e.g.
        # skip_bbsplit=False with bbsplit_fasta_list={...}, with_umi=True with
        # umitools_bc_pattern="NNNNNN", or skip_qualimap=True.
    )
    assert result.outdir is not None
    return result.outdir
