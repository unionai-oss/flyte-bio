"""Run the scrnaseq pipeline (simpleaf / alevin-fry) on a samplesheet.

The pipeline functions run inside your own task, which fans out to the
flyte-bio module tasks. So the task's image needs ``flyte-bio`` installed,
and its environment must ``depends_on`` the flyte-bio module environments
(which deploys every tool's task and image alongside it).

Run it on the upstream test profile's data (two 10x v2 samples, one with two
lanes, against mouse chr19)::

    DATA=https://github.com/nf-core/test-datasets/raw/scrnaseq
    curl -sLO $DATA/samplesheet-2-0.csv
    flyte run examples/scrnaseq.py scrnaseq_example \\
        --samplesheet samplesheet-2-0.csv \\
        --fasta $DATA/reference/GRCm38.p6.genome.chr19.fa \\
        --gtf $DATA/reference/gencode.vM19.annotation.chr19.gtf \\
        --protocol 10XV2 \\
        --skip_qcatch --skip_cellbender

(QCatch and CellBender need a realistic number of cells, so upstream's test
profile skips them; on real data leave them on.) For your own data, give
FASTQ paths the cluster can read (``s3://…``, ``gs://…``) and add
``--outdir s3://my-bucket/scrnaseq`` to publish the results there.
"""

import flyte
from flyte.io import Dir, File

from flyte_bio import env as bio_env
from flyte_bio.pipelines.scrnaseq import scrnaseq

# flyte-bio isn't on PyPI yet: install it from GitHub. To run your local
# checkout instead, use `.with_uv_project("pyproject.toml",
# project_install_mode="install_project")` from the repo root.
FLYTE_BIO = "flyte-bio @ https://github.com/unionai-oss/flyte-bio/archive/refs/heads/main.zip"

env = flyte.TaskEnvironment(
    name="scrnaseq_example",
    image=flyte.Image.from_debian_base().with_pip_packages(FLYTE_BIO),
    resources=flyte.Resources(cpu=2, memory="4Gi"),
    depends_on=[bio_env],
)


@env.task
async def scrnaseq_example(
    samplesheet: File,
    fasta: File,
    gtf: File,
    protocol: str,
    skip_qcatch: bool = False,
    skip_cellbender: bool = False,
    outdir: str | None = None,
) -> Dir:
    """Quantify every sample to count matrices; return the results in upstream's --outdir layout."""
    result = await scrnaseq(
        samplesheet,
        fasta=fasta,
        gtf=gtf,
        protocol=protocol,
        skip_qcatch=skip_qcatch,
        skip_cellbender=skip_cellbender,
        outdir=outdir,
        publish_results=True,
        # Any other upstream option is a keyword argument, e.g. cellbender_gpu=True to run CellBender
        # on a GPU, simpleaf_index=... to reuse an index, or save_reference=True.
    )
    assert result.outdir is not None
    return result.outdir
