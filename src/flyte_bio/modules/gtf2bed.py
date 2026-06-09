"""gtf2bed — derive a BED12 gene model from a GTF annotation.

The converter is a Perl helper, so it runs in a Perl biocontainer. Rather
than baking a custom image, we ship the real ``gtf2bed.pl`` as a Flyte File
input to a stock Perl image — Flyte stages input files into the container
just like any tool input. Input is a plain (uncompressed) GTF path; output
BED is written to the second argument.

The vendored script rides Flyte's default code bundle via the stub-module
trick in :mod:`flyte_bio.scripts` (no ``--copy-style all`` needed).
"""



import flyte
from flyte.extras import shell
from flyte.io import File

from flyte_bio.scripts import path

# Pinned biocontainer URI (a bare Perl interpreter). Upstream's main.nf names
# the bare `biocontainers/perl:5.26.2`, which resolves to Docker Hub and isn't
# anonymously pullable — use the quay.io biocontainers mirror instead.
GTF2BED_IMAGE = "quay.io/biocontainers/perl:5.26.2"

DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory="6Gi")

gtf2bed_cmd = shell.create(
    name="gtf2bed",
    image=GTF2BED_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"script": File, "gtf": File},
    outputs={"bed": File},
    script="perl {inputs.script} {inputs.gtf} {outputs.bed}\n",
)


env = flyte.TaskEnvironment.from_task(
    "gtf2bed",
    gtf2bed_cmd.as_task(),
)


async def gtf2bed(gtf: File) -> File:
    """Convert ``gtf`` to a BED12 gene model.

    Stages the vendored ``gtf2bed.pl`` as a File input, so the Perl image
    stays stock. Awaited from a task; the script is read from the code
    bundle in the caller's pod, so the run needs ``--copy-style all``.
    """
    script = await File.from_local(str(path("gtf2bed.pl")))
    return await gtf2bed_cmd(script=script, gtf=gtf)
