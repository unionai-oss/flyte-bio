"""untar — extract a tar archive.

Replicates the upstream strip-components heuristic: if every entry in the
archive shares a single top-level directory prefix, strip it; otherwise
extract as-is. That way the output is a flat list of payload files in
either layout.
"""



import flyte
from flyte.extras import shell
from flyte.extras.shell import Glob
from flyte.io import File

# Seqera wave coreutils image.
UNTAR_IMAGE = "community.wave.seqera.io/library/coreutils_grep_gzip_lbzip2_pruned:838ba80435a629f8"

# Single-threaded; fixed modest resources.
DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory="6Gi")


untar = shell.create(
    name="untar",
    image=UNTAR_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"archive": File},
    outputs={"files": Glob("**/*")},
    script=r"""
        if [[ $(tar -taf {inputs.archive} | grep -o -P "^.*?/" | uniq | wc -l) -eq 1 ]]; then
            tar -C {outputs.files} --strip-components 1 -xavf {inputs.archive}
        else
            tar -C {outputs.files} -xavf {inputs.archive}
        fi
    """,
)


env = flyte.TaskEnvironment.from_task(
    "untar",
    untar.as_task(),
)
