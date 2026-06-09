"""cat — file concatenation utilities.

Currently exposes:

- :data:`cat_fastq` — concatenate gzipped FASTQ shards into a single
  merged file, preserving input order.

The wrapper assumes inputs are already gzipped (.fastq.gz): when the
first input ends in ``.gz`` it just calls ``cat`` byte-for-byte without
recompression — which keeps the output MD5 stable and snapshot-comparable.
"""



import flyte
from flyte.extras import shell
from flyte.io import File

# Seqera wave coreutils image.
CAT_IMAGE = "community.wave.seqera.io/library/coreutils_grep_gzip_lbzip2_pruned:838ba80435a629f8"

# Single-threaded; fixed modest resources.
DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory="6Gi")


cat_fastq = shell.create(
    name="cat_fastq",
    image=CAT_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"reads": list[File]},
    outputs={"merged": File},
    script=r"""
        ls -1 -v {inputs.reads} | xargs cat > {outputs.merged}
    """,
)


env = flyte.TaskEnvironment.from_task(
    "cat",
    cat_fastq.as_task(),
)
