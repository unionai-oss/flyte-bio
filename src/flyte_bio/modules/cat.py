"""cat — file concatenation utilities.

Currently exposes:

- :func:`cat_fastq` — concatenate gzipped FASTQ shards into a single
  merged file, preserving input order.

Inputs are assumed already gzipped (.fastq.gz) and are joined with ``cat``
byte-for-byte without recompression (concatenated gzip members are a valid
gzip stream) — which keeps the output MD5 stable and snapshot-comparable.

Order matters: a sample's R1 and R2 shards must be merged in the same order
or read pairs stop lining up. A ``list[File]`` input can't carry that order
into the container — list items are staged under their original basenames
and a glob lists them alphabetically — so :func:`cat_fastq` folds the list
pairwise through :data:`cat_fastq_cmd`, whose two inputs are explicit.
"""



import flyte
from flyte.extras import shell
from flyte.io import File

# Seqera wave coreutils image.
CAT_IMAGE = "community.wave.seqera.io/library/coreutils_grep_gzip_lbzip2_pruned:838ba80435a629f8"

# Single-threaded; fixed modest resources.
DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory="6Gi")


cat_fastq_cmd = shell.create(
    name="cat_fastq",
    image=CAT_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"first": File, "second": File},
    outputs={"merged": File},
    script=r"""
        cat {inputs.first} {inputs.second} > {outputs.merged}
    """,
)


env = flyte.TaskEnvironment.from_task(
    "cat",
    cat_fastq_cmd.as_task(),
)


async def cat_fastq(reads: list[File]) -> File:
    """Concatenate ``reads`` in list order into one file.

    A single shard is returned as-is (no copy), matching upstream, which
    only merges samples that have more than one run.
    """
    if not reads:
        raise ValueError("cat_fastq needs at least one file")
    merged = reads[0]
    for shard in reads[1:]:
        merged = await cat_fastq_cmd(first=merged, second=shard)
    return merged
