"""samtools — read/write and query high-throughput sequencing alignments.

Exposes the subset of the samtools suite the RNA-seq pipeline leans on:

- :data:`samtools_faidx` — index a FASTA (and emit chromosome sizes).
- :data:`samtools_sort` — coordinate-sort an alignment.
- :data:`samtools_index` — build a BAM index (``.bai``).
- :data:`samtools_stats` — full alignment statistics.
- :data:`samtools_flagstat` — FLAG-field tallies.
- :data:`samtools_idxstats` — per-reference mapped/unmapped counts.

All tasks share one htslib+samtools biocontainer, so they live in a single
:class:`flyte.TaskEnvironment`.
"""



import flyte
from flyte.extras import shell
from flyte.io import File

# Pinned biocontainer URI (matches the version the RNA-seq pipeline pins).
SAMTOOLS_IMAGE = "community.wave.seqera.io/library/htslib_samtools:1.23.1--5b6bb4ede7e612e5"

# Most of these are single-threaded I/O over modest inputs; sort is the one
# memory/throughput-sensitive step.
DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory="6Gi")
SORT_RESOURCES = flyte.Resources(cpu=4, memory="12Gi")


# `faidx` writes the index next to its input, so symlink the FASTA into the
# (writable) work dir first — the inputs mount may be read-only. The `.fai`
# columns are derived from sequence content + line wrapping, not the file
# name, so the md5 is reproducible. Chromosome sizes are just its first two
# columns; emit them too (the pipeline needs them for coverage tracks).
samtools_faidx = shell.create(
    name="samtools_faidx",
    image=SAMTOOLS_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"fasta": File},
    outputs={"fai": File, "sizes": File},
    script=r"""
        ln -s {inputs.fasta} genome.fa
        samtools faidx genome.fa
        cut -f 1,2 genome.fa.fai > genome.sizes
        mv genome.fa.fai {outputs.fai}
        mv genome.sizes {outputs.sizes}
    """,
)


samtools_sort = shell.create(
    name="samtools_sort",
    image=SAMTOOLS_IMAGE,
    resources=SORT_RESOURCES,
    inputs={"bam": File},
    outputs={"bam_sorted": File},
    script=r"""
        samtools sort -o {outputs.bam_sorted} {inputs.bam}
    """,
)


# Index reads the (coordinate-sorted) BAM and writes a sidecar. We name the
# output explicitly so it lands in the outputs mount rather than next to the
# read-only input.
samtools_index = shell.create(
    name="samtools_index",
    image=SAMTOOLS_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"bam": File},
    outputs={"bai": File},
    script=r"""
        samtools index {inputs.bam} {outputs.bai}
    """,
)


samtools_stats = shell.create(
    name="samtools_stats",
    image=SAMTOOLS_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"bam": File},
    outputs={"stats": File},
    script=r"""
        samtools stats {inputs.bam} > {outputs.stats}
    """,
)


samtools_flagstat = shell.create(
    name="samtools_flagstat",
    image=SAMTOOLS_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"bam": File},
    outputs={"flagstat": File},
    script=r"""
        samtools flagstat {inputs.bam} > {outputs.flagstat}
    """,
)


# idxstats needs the index present alongside the BAM under a matching name,
# so stage both under one basename before querying.
samtools_idxstats = shell.create(
    name="samtools_idxstats",
    image=SAMTOOLS_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"bam": File, "bai": File},
    outputs={"idxstats": File},
    script=r"""
        ln -s {inputs.bam} aln.bam
        ln -s {inputs.bai} aln.bam.bai
        samtools idxstats aln.bam > {outputs.idxstats}
    """,
)


env = flyte.TaskEnvironment.from_task(
    "samtools",
    samtools_faidx.as_task(),
    samtools_sort.as_task(),
    samtools_index.as_task(),
    samtools_stats.as_task(),
    samtools_flagstat.as_task(),
    samtools_idxstats.as_task(),
)
