"""rasusa — Rasusa, random subsampling of reads to a target depth.

Exposes :func:`rasusa`, which subsamples one sample's reads (single-end, or
an R1/R2 pair kept in step) to ``depth_cutoff``-fold coverage of a genome of
``genome_size`` (e.g. ``5mb``, or a FASTA/FAI to size it from). Further
Rasusa options (``--seed``, …) go in ``args``.

The outputs keep upstream's names: ``<prefix>.fastq.gz`` (single-end) or
``<prefix>_1.fastq.gz`` and ``<prefix>_2.fastq.gz`` (paired-end).

Rasusa exits 255 on error, which the platform would retry forever, so its
failures exit 1.
"""

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

RASUSA_IMAGE = "community.wave.seqera.io/library/rasusa:0.3.0--e22d5c77ebcbec79"
DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory="4Gi")

rasusa_cmd = shell.create(
    name="rasusa",
    image=RASUSA_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "reads_1": File,
        "reads_2": File | None,
        "genome_size": str,
        "depth_cutoff": str,
        "prefix": str,
        "args": str,
    },
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        R1=({inputs.reads_1})
        shopt -s nullglob; R2=({inputs.reads_2}); shopt -u nullglob
        GENOME_SIZE={inputs.genome_size}
        DEPTH={inputs.depth_cutoff}
        P={inputs.prefix}
        ARGS={inputs.args}
        OUT={outputs.results}
        if [ ${#R2[@]} -gt 0 ]; then
            OUTPUTS=(--output "$OUT/${P}_1.fastq.gz" "$OUT/${P}_2.fastq.gz")
        else
            OUTPUTS=(--output "$OUT/$P.fastq.gz")
        fi
        rasusa \
            $ARGS \
            --coverage "$DEPTH" \
            --genome-size "$GENOME_SIZE" \
            --input "${R1[0]}" "${R2[@]}" \
            "${OUTPUTS[@]}" \
            || exit 1
    """,
)

env = flyte.TaskEnvironment.from_task("rasusa", rasusa_cmd.as_task())


async def rasusa(
    reads: File | list[File], genome_size: str, depth_cutoff: int | float, prefix: str = "sample", args: str = ""
) -> list[File]:
    """Subsample ``reads`` (one file, or an R1/R2 pair) to ``depth_cutoff``-fold coverage of ``genome_size``.

    Returns the subsampled reads in the same order.
    """
    reads = reads if isinstance(reads, list) else [reads]
    if len(reads) not in (1, 2):
        raise ValueError(f"rasusa takes one or two read files, got {len(reads)}")
    results = await rasusa_cmd(
        reads_1=reads[0],
        reads_2=reads[1] if len(reads) == 2 else None,
        genome_size=genome_size,
        depth_cutoff=str(depth_cutoff),
        prefix=prefix,
        args=args,
    )
    names = [f"{prefix}_1.fastq.gz", f"{prefix}_2.fastq.gz"] if len(reads) == 2 else [f"{prefix}.fastq.gz"]
    out = []
    for name in names:
        found = await results.get_file(name)
        if found is None:
            raise FileNotFoundError(f"rasusa wrote no {name}")
        out.append(found)
    return out
