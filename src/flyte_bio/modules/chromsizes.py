"""chromsizes — chromosome sizes of a FASTA, via ``samtools faidx``.

Exposes :func:`chrom_sizes`, which indexes a FASTA (plain or bgzipped) and
returns its chromosome sizes (the index's first two columns) with the index
itself. Outputs keep upstream's names, after the FASTA's own name:

- ``<fasta name>.sizes``: one ``<name>\\t<length>`` line per sequence;
- ``<fasta name>.fai``: the FASTA index;
- ``<fasta name>.gzi``: the bgzip index, for a bgzipped FASTA only.

This follows the upstream module's last working version; later upstream
versions only fail with a deprecation notice pointing at ``samtools faidx``.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

# The image the upstream module pins.
CHROMSIZES_IMAGE = "community.wave.seqera.io/library/htslib_samtools:1.22.1--ef1e734834ad9fa5"
# Upstream's process_single.
DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory="2Gi")

# faidx writes its index next to the FASTA, so link it into a writable work dir.
chrom_sizes_cmd = shell.create(
    name="chrom_sizes",
    image=CHROMSIZES_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"fasta": File},
    outputs={"results": Dir},
    script=r"""
        FA=({inputs.fasta})
        NAME=$(basename "${FA[0]}")
        W=$(mktemp -d); cd "$W"
        ln -s "${FA[0]}" "$NAME"
        samtools faidx "$NAME"
        cut -f 1,2 "$NAME.fai" > "$NAME.sizes"
        shopt -s nullglob
        mv ./*.sizes ./*.fai ./*.gzi {outputs.results}/
    """,
)

env = flyte.TaskEnvironment.from_task("chromsizes", chrom_sizes_cmd.as_task())


@dataclass
class ChromSizesResult:
    sizes: File
    fai: File
    gzi: File | None  # for a bgzipped FASTA
    results: Dir  # everything faidx wrote


async def chrom_sizes(fasta: File) -> ChromSizesResult:
    """Chromosome sizes and the ``.fai`` (and, if bgzipped, ``.gzi``) index of ``fasta``."""
    results = await chrom_sizes_cmd(fasta=fasta)
    files = [f async for f in results.walk()]
    found = {suffix: [f for f in files if f.path.endswith(suffix)] for suffix in (".sizes", ".fai", ".gzi")}
    for suffix in (".sizes", ".fai"):
        if len(found[suffix]) != 1:
            raise FileNotFoundError(f"samtools faidx wrote {len(found[suffix])} *{suffix} files, expected one")
    return ChromSizesResult(
        sizes=found[".sizes"][0],
        fai=found[".fai"][0],
        gzi=found[".gzi"][0] if found[".gzi"] else None,
        results=results,
    )
