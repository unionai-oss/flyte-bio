"""minia — de Bruijn graph assembly of short reads with a small memory footprint.

Exposes :func:`minia`, which assembles one sample's reads (one or more
FASTQ files, listed for minia in ``input_files.txt`` as upstream) and
returns the contigs ``<prefix>.contigs.fa.gz`` and unitigs
``<prefix>.unitigs.fa.gz`` (gzipped without names or dates), the graph
``<prefix>.h5`` and minia's log ``<prefix>-minia.log``. ``args`` carries
minia's options (e.g. ``-kmer-size 31``); its cores follow the task's CPUs.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

MINIA_IMAGE = "community.wave.seqera.io/library/minia:3.2.6--df502ab09998dab4"
# Upstream: process_high.
DEFAULT_RESOURCES = flyte.Resources(cpu=4, memory="8Gi")

minia_cmd = shell.create(
    name="minia",
    image=MINIA_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"reads": list[File], "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        READS=({{inputs.reads}})
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        printf '%s\n' "${{READS[@]}}" > input_files.txt
        minia \
            $ARGS \
            -nb-cores $TASK_CPUS \
            -in input_files.txt \
            -out "$P" > "$P-minia.log" 2>&1 || {{ tail -n 20 "$P-minia.log" >&2; exit 1; }}
        if [ -f "$P.contigs.fa" ]; then gzip -cn "$P.contigs.fa" > "$OUT/$P.contigs.fa.gz"; fi
        if [ -f "$P.unitigs.fa" ]; then gzip -cn "$P.unitigs.fa" > "$OUT/$P.unitigs.fa.gz"; fi
        mv "$P.h5" "$P-minia.log" "$OUT"/
    """,
)

env = flyte.TaskEnvironment.from_task("minia", minia_cmd.as_task())


@dataclass
class MiniaResult:
    contigs: File  # <prefix>.contigs.fa.gz
    unitigs: File  # <prefix>.unitigs.fa.gz
    h5: File  # <prefix>.h5
    log: File  # <prefix>-minia.log


async def minia(reads: list[File], prefix: str = "sample", args: str = "") -> MiniaResult:
    """An assembly of ``reads`` (any number of FASTQ files)."""
    results = await minia_cmd(reads=reads, prefix=prefix, args=args)
    picked = []
    for name in (f"{prefix}.contigs.fa.gz", f"{prefix}.unitigs.fa.gz", f"{prefix}.h5", f"{prefix}-minia.log"):
        found = await results.get_file(name)
        if found is None:
            raise FileNotFoundError(f"minia wrote no {name}")
        picked.append(found)
    return MiniaResult(*picked)
