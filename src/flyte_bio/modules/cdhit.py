"""cdhit — clustering of nucleotide sequences with CD-HIT.

Exposes :func:`cdhit_cdhitest`, which clusters a FASTA or FASTQ of
nucleotide sequences with ``cd-hit-est`` and returns the representative
sequences, ``<prefix>.fa`` (``<prefix>.fq`` for FASTQ input), and the
clusters, ``<prefix>.fa.clstr`` (``<prefix>.fq.clstr``). As upstream, input
whose name ends in ``.fa``/``.fasta``/``.fna``/``.fas`` (optionally
``.gz``) is taken as FASTA.

As upstream, ``-T`` is the task's CPUs and ``-M`` 80% of its memory, in MB;
options such as the identity threshold (``-c``) go in ``args``.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

CDHIT_IMAGE = "community.wave.seqera.io/library/cd-hit:4.8.1--2cbab405453b19db"
# Upstream: process_medium.
DEFAULT_RESOURCES = flyte.Resources(cpu=4, memory="8Gi")

cdhit_cdhitest_cmd = shell.create(
    name="cdhit_cdhitest",
    image=CDHIT_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"sequences": File, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        SEQ=({{inputs.sequences}})
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        SUFFIX=fq
        if [[ "$(basename "${{SEQ[0]}}")" =~ f[astn]*a(\.gz)?$ ]]; then SUFFIX=fa; fi
        W=$(mktemp -d); cd "$W"
        cd-hit-est \
            $ARGS \
            -i "${{SEQ[0]}}" \
            -o "$P.$SUFFIX" \
            -M $(( TASK_MEMORY_MB * 8 / 10 )) \
            -T $TASK_CPUS \
            || exit 1
        mv "$P.$SUFFIX" "$P.$SUFFIX.clstr" "$OUT"/
    """,
)

env = flyte.TaskEnvironment.from_task("cdhit", cdhit_cdhitest_cmd.as_task())


@dataclass
class CdhitResult:
    fasta: File  # the representative sequences: <prefix>.fa (or .fq)
    clusters: File  # <prefix>.fa.clstr (or .fq.clstr)


async def cdhit_cdhitest(sequences: File, prefix: str = "sample", args: str = "") -> CdhitResult:
    """Cluster the nucleotide ``sequences``: their representatives and the clusters."""
    results = await cdhit_cdhitest_cmd(sequences=sequences, prefix=prefix, args=args)
    files = {f.name: f for f in await results.list_files()}
    for suffix in ("fa", "fq"):
        if f"{prefix}.{suffix}" in files and f"{prefix}.{suffix}.clstr" in files:
            return CdhitResult(fasta=files[f"{prefix}.{suffix}"], clusters=files[f"{prefix}.{suffix}.clstr"])
    raise FileNotFoundError(f"cd-hit-est wrote no {prefix}.fa or {prefix}.fq with its clusters")
