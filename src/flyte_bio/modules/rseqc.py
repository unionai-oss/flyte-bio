"""rseqc — RNA-seq BAM quality control (RSeQC 5).

Exposes one task per RSeQC script, as upstream's modules: :data:`bam_stat`,
:data:`infer_experiment`, :data:`inner_distance`,
:data:`junction_annotation`, :data:`junction_saturation`,
:data:`read_distribution` and :data:`read_duplication`. Each writes its
outputs under upstream's ``<prefix>.<suffix>`` names into an output Dir.

As upstream's staging does, every task lays its inputs out in a scratch
directory under their original basenames — the BAM with its index beside
it as ``<bam>.bai`` (where pysam looks for it) and the BED — runs the
script there with relative paths, and then moves out everything except the
inputs. RSeQC's logs and generated R scripts therefore name files exactly as
upstream's do.
"""

import asyncio

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

RSEQC_IMAGE = "community.wave.seqera.io/library/rseqc_r-base:2e29d2dfda9cef15"

DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory="8Gi")

# Stage BAM (+ .bai) and, when given, the BED under their basenames in a
# scratch dir; $B / $BED name them for the tool, $P is the output prefix.
STAGE = r"""
        BAM=({inputs.bam}); BAI=({inputs.bai}); P={inputs.prefix}
        W=$(mktemp -d); cd "$W"
        B=$(basename "${BAM[0]}")
        ln -s "${BAM[0]}" "$B"; ln -s "${BAI[0]}" "$B.bai"
        STAGED=("$B" "$B.bai")
"""
STAGE_BED = r"""
        BEDS=({inputs.bed}); BED=$(basename "${BEDS[0]}")
        ln -s "${BEDS[0]}" "$BED"; STAGED+=("$BED")
"""
# Move every output (anything that isn't a staged input) to the results Dir.
COLLECT = r"""
        for f in *; do
            keep=1; for s in "${STAGED[@]}"; do [ "$f" = "$s" ] && keep=0; done
            [ $keep = 1 ] && mv "$f" {outputs.results}/
        done
        true  # the loop's last test can be false (a staged input); that's not a failure
"""

BAM_INPUTS = {"bam": File, "bai": File, "prefix": str}
BED_INPUTS = {**BAM_INPUTS, "bed": File}


def rseqc_task(name: str, inputs: dict, body: str):
    return shell.create(
        name=name,
        image=RSEQC_IMAGE,
        resources=DEFAULT_RESOURCES,
        inputs=inputs,
        outputs={"results": Dir},
        script=STAGE + (STAGE_BED if "bed" in inputs else "") + body + COLLECT,
    )


bam_stat = rseqc_task(
    "rseqc_bam_stat",
    BAM_INPUTS,
    r"""
        bam_stat.py -i "$B" > "$P.bam_stat.txt"
""",
)

infer_experiment = rseqc_task(
    "rseqc_infer_experiment",
    BED_INPUTS,
    r"""
        infer_experiment.py -i "$B" -r "$BED" > "$P.infer_experiment.txt"
""",
)

# As upstream: inner distance is only defined for paired-end data.
inner_distance = rseqc_task(
    "rseqc_inner_distance",
    {**BED_INPUTS, "single_end": bool},
    r"""
        if [ {inputs.single_end} = true ]; then
            echo "inner_distance.py doesn't support single-end data" > "$P.inner_distance.txt"
        else
            inner_distance.py -i "$B" -r "$BED" -o "$P" > stdout.txt
            head -n 2 stdout.txt > "$P.inner_distance_mean.txt"
            STAGED+=(stdout.txt)
        fi
""",
)

junction_annotation = rseqc_task(
    "rseqc_junction_annotation",
    BED_INPUTS,
    r"""
        junction_annotation.py -i "$B" -r "$BED" -o "$P" \
            2>| >(grep -v 'E::idx_find_and_load' | tee "$P.junction_annotation.log" >&2)
        sleep 1  # let the tee'd log finish writing
""",
)

junction_saturation = rseqc_task(
    "rseqc_junction_saturation",
    BED_INPUTS,
    r"""
        junction_saturation.py -i "$B" -r "$BED" -o "$P"
""",
)

read_distribution = rseqc_task(
    "rseqc_read_distribution",
    BED_INPUTS,
    r"""
        read_distribution.py -i "$B" -r "$BED" > "$P.read_distribution.txt"
""",
)

read_duplication = rseqc_task(
    "rseqc_read_duplication",
    BAM_INPUTS,
    r"""
        read_duplication.py -i "$B" -o "$P"
""",
)


env = flyte.TaskEnvironment.from_task(
    "rseqc",
    bam_stat.as_task(),
    infer_experiment.as_task(),
    inner_distance.as_task(),
    junction_annotation.as_task(),
    junction_saturation.as_task(),
    read_distribution.as_task(),
    read_duplication.as_task(),
)


# Upstream's default `rseqc_modules`.
DEFAULT_MODULES = (
    "bam_stat",
    "inner_distance",
    "infer_experiment",
    "junction_annotation",
    "junction_saturation",
    "read_distribution",
    "read_duplication",
)


async def rseqc(
    bam: File,
    bai: File,
    bed: File,
    prefix: str,
    single_end: bool,
    modules: tuple[str, ...] = DEFAULT_MODULES,
) -> dict[str, Dir]:
    """Run the requested RSeQC ``modules`` in parallel; returns module name → output Dir."""
    unknown = set(modules) - set(DEFAULT_MODULES)
    if unknown:
        raise ValueError(f"unknown RSeQC modules: {sorted(unknown)}")
    calls = {
        "bam_stat": lambda: bam_stat(bam=bam, bai=bai, prefix=prefix),
        "infer_experiment": lambda: infer_experiment(bam=bam, bai=bai, bed=bed, prefix=prefix),
        "inner_distance": lambda: inner_distance(bam=bam, bai=bai, bed=bed, prefix=prefix, single_end=single_end),
        "junction_annotation": lambda: junction_annotation(bam=bam, bai=bai, bed=bed, prefix=prefix),
        "junction_saturation": lambda: junction_saturation(bam=bam, bai=bai, bed=bed, prefix=prefix),
        "read_distribution": lambda: read_distribution(bam=bam, bai=bai, bed=bed, prefix=prefix),
        "read_duplication": lambda: read_duplication(bam=bam, bai=bai, prefix=prefix),
    }
    results = await asyncio.gather(*(calls[m]() for m in modules))
    return dict(zip(modules, results))
