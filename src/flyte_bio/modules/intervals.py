"""intervals — genomic intervals for scatter-gather.

Exposes:

- :func:`build_intervals` — one BED interval per sequence of a FASTA index
  (``.fai``), covering it whole, as upstream builds intervals when none are
  given.
- :func:`create_intervals_bed` — split intervals into BED files for
  scatter-gather: one file per interval, except that short ones are grouped
  so each file's estimated runtime (length / ``nucleotides_per_second``, or a
  BED's fifth column) is comparable, as upstream does. Takes a BED, a Picard
  ``.interval_list`` or ``chr:start-end`` lines.

Both are upstream's awk programs, unchanged.
"""

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

GAWK_IMAGE = "community.wave.seqera.io/library/gawk:5.3.1--e09efb5dfc4b8156"
# Upstream's create_intervals_bed names the bare `biocontainers/gawk:5.3.0`
# (Docker Hub); use the quay.io mirror.
CREATE_INTERVALS_IMAGE = "quay.io/biocontainers/gawk:5.3.0"

DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory="2Gi")

build_intervals_cmd = shell.create(
    name="build_intervals",
    image=GAWK_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"fai": File, "prefix": str},
    outputs={"results": Dir},
    script=r"""
        FAI=({inputs.fai})
        awk -v FS='\t' -v OFS='\t' '{ print $1, "0", $2 }' "${FAI[0]}" > {outputs.results}/{inputs.prefix}.bed
    """,
)

create_intervals_bed_cmd = shell.create(
    name="create_intervals_bed",
    image=CREATE_INTERVALS_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"intervals": File, "nucleotides_per_second": int},
    outputs={"results": Dir},
    script=r"""
        IN=({inputs.intervals})
        NPS={inputs.nucleotides_per_second}
        NAME=$(basename "${IN[0]}" | tr '[:upper:]' '[:lower:]')
        cd {outputs.results}
        case "$NAME" in
        *bed)
            awk -vFS="\t" -v nps="$NPS" '{
                t = $5  # runtime estimate
                if (t == "") {
                    # no runtime estimate in this row, assume default value
                    t = ($3 - $2) / nps
                }
                if (name == "" || (chunk > 600 && (chunk + t) > longest * 1.05)) {
                    # start a new chunk
                    name = sprintf("%s_%d-%d.bed", $1, $2+1, $3)
                    chunk = 0
                    longest = 0
                }
                if (t > longest)
                    longest = t
                chunk += t
                print $0 > name
            }' "${IN[0]}"
            ;;
        *interval_list)
            grep -v '^@' "${IN[0]}" | awk -vFS="\t" '{
                name = sprintf("%s_%d-%d", $1, $2, $3);
                printf("%s\t%d\t%d\n", $1, $2-1, $3) > name ".bed"
            }'
            ;;
        *)
            awk -vFS="[:-]" '{
                name = sprintf("%s_%d-%d", $1, $2, $3);
                printf("%s\t%d\t%d\n", $1, $2-1, $3) > name ".bed"
            }' "${IN[0]}"
            ;;
        esac
    """,
)


# One image per tool, so one env each; `env` aggregates them.
build_intervals_env = flyte.TaskEnvironment.from_task("intervals_build", build_intervals_cmd.as_task())
create_intervals_bed_env = flyte.TaskEnvironment.from_task("intervals_create_bed", create_intervals_bed_cmd.as_task())
env = flyte.TaskEnvironment(name="intervals", depends_on=[build_intervals_env, create_intervals_bed_env])


async def build_intervals(fai: File, prefix: str) -> File:
    """``<prefix>.bed``: each sequence in ``fai`` as one whole interval."""
    results = await build_intervals_cmd(fai=fai, prefix=prefix)
    bed = await results.get_file(f"{prefix}.bed")
    if bed is None:
        raise FileNotFoundError(f"build_intervals wrote no {prefix}.bed")
    return bed


async def create_intervals_bed(intervals: File, nucleotides_per_second: int = 200000) -> list[File]:
    """Split ``intervals`` into scatter-gather BED files, sorted longest-running first.

    Each file's estimated runtime is the sum of its intervals' (a BED's fifth
    column, or length / ``nucleotides_per_second``). Upstream processes the
    biggest chunks first, so they come first here too.
    """
    results = await create_intervals_bed_cmd(intervals=intervals, nucleotides_per_second=nucleotides_per_second)
    beds = [f async for f in results.walk() if f.path.endswith(".bed")]

    async def runtime(bed: File) -> float:
        async with bed.open("rb") as fh:
            text = bytes(await fh.read()).decode()
        total = 0.0
        for line in text.splitlines():
            fields = line.split("\t")
            if len(fields) >= 5:
                total += float(fields[4])
            elif len(fields) >= 3:
                total += (int(fields[2]) - int(fields[1])) / nucleotides_per_second
        return total

    runtimes = [await runtime(b) for b in beds]
    return [b for _, b in sorted(zip(runtimes, beds), key=lambda pair: -pair[0])]
