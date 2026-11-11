"""repeatmasker — mask interspersed repeats and low-complexity DNA with RepeatMasker.

Exposes :func:`repeatmasker`, which masks the repeats of a FASTA (may be
gzipped; decompressed first, as upstream) found with a repeat library
(``-lib``: a FASTA of consensus sequences, or with ``-engine hmmer`` in
``args`` an HMM library such as a Dfam export). Outputs keep upstream's names:

- ``<prefix>.masked``: the masked FASTA;
- ``<prefix>.out``: the repeats annotated;
- ``<prefix>.tbl``: a summary table (it records the run's date);
- ``<prefix>.gff``: with ``-gff`` in ``args``, the repeats as GFF.

``-pa`` follows the task's CPUs.

Libraries: the image carries RepeatMasker 4.1.5's bundled libraries only:
``RepeatMasker.lib`` (Dfam's curated consensus sequences, CC0),
``RepeatPeps.lib`` (repeat proteins, distributed with RepeatMasker) and
``RMRBMeta.embl`` (metadata); its ``Dfam.h5`` is a placeholder, so species
names (``-species``) don't resolve to the full Dfam. No RepBase library (which
is licensed) is included. Pass a library (``lib``) for real work: a Dfam
partition export (CC0), a RepeatModeler library, or RepBase under your own
licence.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

REPEATMASKER_IMAGE = "community.wave.seqera.io/library/repeatmasker:4.1.5--3cc5e72bc3d6ed07"
# Upstream's label is process_high. These fit test-sized genomes, and the 4
# CPUs match upstream's tests; raise them with repeatmasker_cmd.override(resources=...).
DEFAULT_RESOURCES = flyte.Resources(cpu=4, memory="16Gi")

repeatmasker_cmd = shell.create(
    name="repeatmasker",
    image=REPEATMASKER_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"fasta": File, "lib": File | None, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        FA=({{inputs.fasta}})
        shopt -s nullglob; LIB=({{inputs.lib}}); shopt -u nullglob
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"

        NAME=$(basename "${{FA[0]}}")
        if [[ $NAME == *.gz ]]; then
            NAME=${{NAME%.gz}}
            gunzip -c "${{FA[0]}}" > "$NAME" || exit 1
        else
            ln -s "${{FA[0]}}" "$NAME"
        fi
        LIB_ARGS=()
        if [ ${{#LIB[@]}} -gt 0 ]; then
            ln -s "${{LIB[0]}}" "$(basename "${{LIB[0]}}")"
            LIB_ARGS=(-lib "$(basename "${{LIB[0]}}")")
        fi

        RepeatMasker "${{LIB_ARGS[@]}}" -pa $TASK_CPUS -dir "$P" $ARGS "$NAME" || exit 1

        mv "$P/$NAME.masked" "$OUT/$P.masked"
        mv "$P/$NAME.out" "$OUT/$P.out"
        mv "$P/$NAME.tbl" "$OUT/$P.tbl"
        mv "$P/$NAME.out.gff" "$OUT/$P.gff" 2>/dev/null || echo "GFF is not produced"
    """,
)

env = flyte.TaskEnvironment.from_task("repeatmasker", repeatmasker_cmd.as_task())


@dataclass
class RepeatmaskerResult:
    masked: File  # <prefix>.masked
    out: File  # <prefix>.out
    tbl: File  # <prefix>.tbl
    gff: File | None  # <prefix>.gff, with -gff in args


async def repeatmasker(
    fasta: File, lib: File | None = None, prefix: str = "sample", args: str = ""
) -> RepeatmaskerResult:
    """Mask the repeats of ``fasta`` found with the library ``lib`` (or RepeatMasker's bundled one)."""
    results = await repeatmasker_cmd(fasta=fasta, lib=lib, prefix=prefix, args=args)

    async def pick(name: str) -> File:
        found = await results.get_file(name)
        if found is None:
            raise FileNotFoundError(f"RepeatMasker wrote no {name}")
        return found

    return RepeatmaskerResult(
        masked=await pick(f"{prefix}.masked"),
        out=await pick(f"{prefix}.out"),
        tbl=await pick(f"{prefix}.tbl"),
        gff=await results.get_file(f"{prefix}.gff"),
    )
