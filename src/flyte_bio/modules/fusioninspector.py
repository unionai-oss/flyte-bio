"""fusioninspector — in-depth validation of candidate gene fusions with FusionInspector.

Exposes :func:`fusioninspector`, which realigns one sample's reads to
fusion-gene contigs built from a CTAT genome library for a list of candidate
fusions (``GENE1--GENE2``, one per line, or a STAR-Fusion / fusion-report
table), and scores the evidence. With ``--vis`` (always on, as upstream) it
also writes an IGV web report. Outputs, under upstream's names:

- ``<prefix>.FusionInspector.fusions.tsv`` and
  ``<prefix>.FusionInspector.fusions.abridged.tsv``;
- ``FusionInspector.log``;
- the HTML report(s) (``*html``);
- the fusion contigs' annotation, ``fi_workdir/*.gtf``, as ``<name>.gtf``;
- the ``IGV_inputs``, ``fi_workdir`` and ``chckpts_dir`` directories, as
  ``<dir>.tar.gz`` (unlike upstream, which emits the directories: they hold
  many files, which Flyte stages poorly).

The CTAT genome library travels as a ``.tar.gz`` of its directory, as in
:mod:`flyte_bio.modules.starfusion`. Resources default to what upstream's
test runs with (4 CPUs); a human library needs more memory (override
``resources``).
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS
from flyte_bio.modules.starfusion import ARCHIVE_DIR, EXTRACT_GENOME_LIB

FUSIONINSPECTOR_IMAGE = (
    "community.wave.seqera.io/library/fusion-inspector_perl-json-xs_perl-carp-assert_pip_pruned:367be466d24aba4a"
)
DEFAULT_RESOURCES = flyte.Resources(cpu=4, memory="16Gi")

fusioninspector_cmd = shell.create(
    name="fusioninspector",
    image=FUSIONINSPECTOR_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "reads_1": File,
        "reads_2": File | None,
        "fusion_list": File,
        "genome_lib": File,
        "prefix": str,
        "args": str,
    },
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        {ARCHIVE_DIR}
        R1=({{inputs.reads_1}}); FUSIONS=({{inputs.fusion_list}}); LIB=({{inputs.genome_lib}})
        shopt -s nullglob; R2=({{inputs.reads_2}}); shopt -u nullglob
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        {EXTRACT_GENOME_LIB}
        W=$(mktemp -d); cd "$W"
        READS=(--left_fq "${{R1[0]}}")
        if [ ${{#R2[@]}} -gt 0 ]; then READS+=(--right_fq "${{R2[0]}}"); fi
        FusionInspector \
            --fusions "${{FUSIONS[0]}}" \
            --genome_lib "$GENOME_LIB" \
            "${{READS[@]}}" \
            --CPU $TASK_CPUS \
            -O . \
            --out_prefix "$P" \
            --vis $ARGS || exit 1
        shopt -s nullglob
        mv ./*FusionInspector.fusions.tsv ./*FusionInspector.log ./*html ./*abridged.tsv "$OUT"/
        for f in fi_workdir/*.gtf; do cp "$f" "$OUT"/; done
        for d in IGV_inputs fi_workdir chckpts_dir; do
            if [ -d "$d" ]; then archive_dir "$OUT/$d.tar.gz" "$d"; fi
        done
    """,
)

env = flyte.TaskEnvironment.from_task("fusioninspector", fusioninspector_cmd.as_task())


@dataclass
class FusionInspectorResult:
    tsv: File | None  # <prefix>.FusionInspector.fusions.tsv
    abridged_tsv: File | None  # <prefix>.FusionInspector.fusions.abridged.tsv
    log: File | None  # *FusionInspector.log
    gtfs: list[File]  # fi_workdir/*.gtf: the fusion contigs' annotation
    html: list[File]  # the IGV web report(s)
    results: Dir  # everything, incl. IGV_inputs.tar.gz, fi_workdir.tar.gz and chckpts_dir.tar.gz


async def fusioninspector(
    reads_1: File,
    fusion_list: File,
    genome_lib: File,
    prefix: str,
    reads_2: File | None = None,
    args: str = "",
) -> FusionInspectorResult:
    """Inspect the candidate fusions in ``fusion_list`` in one sample's reads (``reads_2`` None for single-end).

    ``genome_lib`` is a CTAT genome library archive (``.tar.gz``). Options such as ``--annotate`` and
    ``--examine_coding_effect`` go in ``args``.
    """
    results = await fusioninspector_cmd(
        reads_1=reads_1,
        reads_2=reads_2,
        fusion_list=fusion_list,
        genome_lib=genome_lib,
        prefix=prefix,
        args=args,
    )
    files = sorted([f async for f in results.walk()], key=lambda f: f.name)
    return FusionInspectorResult(
        tsv=await results.get_file(f"{prefix}.FusionInspector.fusions.tsv"),
        abridged_tsv=await results.get_file(f"{prefix}.FusionInspector.fusions.abridged.tsv"),
        log=next((f for f in files if f.name.endswith("FusionInspector.log")), None),
        gtfs=[f for f in files if f.name.endswith(".gtf")],
        html=[f for f in files if f.name.endswith("html")],
        results=results,
    )
