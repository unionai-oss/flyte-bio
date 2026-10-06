"""qualimap — RNA-seq alignment QC (``qualimap rnaseq``).

Exposes :func:`qualimap_rnaseq`, which writes Qualimap's RNA-seq QC report
directory (``qualimapReport.html``, ``rnaseq_qc_results.txt``, raw data and
plots) for one BAM.

Qualimap records the BAM and GTF paths it was given in
``rnaseq_qc_results.txt``. As upstream's staging does, the inputs are
symlinked into a scratch directory under their (original) basenames and
passed as relative paths, so the report names them the same way.
"""

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

# Pinned biocontainer URI. Upstream's main.nf names the bare
# `biocontainers/qualimap:...` (Docker Hub) — use the quay.io mirror.
QUALIMAP_IMAGE = "quay.io/biocontainers/qualimap:2.3--hdfd78af_0"

QUALIMAP_MEMORY_MB = 8192
# As upstream: give the JVM 80% of the task's memory.
QUALIMAP_JAVA_MB = int(QUALIMAP_MEMORY_MB * 0.8)
DEFAULT_RESOURCES = flyte.Resources(cpu=2, memory=f"{QUALIMAP_MEMORY_MB}Mi")

# Upstream maps the sample's strandedness to Qualimap's protocol names.
qualimap_rnaseq_cmd = shell.create(
    name="qualimap_rnaseq",
    image=QUALIMAP_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"bam": File, "gtf": File, "prefix": str, "strandedness": str, "single_end": bool, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        BAM=({{inputs.bam}}); GTF=({{inputs.gtf}})
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        case {{inputs.strandedness}} in
            forward) STRAND=strand-specific-forward ;;
            reverse) STRAND=strand-specific-reverse ;;
            *) STRAND=non-strand-specific ;;
        esac
        PE=()
        [ {{inputs.single_end}} = true ] || PE=(-pe)
        W=$(mktemp -d); cd "$W"
        ln -s "${{BAM[0]}}" "$(basename "${{BAM[0]}}")"
        ln -s "${{GTF[0]}}" "$(basename "${{GTF[0]}}")"
        unset DISPLAY
        mkdir -p tmp
        export _JAVA_OPTIONS=-Djava.io.tmpdir=./tmp
        qualimap \
            --java-mem-size={QUALIMAP_JAVA_MB}M \
            rnaseq \
            $ARGS \
            -bam "$(basename "${{BAM[0]}}")" \
            -gtf "$(basename "${{GTF[0]}}")" \
            -p $STRAND \
            "${{PE[@]}}" \
            -outdir "$P"
        cp -r "$P"/. {{outputs.results}}/
    """,
)


env = flyte.TaskEnvironment.from_task(
    "qualimap",
    qualimap_rnaseq_cmd.as_task(),
)


async def qualimap_rnaseq(
    bam: File, gtf: File, prefix: str, strandedness: str = "unstranded", single_end: bool = False, args: str = ""
) -> Dir:
    """Run ``qualimap rnaseq`` on ``bam``; returns the report directory.

    ``strandedness`` is forward / reverse / anything else (non-strand-specific).
    Pass ``args="--sorted"`` for a name-sorted BAM, as upstream does.
    """
    return await qualimap_rnaseq_cmd(
        bam=bam, gtf=gtf, prefix=prefix, strandedness=strandedness, single_end=single_end, args=args
    )
