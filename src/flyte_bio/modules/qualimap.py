"""qualimap — alignment QC (``qualimap rnaseq``, ``qualimap bamqc``).

Exposes:

- :func:`qualimap_rnaseq`, which writes Qualimap's RNA-seq QC report
  directory (``qualimapReport.html``, ``rnaseq_qc_results.txt``, raw data and
  plots) for one BAM.
- :func:`qualimap_bamqc`, which writes Qualimap's BAM QC report directory
  (``qualimapReport.html``, ``genome_results.txt``, ``css/``,
  ``images_qualimapReport/``, ``raw_data_qualimapReport/``) for one BAM,
  optionally restricted to the regions of a GFF/BED (``gff``). Paired-end
  samples add ``--collect-overlap-pairs``, as upstream; threads follow the
  task's CPUs. Pass the protocol in ``args`` (``-p non-strand-specific``, ...),
  as upstream's ``ext.args``.
- :func:`qualimap_bamqccram`, which writes the same BAM QC report directory
  for a CRAM, decoded with ``samtools view -T <fasta>`` and streamed to
  ``qualimap bamqc -bam /dev/stdin`` (so the report names the BAM
  ``/dev/stdin``), as upstream.

Qualimap records the BAM and GTF paths it was given in
``rnaseq_qc_results.txt``. As upstream's staging does, the inputs are
symlinked into a scratch directory under their (original) basenames and
passed as relative paths, so the report names them the same way.
"""

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import JAVA_HEAP_MB, TASK_LIMITS

# Pinned biocontainer URI. Upstream's main.nf names the bare
# `biocontainers/qualimap:...` (Docker Hub) — use the quay.io mirror.
QUALIMAP_IMAGE = "quay.io/biocontainers/qualimap:2.3--hdfd78af_0"
# bamqc is pinned upstream to a newer build of the same release.
QUALIMAP_BAMQC_IMAGE = "community.wave.seqera.io/library/qualimap:2.3--e70af3f76915782f"
# bamqccram's image adds samtools, to decode the CRAM.
QUALIMAP_SAMTOOLS_IMAGE = "community.wave.seqera.io/library/qualimap_samtools:36c0048b499c176d"

QUALIMAP_MEMORY_MB = 8192
DEFAULT_RESOURCES = flyte.Resources(cpu=2, memory=f"{QUALIMAP_MEMORY_MB}Mi")
# Upstream labels bamqc process_medium.
BAMQC_RESOURCES = flyte.Resources(cpu=4, memory="16Gi")

# Upstream maps the sample's strandedness to Qualimap's protocol names.
qualimap_rnaseq_cmd = shell.create(
    name="qualimap_rnaseq",
    image=QUALIMAP_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "bam": File,
        "gtf": File,
        "prefix": str,
        "strandedness": str,
        "single_end": bool,
        "bam_name": str,
        "args": str,
    },
    defaults={"bam_name": "", "args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
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
        BAM_NAME={{inputs.bam_name}}
        [ -n "$BAM_NAME" ] || BAM_NAME=$(basename "${{BAM[0]}}")
        ln -s "${{BAM[0]}}" "$BAM_NAME"
        ln -s "${{GTF[0]}}" "$(basename "${{GTF[0]}}")"
        unset DISPLAY
        mkdir -p tmp
        export _JAVA_OPTIONS=-Djava.io.tmpdir=./tmp
        qualimap \
            --java-mem-size={JAVA_HEAP_MB}M \
            rnaseq \
            $ARGS \
            -bam "$BAM_NAME" \
            -gtf "$(basename "${{GTF[0]}}")" \
            -p $STRAND \
            "${{PE[@]}}" \
            -outdir "$P"
        cp -r "$P"/. {{outputs.results}}/
    """,
)


# As for rnaseq, the BAM (and GFF) are linked into the working directory under
# their own names, which Qualimap records in genome_results.txt.
qualimap_bamqc_cmd = shell.create(
    name="qualimap_bamqc",
    image=QUALIMAP_BAMQC_IMAGE,
    resources=BAMQC_RESOURCES,
    inputs={"bam": File, "gff": File | None, "prefix": str, "single_end": bool, "args": str},
    defaults={"single_end": False, "args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        BAM=({{inputs.bam}})
        shopt -s nullglob; GFF=({{inputs.gff}}); shopt -u nullglob
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        W=$(mktemp -d); cd "$W"
        BAM_NAME=$(basename "${{BAM[0]}}")
        ln -s "${{BAM[0]}}" "$BAM_NAME"
        REGIONS=()
        if [ ${{#GFF[@]}} -gt 0 ]; then
            ln -s "${{GFF[0]}}" "$(basename "${{GFF[0]}}")"
            REGIONS=(--gff "$(basename "${{GFF[0]}}")")
        fi
        PAIRS=()
        [ {{inputs.single_end}} = true ] || PAIRS=(--collect-overlap-pairs)
        unset DISPLAY
        mkdir -p tmp
        export _JAVA_OPTIONS=-Djava.io.tmpdir=./tmp
        qualimap \
            --java-mem-size={JAVA_HEAP_MB}M \
            bamqc \
            $ARGS \
            -bam "$BAM_NAME" \
            "${{REGIONS[@]}}" \
            "${{PAIRS[@]}}" \
            -outdir "$P" \
            -nt $TASK_CPUS
        cp -r "$P"/. {{outputs.results}}/
    """,
)

# The CRAM is decoded against the FASTA (linked next to its .fai) and streamed
# to Qualimap; upstream's label is process_medium, as for bamqc.
qualimap_bamqccram_cmd = shell.create(
    name="qualimap_bamqccram",
    image=QUALIMAP_SAMTOOLS_IMAGE,
    resources=BAMQC_RESOURCES,
    inputs={
        "cram": File,
        "crai": File,
        "gff": File | None,
        "fasta": File,
        "fai": File,
        "prefix": str,
        "single_end": bool,
        "args": str,
    },
    defaults={"single_end": False, "args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        set -o pipefail
        CRAM=({{inputs.cram}}); CRAI=({{inputs.crai}}); FA=({{inputs.fasta}}); FAI=({{inputs.fai}})
        shopt -s nullglob; GFF=({{inputs.gff}}); shopt -u nullglob
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        W=$(mktemp -d); cd "$W"
        for F in "${{CRAM[0]}}" "${{CRAI[0]}}" "${{FA[0]}}" "${{FAI[0]}}"; do ln -s "$F" "$(basename "$F")"; done
        REGIONS=()
        if [ ${{#GFF[@]}} -gt 0 ]; then
            ln -s "${{GFF[0]}}" "$(basename "${{GFF[0]}}")"
            REGIONS=(--gff "$(basename "${{GFF[0]}}")")
        fi
        PAIRS=()
        [ {{inputs.single_end}} = true ] || PAIRS=(--collect-overlap-pairs)
        unset DISPLAY
        mkdir -p tmp
        export _JAVA_OPTIONS=-Djava.io.tmpdir=./tmp
        samtools view -hb -T "$(basename "${{FA[0]}}")" "$(basename "${{CRAM[0]}}")" |
        qualimap \
            --java-mem-size={JAVA_HEAP_MB}M \
            bamqc \
            $ARGS \
            -bam /dev/stdin \
            "${{REGIONS[@]}}" \
            "${{PAIRS[@]}}" \
            -outdir "$P" \
            -nt $TASK_CPUS
        cp -r "$P"/. {{outputs.results}}/
    """,
)

# One image per environment: bamqc's and bamqccram's differ from rnaseq's, so `env` has three.
qualimap_rnaseq_env = flyte.TaskEnvironment.from_task("qualimap_rnaseq", qualimap_rnaseq_cmd.as_task())
qualimap_bamqc_env = flyte.TaskEnvironment.from_task("qualimap_bamqc", qualimap_bamqc_cmd.as_task())
qualimap_bamqccram_env = flyte.TaskEnvironment.from_task("qualimap_bamqccram", qualimap_bamqccram_cmd.as_task())
env = flyte.TaskEnvironment(
    name="qualimap", depends_on=[qualimap_rnaseq_env, qualimap_bamqc_env, qualimap_bamqccram_env]
)


async def qualimap_rnaseq(
    bam: File,
    gtf: File,
    prefix: str,
    strandedness: str = "unstranded",
    single_end: bool = False,
    args: str = "",
    bam_name: str = "",
) -> Dir:
    """Run ``qualimap rnaseq`` on ``bam``; returns the report directory.

    ``strandedness`` is forward / reverse / anything else (non-strand-specific).
    Pass ``args="--sorted"`` for a name-sorted BAM, as upstream does.
    ``bam_name`` stages the BAM under that name (Qualimap records it, and
    MultiQC takes the sample name from it); default: its own basename.
    """
    return await qualimap_rnaseq_cmd(
        bam=bam,
        gtf=gtf,
        prefix=prefix,
        strandedness=strandedness,
        single_end=single_end,
        bam_name=bam_name,
        args=args,
    )


async def qualimap_bamqc(
    bam: File, prefix: str, gff: File | None = None, single_end: bool = False, args: str = ""
) -> Dir:
    """Run ``qualimap bamqc`` on ``bam`` (within ``gff``'s regions, if given); returns the report directory."""
    return await qualimap_bamqc_cmd(bam=bam, gff=gff, prefix=prefix, single_end=single_end, args=args)


async def qualimap_bamqccram(
    cram: File,
    crai: File,
    fasta: File,
    fai: File,
    prefix: str,
    gff: File | None = None,
    single_end: bool = False,
    args: str = "",
) -> Dir:
    """Run ``qualimap bamqc`` on ``cram``, decoded against ``fasta`` (within ``gff``'s regions, if given)."""
    return await qualimap_bamqccram_cmd(
        cram=cram, crai=crai, gff=gff, fasta=fasta, fai=fai, prefix=prefix, single_end=single_end, args=args
    )
