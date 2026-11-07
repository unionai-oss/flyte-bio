"""spaceranger — 10x Genomics Space Ranger: spatial gene expression (Visium, Visium HD).

Exposes:

- :func:`spaceranger_mkref`, which builds a Space Ranger reference from a
  genome FASTA and a GTF (``spaceranger mkref``), and returns it as
  ``<reference_name>.tar.gz``, an archive of the ``<reference_name>`` folder
  (a STAR index plus gene annotations). It's an archive rather than a Dir as
  a reference holds too many files for a Dir input to stage reliably.
- :func:`spaceranger_count`, which counts one Visium or Visium HD sample
  (``spaceranger count``) against such a reference, extracting it first, and
  returns Space Ranger's ``outs`` folder: ``web_summary.html``,
  ``metrics_summary.csv``, the filtered and raw feature-barcode matrices,
  ``spatial/`` (images, tissue positions, scale factors), ``analysis/``
  and, for Visium HD, ``binned_outputs/`` and ``segmented_outputs/``.

As upstream, :func:`spaceranger_count`:

- takes the sample's FASTQs (``reads``, named as ``bcl2fastq`` names them,
  ``<sample>_S1_L001_R1_001.fastq.gz``) and its ``sample`` name;
- passes each image given (``image``, ``cytaimage``, ``darkimage``,
  ``colorizedimage``), a manual alignment (``alignment``), a slide layout
  file (``slidefile``) and a probe set (``probeset``) only when given;
- passes ``--slide`` and ``--area``, except for a ``visium-*`` slide with no
  area and no slide file, which is passed as ``--unknown-slide``.

Space Ranger 4 requires ``--create-bam`` (``true`` or ``false``), so ``args``
must set it. Its cores and memory (``--localcores``, ``--localmem``, and
mkref's ``--nthreads``) are the task's. Telemetry is disabled
(``TENX_DISABLE_TELEMETRY``), as upstream's tests do.

The container is upstream's build of 10x's Space Ranger, distributed under
10x Genomics' end-user license.
"""

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

SPACERANGER_IMAGE = "quay.io/nf-core/spaceranger:9c5e7dc93c32448e"
NO_TELEMETRY = {"TENX_DISABLE_TELEMETRY": "1"}
# Upstream's process_high asks for 12 CPUs and 72 GB. Its tests' tiny genome
# and reads run in 4 CPUs and 7 GB; real genomes and samples want 32-64 GiB
# (10x recommends 64 GB or more): raise with
# <cmd>.as_task().override(resources=...).
MKREF_RESOURCES = flyte.Resources(cpu=4, memory="16Gi", disk="50Gi")
COUNT_RESOURCES = flyte.Resources(cpu=4, memory="16Gi", disk="50Gi")

spaceranger_mkref_cmd = shell.create(
    name="spaceranger_mkref",
    image=SPACERANGER_IMAGE,
    resources=MKREF_RESOURCES,
    env_vars=NO_TELEMETRY,
    inputs={"fasta": File, "gtf": File, "reference_name": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        FA=({{inputs.fasta}})
        GTF=({{inputs.gtf}})
        NAME={{inputs.reference_name}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        export HOME="$W"
        spaceranger \
            mkref \
            --genome="$NAME" \
            --fasta="${{FA[0]}}" \
            --genes="${{GTF[0]}}" \
            --localcores=$TASK_CPUS \
            --localmem=$(( TASK_MEMORY_MB / 1024 )) \
            --nthreads=$TASK_CPUS \
            $ARGS || exit 1
        tar -czf "$OUT/$NAME.tar.gz" "$NAME"
    """,
)

spaceranger_count_cmd = shell.create(
    name="spaceranger_count",
    image=SPACERANGER_IMAGE,
    resources=COUNT_RESOURCES,
    env_vars=NO_TELEMETRY,
    inputs={
        "reads": list[File],
        "reference": File,
        "sample": str,
        "prefix": str,
        "slide": str,
        "area": str,
        "image": File | None,
        "cytaimage": File | None,
        "darkimage": File | None,
        "colorizedimage": File | None,
        "alignment": File | None,
        "slidefile": File | None,
        "probeset": File | None,
        "args": str,
    },
    defaults={"area": "", "args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        READS=({{inputs.reads}})
        REF_ARCHIVE=({{inputs.reference}})
        shopt -s nullglob
        IMAGE=({{inputs.image}}); CYTA=({{inputs.cytaimage}}); DARK=({{inputs.darkimage}})
        COLORIZED=({{inputs.colorizedimage}}); ALIGNMENT=({{inputs.alignment}})
        SLIDEFILE=({{inputs.slidefile}}); PROBESET=({{inputs.probeset}})
        shopt -u nullglob
        SAMPLE={{inputs.sample}}
        P={{inputs.prefix}}
        SLIDE={{inputs.slide}}
        AREA={{inputs.area}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        export HOME="$W"

        mkdir reference fastqs
        tar -xzf "${{REF_ARCHIVE[0]}}" -C reference
        REF=(reference/*/)
        for R in "${{READS[@]}}"; do ln -s "$R" fastqs/; done

        EXTRA=()
        if [ ${{#IMAGE[@]}} -gt 0 ]; then EXTRA+=(--image="${{IMAGE[0]}}"); fi
        if [ ${{#CYTA[@]}} -gt 0 ]; then EXTRA+=(--cytaimage="${{CYTA[0]}}"); fi
        if [ ${{#DARK[@]}} -gt 0 ]; then EXTRA+=(--darkimage="${{DARK[0]}}"); fi
        if [ ${{#COLORIZED[@]}} -gt 0 ]; then EXTRA+=(--colorizedimage="${{COLORIZED[0]}}"); fi
        if [[ "$SLIDE" == visium-* ]] && [ -z "$AREA" ] && [ ${{#SLIDEFILE[@]}} -eq 0 ]; then
            EXTRA+=(--unknown-slide="$SLIDE")
        else
            EXTRA+=(--slide="$SLIDE" --area="$AREA")
        fi
        if [ ${{#PROBESET[@]}} -gt 0 ]; then EXTRA+=(--probe-set="${{PROBESET[0]}}"); fi
        if [ ${{#ALIGNMENT[@]}} -gt 0 ]; then EXTRA+=(--loupe-alignment="${{ALIGNMENT[0]}}"); fi
        if [ ${{#SLIDEFILE[@]}} -gt 0 ]; then EXTRA+=(--slidefile="${{SLIDEFILE[0]}}"); fi

        spaceranger count \
            --id="$P" \
            --sample="$SAMPLE" \
            --fastqs="$W/fastqs" \
            --transcriptome="$W/${{REF[0]%/}}" \
            --localcores=$TASK_CPUS \
            --localmem=$(( TASK_MEMORY_MB / 1024 )) \
            "${{EXTRA[@]}}" \
            $ARGS || exit 1
        mv "$P"/outs/* "$OUT"/
    """,
)

env = flyte.TaskEnvironment.from_task("spaceranger", spaceranger_mkref_cmd.as_task(), spaceranger_count_cmd.as_task())


async def spaceranger_mkref(fasta: File, gtf: File, reference_name: str, args: str = "") -> File:
    """A Space Ranger reference built from ``fasta`` and ``gtf``: ``<reference_name>.tar.gz``."""
    results = await spaceranger_mkref_cmd(fasta=fasta, gtf=gtf, reference_name=reference_name, args=args)
    archive = await results.get_file(f"{reference_name}.tar.gz")
    if archive is None:
        raise FileNotFoundError(f"spaceranger mkref wrote no {reference_name}.tar.gz")
    return archive


async def spaceranger_count(
    reads: list[File],
    reference: File,
    sample: str,
    slide: str,
    area: str = "",
    image: File | None = None,
    cytaimage: File | None = None,
    darkimage: File | None = None,
    colorizedimage: File | None = None,
    alignment: File | None = None,
    slidefile: File | None = None,
    probeset: File | None = None,
    prefix: str | None = None,
    args: str = "",
) -> Dir:
    """Space Ranger's ``outs`` for ``sample``'s ``reads`` against ``reference``, from :func:`spaceranger_mkref`."""
    if not reads:
        raise ValueError("spaceranger_count needs the sample's FASTQs")
    return await spaceranger_count_cmd(
        reads=reads,
        reference=reference,
        sample=sample,
        prefix=prefix or sample,
        slide=slide,
        area=area,
        image=image,
        cytaimage=cytaimage,
        darkimage=darkimage,
        colorizedimage=colorizedimage,
        alignment=alignment,
        slidefile=slidefile,
        probeset=probeset,
        args=args,
    )
