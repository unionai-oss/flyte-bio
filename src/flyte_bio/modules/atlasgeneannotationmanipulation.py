"""atlasgeneannotationmanipulation — feature annotation tables from a GTF.

Exposes :func:`atlasgeneannotationmanipulation_gtf2featureannotation`, which
runs ``gtf2featureAnnotation.R`` (Expression Atlas' gene annotation
manipulation) on a GTF, writing ``<prefix>.anno.tsv``: one row per feature
(genes by default; ``--feature-type transcript`` in ``args`` for
transcripts) with its GTF attributes. Given a cDNA FASTA, it is parsed
(``--parse-cdnas``) so ``args`` can ask for the cDNAs matching the
annotation, e.g. ``--filter-cdnas-output <prefix>.filtered.fa.gz``; any
``*.fa.gz`` written is returned as ``filtered_cdna``.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

ATLAS_IMAGE = "community.wave.seqera.io/library/atlas-gene-annotation-manipulation:1.1.1--b9e8c32c92709512"
# Upstream: process_single.
DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory="6Gi")

gtf2featureannotation_cmd = shell.create(
    name="atlasgeneannotationmanipulation_gtf2featureannotation",
    image=ATLAS_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"gtf": File, "fasta": File | None, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        GTF=({inputs.gtf})
        shopt -s nullglob; FASTA=({inputs.fasta}); shopt -u nullglob
        P={inputs.prefix}
        ARGS={inputs.args}
        OUT={outputs.results}
        W=$(mktemp -d); cd "$W"
        GTF_NAME=$(basename "${GTF[0]}"); ln -s "${GTF[0]}" "$GTF_NAME"
        CDNA=()
        if [ ${#FASTA[@]} -gt 0 ]; then
            FASTA_NAME=$(basename "${FASTA[0]}"); ln -s "${FASTA[0]}" "$FASTA_NAME"
            CDNA=(--parse-cdnas "$FASTA_NAME")
        fi
        gtf2featureAnnotation.R --gtf-file "$GTF_NAME" --output-file "$P.anno.tsv" "${CDNA[@]}" $ARGS || exit 1
        shopt -s nullglob
        FILTERED=(./*.fa.gz)
        for f in "${FILTERED[@]}"; do [ -L "$f" ] || mv "$f" "$OUT"/; done
        mv "$P.anno.tsv" "$OUT"/
    """,
)

env = flyte.TaskEnvironment.from_task("atlasgeneannotationmanipulation", gtf2featureannotation_cmd.as_task())


@dataclass
class FeatureAnnotationResult:
    feature_annotation: File  # <prefix>.anno.tsv
    filtered_cdna: File | None  # the *.fa.gz that args ask for


async def atlasgeneannotationmanipulation_gtf2featureannotation(
    gtf: File, prefix: str, fasta: File | None = None, args: str = ""
) -> FeatureAnnotationResult:
    """``<prefix>.anno.tsv`` from ``gtf``, parsing the cDNA ``fasta`` when given."""
    results = await gtf2featureannotation_cmd(gtf=gtf, fasta=fasta, prefix=prefix, args=args)
    annotation = await results.get_file(f"{prefix}.anno.tsv")
    if annotation is None:
        raise FileNotFoundError(f"gtf2featureAnnotation.R wrote no {prefix}.anno.tsv")
    cdnas = [f for f in await results.list_files() if f.name.endswith(".fa.gz")]
    return FeatureAnnotationResult(feature_annotation=annotation, filtered_cdna=cdnas[0] if cdnas else None)
