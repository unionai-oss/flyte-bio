"""annotsv — AnnotSV annotation and ranking of structural variants.

Exposes:

- :func:`annotsv_installannotations` — download AnnotSV's human annotation
  bundle (and, as the tool's install script intends, Exomiser's phenotype
  data) with ``INSTALL_annotations.sh``: ``AnnotSV_annotations.tar.gz``, an
  archive of the ``AnnotSV_annotations/`` directory.
- :func:`annotsv_annotsv` — annotate the structural variants of a VCF (or
  BED) against an annotation archive, optionally with candidate genes,
  candidate and false-positive small variants (for compound-heterozygosity
  and breakpoint filtering) and preferred gene transcripts:
  ``<prefix>.tsv`` (full and split annotation lines), and when written
  ``<prefix>.unannotated.tsv`` and, with ``-vcf 1`` in ``args``,
  ``<prefix>.vcf`` (converted by variantconvert, whose log is
  ``<prefix>.variantconvert.log``; see below).

Divergences from upstream:

- The annotations travel as one ``.tar.gz`` File, not a directory: the human
  bundle alone is a 5.3 GB download of large files. :func:`annotsv_annotsv`
  extracts it into its working directory, where AnnotSV may also write the
  overlap files it caches next to the annotations (upstream's Singularity
  runs need ``--writable-tmpfs`` for the same reason). Size the tasks' disk
  for the extracted bundle.
- ``-vcf 1`` (VCF output) doesn't work with upstream's container, so the task
  fails when it is asked for and no VCF appears, with the output listing and
  variantconvert's log; upstream silently emits no VCF. AnnotSV converts with
  the bundled variantconvert 2.0.1, which the image doesn't install and has no
  pip to install at run time, and which, run from its source, fails under the
  image's pandas 3.1.0 (``TypeError: Invalid value ... for dtype 'float64'``
  in ``vcf_from_annotsv.py``; pandas 3 turned that FutureWarning into an
  error).
- :func:`annotsv_installannotations` fails when the script installed no
  human annotations. The container's script (AnnotSV 3.5.10) only finds
  default versions from an AnnotSV Makefile it can't locate, so without
  versions it downloads nothing and still exits 0; pass both versions.
  Its Exomiser download names an unset variable, so no Exomiser data is
  installed with either form, and AnnotSV runs without Exomiser scores.

Licensing: the annotation bundle gathers third-party sources under their
own terms, among them OMIM (whose data is licensed for non-commercial use),
ClinGen, gnomAD, DECIPHER, COSMIC, GenCC and ENCODE; check them before use.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

ANNOTSV_IMAGE = "community.wave.seqera.io/library/annotsv:3.5.10--28ac84a07bcbf122"

# The bundle is downloaded, extracted and archived again: room for all three.
INSTALL_RESOURCES = flyte.Resources(cpu=2, memory="4Gi", disk="100Gi")
# Upstream labels annotation process_low; the disk holds the extracted bundle.
ANNOTATE_RESOURCES = flyte.Resources(cpu=2, memory="12Gi", disk="60Gi")

annotsv_installannotations_cmd = shell.create(
    name="annotsv_installannotations",
    image=ANNOTSV_IMAGE,
    resources=INSTALL_RESOURCES,
    inputs={"annotsv_version": str, "exomiser_version": str},
    defaults={"annotsv_version": "", "exomiser_version": ""},
    outputs={"results": Dir},
    script=r"""
        ANNOTSV_VERSION={inputs.annotsv_version}
        EXOMISER_VERSION={inputs.exomiser_version}
        W=$(mktemp -d); cd "$W"
        INSTALL_annotations.sh $ANNOTSV_VERSION $EXOMISER_VERSION || exit 1
        if [ ! -d AnnotSV_annotations/Annotations_Human ]; then
            echo "INSTALL_annotations.sh installed no human annotations" >&2
            exit 1
        fi
        tar -czf {outputs.results}/AnnotSV_annotations.tar.gz AnnotSV_annotations || exit 1
    """,
)

annotsv_annotsv_cmd = shell.create(
    name="annotsv_annotsv",
    image=ANNOTSV_IMAGE,
    resources=ANNOTATE_RESOURCES,
    inputs={
        "sv_vcf": File,
        "sv_vcf_index": File | None,
        "candidate_small_variants": File | None,
        "annotations": File,
        "candidate_genes": File | None,
        "false_positive_snv": File | None,
        "gene_transcripts": File | None,
        "prefix": str,
        "args": str,
    },
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        SV=({{inputs.sv_vcf}}); ANNOT=({{inputs.annotations}})
        shopt -s nullglob
        IDX=({{inputs.sv_vcf_index}}); SMALL=({{inputs.candidate_small_variants}}); GENES=({{inputs.candidate_genes}})
        FP=({{inputs.false_positive_snv}}); TX=({{inputs.gene_transcripts}})
        shopt -u nullglob
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        W=$(mktemp -d); cd "$W"
        for F in "${{SV[0]}}" "${{IDX[@]}}" "${{SMALL[@]}}" "${{GENES[@]}}" "${{FP[@]}}" "${{TX[@]}}"; do
            ln -s "$F" "$(basename "$F")"
        done
        # The annotation directory is the one holding Annotations_Human.
        mkdir annotations
        tar -xzf "${{ANNOT[0]}}" -C annotations || exit 1
        HUMAN=$(find annotations -maxdepth 3 -type d -name Annotations_Human | head -n 1)
        if [ -z "$HUMAN" ]; then echo "the annotations archive has no Annotations_Human directory" >&2; exit 1; fi
        ANNOTATIONS_DIR="$W/$(dirname "$HUMAN")"
        OPTIONS=()
        if [ ${{#GENES[@]}} -gt 0 ]; then OPTIONS+=(-candidateGenesFile "$(basename "${{GENES[0]}}")"); fi
        if [ ${{#SMALL[@]}} -gt 0 ]; then OPTIONS+=(-candidateSnvIndelFiles "$(basename "${{SMALL[0]}}")"); fi
        if [ ${{#FP[@]}} -gt 0 ]; then OPTIONS+=(-snvIndelFiles "$(basename "${{FP[0]}}")"); fi
        if [ ${{#TX[@]}} -gt 0 ]; then OPTIONS+=(-txFile "$(basename "${{TX[0]}}")"); fi
        AnnotSV \
            -annotationsDir "$ANNOTATIONS_DIR" \
            "${{OPTIONS[@]}}" \
            -outputFile "$P.tsv" \
            -SVinputFile "$(basename "${{SV[0]}}")" \
            $ARGS > annotsv.log 2>&1 || {{ tail -n 40 annotsv.log | tee -a /proc/1/fd/2 >&2; exit 1; }}
        cat annotsv.log
        mv ./*_AnnotSV/* .
        if [ ! -f "$P.tsv" ]; then echo "AnnotSV wrote no $P.tsv" >&2; exit 1; fi
        if [[ " $ARGS " == *" -vcf 1 "* ]] && [ ! -f "$P.vcf" ]; then
            # Also to the container log (the task's stderr is captured to a file), so the failure says why.
            {{
                echo "AnnotSV was asked for a VCF (-vcf 1) but wrote none; it wrote:"
                ls -la
                if [ -f "$P.variantconvert.log" ]; then tail -n 30 "$P.variantconvert.log"; fi
            }} 2>&1 | tee -a /proc/1/fd/2 >&2
            exit 1
        fi
        for F in "$P.tsv" "$P.unannotated.tsv" "$P.variantconvert.log" ./*.vcf; do
            if [ -f "$F" ] && [ ! -L "$F" ]; then mv "$F" {{outputs.results}}/; fi
        done
    """,
)

env = flyte.TaskEnvironment.from_task(
    "annotsv", annotsv_installannotations_cmd.as_task(), annotsv_annotsv_cmd.as_task()
)


@dataclass
class AnnotsvResult:
    tsv: File  # <prefix>.tsv
    unannotated_tsv: File | None  # <prefix>.unannotated.tsv, when AnnotSV wrote one
    vcf: File | None  # <prefix>.vcf, with -vcf 1
    results: Dir


async def annotsv_installannotations(annotsv_version: str | None = None, exomiser_version: str | None = None) -> File:
    """AnnotSV's annotations (human ``annotsv_version``, Exomiser ``exomiser_version``) as a ``.tar.gz``.

    Give both versions (for example ``"3.5"`` and ``"2406"``) or, as upstream
    allows, neither; see the module notes on why neither installs nothing
    with this container.
    """
    if bool(annotsv_version) != bool(exomiser_version):
        raise ValueError("annotsv_installannotations: give both annotsv_version and exomiser_version or neither")
    results = await annotsv_installannotations_cmd(
        annotsv_version=annotsv_version or "", exomiser_version=exomiser_version or ""
    )
    archive = await results.get_file("AnnotSV_annotations.tar.gz")
    if archive is None:
        raise FileNotFoundError("annotsv_installannotations wrote no AnnotSV_annotations.tar.gz")
    return archive


async def annotsv_annotsv(
    sv_vcf: File,
    annotations: File,
    prefix: str,
    sv_vcf_index: File | None = None,
    candidate_small_variants: File | None = None,
    candidate_genes: File | None = None,
    false_positive_snv: File | None = None,
    gene_transcripts: File | None = None,
    args: str = "",
) -> AnnotsvResult:
    """Annotate the structural variants in ``sv_vcf`` against ``annotations``.

    ``annotations`` is the archive from :func:`annotsv_installannotations` (any
    ``.tar.gz`` with an ``Annotations_Human`` directory near its top).
    """
    results = await annotsv_annotsv_cmd(
        sv_vcf=sv_vcf,
        sv_vcf_index=sv_vcf_index,
        candidate_small_variants=candidate_small_variants,
        annotations=annotations,
        candidate_genes=candidate_genes,
        false_positive_snv=false_positive_snv,
        gene_transcripts=gene_transcripts,
        prefix=prefix,
        args=args,
    )
    tsv = await results.get_file(f"{prefix}.tsv")
    if tsv is None:
        raise FileNotFoundError(f"AnnotSV wrote no {prefix}.tsv")
    return AnnotsvResult(
        tsv=tsv,
        unannotated_tsv=await results.get_file(f"{prefix}.unannotated.tsv"),
        vcf=await results.get_file(f"{prefix}.vcf"),
        results=results,
    )
