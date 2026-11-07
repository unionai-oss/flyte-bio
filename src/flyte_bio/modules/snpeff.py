"""snpeff — variant annotation with snpEff, with its database download.

Exposes:

- :func:`snpeff_download`, which downloads one snpEff database (``snpEff
  download <db>``, e.g. ``WBcel235.99`` or ``GRCh38.105``) into a Dir laid
  out as snpEff's ``-dataDir`` expects: ``<db>/snpEffectPredictor.bin``
  (upstream's ``snpeff_cache``). It's a task of its own, cached on its
  inputs, so each database downloads once.
- :func:`snpeff_snpeff`, which annotates one VCF with database ``db`` from
  that cache (``-dataDir``) and writes, under upstream's names:

  - ``<prefix>.ann.vcf``: the annotated VCF (snpEff's stdout);
  - ``<prefix>.csv``: the statistics (``-csvStats``);
  - ``snpEff_summary.html``: the HTML summary;
  - ``<prefix>.genes.txt``: the per-gene counts.

  Without a cache, snpEff uses its configured data directory, and, as
  usual, downloads ``db`` there itself when it's missing.

As upstream, snpEff's Java heap is 80% of the task's memory.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import JAVA_HEAP_MB, TASK_LIMITS

SNPEFF_IMAGE = "community.wave.seqera.io/library/snpeff:5.4.0c--e08ddc54579e82bd"
# Upstream labels both process_medium (36 GB). Small databases (the C. elegans
# test one) fit in a few GB; GRCh38's needs about 16 GB of heap to annotate.
DOWNLOAD_RESOURCES = flyte.Resources(cpu=1, memory="4Gi", disk="20Gi")
DEFAULT_RESOURCES = flyte.Resources(cpu=2, memory="8Gi")

snpeff_download_cmd = shell.create(
    name="snpeff_download",
    image=SNPEFF_IMAGE,
    resources=DOWNLOAD_RESOURCES,
    inputs={"db": str, "args": str},
    defaults={"args": ""},
    outputs={"cache": Dir},
    script=rf"""
        {TASK_LIMITS}
        DB={{inputs.db}}
        ARGS={{inputs.args}}
        OUT={{outputs.cache}}
        W=$(mktemp -d); cd "$W"
        snpEff \
            -Xmx{JAVA_HEAP_MB}M \
            download "$DB" \
            -dataDir "$OUT" \
            $ARGS
        if [ ! -d "$OUT/$DB" ]; then
            echo "snpEff downloaded no $DB" >&2
            exit 1
        fi
    """,
)

snpeff_snpeff_cmd = shell.create(
    name="snpeff_snpeff",
    image=SNPEFF_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"vcf": File, "db": str, "cache": Dir | None, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        VCF=({{inputs.vcf}})
        DB={{inputs.db}}
        CACHE={{inputs.cache}}
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"

        VCF_NAME=$(basename "${{VCF[0]}}"); ln -s "${{VCF[0]}}" "$VCF_NAME"
        CACHE_OPTION=()
        if [ -d "$CACHE" ]; then CACHE_OPTION=(-dataDir "$CACHE"); fi

        snpEff \
            -Xmx{JAVA_HEAP_MB}M \
            -XX:-UsePerfData \
            "$DB" \
            $ARGS \
            -csvStats "$P.csv" \
            "${{CACHE_OPTION[@]}}" \
            "$VCF_NAME" \
            > "$P.ann.vcf"

        shopt -s nullglob
        for f in ./*.ann.vcf ./*.csv ./*.html ./*.genes.txt; do mv "$f" "$OUT"/; done
    """,
)

env = flyte.TaskEnvironment.from_task("snpeff", snpeff_download_cmd.as_task(), snpeff_snpeff_cmd.as_task())


async def snpeff_download(db: str, args: str = "") -> Dir:
    """The snpEff database ``db`` (e.g. ``WBcel235.99``): a Dir to pass as ``cache``."""
    return await snpeff_download_cmd(db=db, args=args)


@dataclass
class SnpEffResult:
    vcf: File  # <prefix>.ann.vcf
    report: File | None  # <prefix>.csv, the statistics
    summary_html: File | None  # snpEff_summary.html; none with -noStats
    genes_txt: File | None  # <prefix>.genes.txt
    results: Dir  # everything snpEff wrote


async def snpeff_snpeff(
    vcf: File,
    db: str,
    cache: Dir | None = None,
    prefix: str = "sample",
    args: str = "",
) -> SnpEffResult:
    """Annotate ``vcf`` with snpEff database ``db``, read from ``cache`` (a :func:`snpeff_download` result)."""
    results = await snpeff_snpeff_cmd(vcf=vcf, db=db, cache=cache, prefix=prefix, args=args)
    files = await results.list_files()

    def ending(suffix: str) -> File | None:
        matches = [f for f in files if f.name.endswith(suffix)]
        return matches[0] if matches else None

    annotated = await results.get_file(f"{prefix}.ann.vcf")
    if annotated is None:
        raise FileNotFoundError(f"snpEff wrote no {prefix}.ann.vcf")
    return SnpEffResult(
        vcf=annotated,
        report=await results.get_file(f"{prefix}.csv"),
        summary_html=ending(".html"),
        genes_txt=ending(".genes.txt"),
        results=results,
    )
