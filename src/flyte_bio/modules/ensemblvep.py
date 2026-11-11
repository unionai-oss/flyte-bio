"""ensemblvep — variant annotation with Ensembl's Variant Effect Predictor, with its cache download.

Exposes:

- :func:`ensemblvep_download`, which downloads the VEP cache of one species,
  assembly and cache version (``vep_install``) into a Dir laid out as VEP's
  ``--dir_cache`` expects: ``<species>/<cache_version>_<assembly>/``. It's a
  task of its own, cached on its inputs, so each cache downloads once. With
  ``preflight_check`` it first checks that Ensembl's ``CHECKSUMS`` list the
  cache's tarball, so a misspelt species fails fast.
- :func:`ensemblvep_vep`, which annotates one VCF against that cache (or a
  sorted, bgzipped GTF with its ``.tbi``, instead of a cache) and writes:

  - ``<prefix>.vcf.gz`` with its tabix index ``<prefix>.vcf.gz.tbi``
    (``--vcf``, the default), ``<prefix>.tab.gz`` (``--tab``) or
    ``<prefix>.json.gz`` (``--json``), compressed with bgzip unless ``args``
    sets ``--compress_output`` itself;
  - VEP's summary report, ``<prefix>.<format>.gz_summary.html``.

  As upstream, the output format is picked by which of ``--vcf``, ``--json``
  and ``--tab`` ``args`` contains (in that order), and ``args2`` goes to tabix.
  ``fasta`` is the reference (``--fasta``). ``custom_extra_files`` (files for
  ``--custom``) and ``extra_files`` (files for plugins) are linked into the
  working directory under their own names, so ``args`` refers to them by
  name. Without a cache or GTF, VEP looks in ``/.vep``. ``--fork`` is the
  task's CPUs.

The download's default ``args`` are upstream's test arguments, which install
only the cache, without prompting. Upstream's module has no default; without
``--AUTO`` ``vep_install`` asks questions on stdin.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

VEP_IMAGE = "community.wave.seqera.io/library/ensembl-vep_perl-math-cdf_htslib_gzip_tar:43f77a48cc79bf38"
DOWNLOAD_ARGS = "--AUTO c --NO_BIOPERL --NO_HTSLIB --NO_TEST --NO_UPDATE"
# Small caches (the C. elegans test one is about 50 MB) fit easily; the human
# cache is a ~25 GB tarball that unpacks to more, so raise disk with
# ensemblvep_download_cmd.override(resources=...).
DOWNLOAD_RESOURCES = flyte.Resources(cpu=1, memory="4Gi", disk="20Gi")
# Upstream's label asks for 6 CPUs and 36 GB. Each VEP fork holds its share of
# the cache in memory, so a human cache with many forks needs more memory.
DEFAULT_RESOURCES = flyte.Resources(cpu=4, memory="8Gi", disk="20Gi")

ensemblvep_download_cmd = shell.create(
    name="ensemblvep_download",
    image=VEP_IMAGE,
    resources=DOWNLOAD_RESOURCES,
    inputs={"assembly": str, "species": str, "cache_version": str, "preflight_check": bool, "args": str},
    defaults={"preflight_check": True, "args": DOWNLOAD_ARGS},
    outputs={"cache": Dir},
    script=r"""
        ASSEMBLY={inputs.assembly}
        SPECIES={inputs.species}
        VERSION={inputs.cache_version}
        ARGS={inputs.args}
        OUT={outputs.cache}
        W=$(mktemp -d); cd "$W"

        FILENAME="${SPECIES}_vep_${VERSION}_${ASSEMBLY}.tar.gz"
        CHECKSUMS_URL="https://ftp.ensembl.org/pub/release-${VERSION}/variation/indexed_vep_cache/CHECKSUMS"
        # Both steps are perl, whose die exits 255, which the cluster takes for an
        # interruption and retries without end; exit 1 instead, a plain failure.
        if [ {inputs.preflight_check} = true ]; then
            FILENAME="$FILENAME" CHECKSUMS_URL="$CHECKSUMS_URL" perl -MHTTP::Tiny -e '
                my $r = HTTP::Tiny->new(timeout => 30)->get($ENV{CHECKSUMS_URL});
                $r->{success} or die "Failed to fetch CHECKSUMS (HTTP $r->{status})\n";
                $r->{content} =~ /\Q$ENV{FILENAME}\E/ or die "$ENV{FILENAME} not found in CHECKSUMS\n";
                print "Pre-flight OK: $ENV{FILENAME} found in CHECKSUMS\n";
            ' || exit 1
        fi

        vep_install \
            --CACHEDIR "$OUT" \
            --SPECIES "$SPECIES" \
            --ASSEMBLY "$ASSEMBLY" \
            --CACHE_VERSION "$VERSION" \
            $ARGS \
            || exit 1
    """,
)

ensemblvep_vep_cmd = shell.create(
    name="ensemblvep_vep",
    image=VEP_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "vcf": File,
        "custom_extra_files": list[File],
        "genome": str,
        "species": str,
        "cache_version": str,
        "cache": Dir | None,
        "fasta": File | None,
        "extra_files": list[File],
        "gtf": File | None,
        "gtf_tbi": File | None,
        "prefix": str,
        "args": str,
        "args2": str,
    },
    defaults={"args": "", "args2": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        VCF=({{inputs.vcf}})
        shopt -s nullglob
        CUSTOM=({{inputs.custom_extra_files}}); EXTRA_FILES=({{inputs.extra_files}})
        FA=({{inputs.fasta}}); GTF=({{inputs.gtf}}); GTF_TBI=({{inputs.gtf_tbi}})
        shopt -u nullglob
        CACHE={{inputs.cache}}
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        ARGS2={{inputs.args2}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"

        VCF_NAME=$(basename "${{VCF[0]}}"); ln -s "${{VCF[0]}}" "$VCF_NAME"
        for f in "${{CUSTOM[@]}}" "${{EXTRA_FILES[@]}}"; do ln -sf "$f" "$(basename "$f")"; done

        EXT=vcf
        if [[ "$ARGS" == *--vcf* ]]; then EXT=vcf
        elif [[ "$ARGS" == *--json* ]]; then EXT=json
        elif [[ "$ARGS" == *--tab* ]]; then EXT=tab
        fi
        EXTRA=()
        if [[ "$ARGS" != *--compress_output* ]]; then EXTRA+=(--compress_output bgzip); fi
        if [ ${{#FA[@]}} -gt 0 ]; then
            FA_NAME=$(basename "${{FA[0]}}"); ln -s "${{FA[0]}}" "$FA_NAME"; EXTRA+=(--fasta "$FA_NAME")
        fi
        if [ ${{#GTF[@]}} -gt 0 ]; then
            GTF_NAME=$(basename "${{GTF[0]}}"); ln -s "${{GTF[0]}}" "$GTF_NAME"
            if [ ${{#GTF_TBI[@]}} -gt 0 ]; then ln -s "${{GTF_TBI[0]}}" "$GTF_NAME.tbi"; fi
            EXTRA+=(--gtf "$GTF_NAME")
        else
            DIR_CACHE=/.vep
            if [ -d "$CACHE" ]; then DIR_CACHE="$CACHE"; fi
            EXTRA+=(--cache --cache_version {{inputs.cache_version}} --dir_cache "$DIR_CACHE")
        fi

        vep \
            -i "$VCF_NAME" \
            -o "$P.$EXT.gz" \
            $ARGS \
            "${{EXTRA[@]}}" \
            --assembly {{inputs.genome}} \
            --species {{inputs.species}} \
            --fork $TASK_CPUS

        if [ $EXT = vcf ]; then tabix $ARGS2 "$P.$EXT.gz"; fi
        shopt -s nullglob
        for f in "$P.$EXT.gz" ./*.tbi ./*.html ./*_warnings.txt; do mv "$f" "$OUT"/; done
    """,
)

env = flyte.TaskEnvironment.from_task("ensemblvep", ensemblvep_download_cmd.as_task(), ensemblvep_vep_cmd.as_task())


async def ensemblvep_download(
    assembly: str,
    species: str,
    cache_version: str,
    preflight_check: bool = True,
    args: str = DOWNLOAD_ARGS,
) -> Dir:
    """The VEP cache of ``species`` / ``assembly`` / ``cache_version``: a Dir to pass as ``cache``."""
    return await ensemblvep_download_cmd(
        assembly=assembly, species=species, cache_version=cache_version, preflight_check=preflight_check, args=args
    )


@dataclass
class EnsemblVepResult:
    vcf: File | None  # with --vcf (the default)
    tbi: File | None  # the tabix index of vcf
    tab: File | None  # with --tab
    json: File | None  # with --json
    report: File | None  # VEP's summary HTML; none with --no_stats
    results: Dir  # everything VEP wrote


async def ensemblvep_vep(
    vcf: File,
    genome: str,
    species: str,
    cache_version: str,
    cache: Dir | None = None,
    fasta: File | None = None,
    custom_extra_files: list[File] | None = None,
    extra_files: list[File] | None = None,
    gtf: File | None = None,
    gtf_tbi: File | None = None,
    prefix: str = "sample",
    args: str = "",
    args2: str = "",
) -> EnsemblVepResult:
    """Annotate ``vcf`` against ``cache`` (or ``gtf``), in the format ``args`` picks (``--vcf`` by default)."""
    results = await ensemblvep_vep_cmd(
        vcf=vcf,
        custom_extra_files=custom_extra_files or [],
        genome=genome,
        species=species,
        cache_version=cache_version,
        cache=cache,
        fasta=fasta,
        extra_files=extra_files or [],
        gtf=gtf,
        gtf_tbi=gtf_tbi,
        prefix=prefix,
        args=args,
        args2=args2,
    )
    if "--vcf" in args:
        ext = "vcf"
    elif "--json" in args:
        ext = "json"
    elif "--tab" in args:
        ext = "tab"
    else:
        ext = "vcf"
    output = await results.get_file(f"{prefix}.{ext}.gz")
    if output is None:
        raise FileNotFoundError(f"VEP wrote no {prefix}.{ext}.gz")
    html = [f for f in await results.list_files() if f.name.endswith(".html")]
    return EnsemblVepResult(
        vcf=output if ext == "vcf" else None,
        tbi=await results.get_file(f"{prefix}.vcf.gz.tbi") if ext == "vcf" else None,
        tab=output if ext == "tab" else None,
        json=output if ext == "json" else None,
        report=html[0] if html else None,
        results=results,
    )
