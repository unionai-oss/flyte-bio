"""sniffles — Sniffles2 structural-variant calling and genotyping from long reads.

Exposes :func:`sniffles`, which calls structural variants from one sample's
sorted, indexed long-read BAM or CRAM, optionally with its reference
(``--reference``, for insertion sequences and CRAM) and a tandem-repeat BED
(``--tandem-repeats``). It writes the calls (``vcf_output``:
``<prefix>.vcf.gz`` and its ``.tbi``), the SV candidates for later
multi-sample calling (``snf_output``: ``<prefix>.snf``), or both. With
``genotype_vcf`` it instead genotypes the sites of a known SV VCF in the
sample (``--genotype-vcf``).

As upstream:

- Inputs keep their names in the working directory; Sniffles records its
  command line and the date in the VCF header.
- A gzipped ``genotype_vcf`` is decompressed to ``<prefix>.sites.vcf`` first,
  as Sniffles reads it as plain text.
- ``-t`` follows the task's CPUs; options such as ``--minsupport`` go in
  ``args``.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

SNIFFLES_IMAGE = "community.wave.seqera.io/library/htslib_sniffles:14a95caca8d8452b"

# Upstream labels it process_high; these fit test-sized and targeted runs.
DEFAULT_RESOURCES = flyte.Resources(cpu=4, memory="16Gi")

sniffles_cmd = shell.create(
    name="sniffles",
    image=SNIFFLES_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "alignment": File,
        "index": File | None,
        "genotype_vcf": File | None,
        "fasta": File | None,
        "tandem_file": File | None,
        "vcf_output": bool,
        "snf_output": bool,
        "prefix": str,
        "args": str,
    },
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        ALN=({{inputs.alignment}})
        shopt -s nullglob
        IDX=({{inputs.index}}); SITES=({{inputs.genotype_vcf}}); FA=({{inputs.fasta}}); TR=({{inputs.tandem_file}})
        shopt -u nullglob
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        W=$(mktemp -d); cd "$W"
        for F in "${{ALN[0]}}" "${{IDX[@]}}" "${{SITES[@]}}" "${{FA[@]}}" "${{TR[@]}}"; do
            ln -s "$F" "$(basename "$F")"
        done
        REF=(); REPEATS=(); GENOTYPE=()
        if [ ${{#FA[@]}} -gt 0 ]; then REF=(--reference "$(basename "${{FA[0]}}")"); fi
        if [ ${{#TR[@]}} -gt 0 ]; then REPEATS=(--tandem-repeats "$(basename "${{TR[0]}}")"); fi
        if [ ${{#SITES[@]}} -gt 0 ]; then
            SITES_NAME=$(basename "${{SITES[0]}}")
            if [[ "$SITES_NAME" == *.gz ]]; then
                # Sniffles reads the sites file as plain text.
                bgzip -dc "$SITES_NAME" > "$P.sites.vcf" || exit 1
                SITES_NAME="$P.sites.vcf"
            fi
            GENOTYPE=(--genotype-vcf "$SITES_NAME")
        fi
        OUTS=()
        if [ {{inputs.vcf_output}} = true ]; then OUTS+=(--vcf "$P.vcf.gz"); fi
        if [ {{inputs.snf_output}} = true ]; then OUTS+=(--snf "$P.snf"); fi
        sniffles \
            --input "$(basename "${{ALN[0]}}")" \
            "${{REF[@]}}" \
            -t $TASK_CPUS \
            "${{REPEATS[@]}}" \
            "${{GENOTYPE[@]}}" \
            "${{OUTS[@]}}" \
            $ARGS || exit 1
        for F in "$P.vcf.gz" "$P.vcf.gz.tbi" "$P.snf"; do
            if [ -f "$F" ] && [ ! -L "$F" ]; then mv "$F" {{outputs.results}}/; fi
        done
    """,
)

env = flyte.TaskEnvironment.from_task("sniffles", sniffles_cmd.as_task())


@dataclass
class SnifflesResult:
    vcf: File | None  # <prefix>.vcf.gz, with vcf_output
    tbi: File | None  # <prefix>.vcf.gz.tbi, with vcf_output
    snf: File | None  # <prefix>.snf, with snf_output
    results: Dir


async def sniffles(
    alignment: File,
    prefix: str,
    index: File | None = None,
    genotype_vcf: File | None = None,
    fasta: File | None = None,
    tandem_file: File | None = None,
    vcf_output: bool = True,
    snf_output: bool = False,
    args: str = "",
) -> SnifflesResult:
    """Call (or, with ``genotype_vcf``, genotype) structural variants in ``alignment``."""
    if not (vcf_output or snf_output):
        raise ValueError("sniffles: ask for vcf_output, snf_output or both")
    results = await sniffles_cmd(
        alignment=alignment,
        index=index,
        genotype_vcf=genotype_vcf,
        fasta=fasta,
        tandem_file=tandem_file,
        vcf_output=vcf_output,
        snf_output=snf_output,
        prefix=prefix,
        args=args,
    )
    found = {}
    for key, name, wanted in (
        ("vcf", f"{prefix}.vcf.gz", vcf_output),
        ("tbi", f"{prefix}.vcf.gz.tbi", vcf_output),
        ("snf", f"{prefix}.snf", snf_output),
    ):
        found[key] = await results.get_file(name) if wanted else None
        if wanted and found[key] is None:
            raise FileNotFoundError(f"sniffles wrote no {name}")
    return SnifflesResult(**found, results=results)
