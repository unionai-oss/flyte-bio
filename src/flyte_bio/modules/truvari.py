"""truvari — Truvari, structural variant benchmarking.

Exposes :func:`truvari_bench`, which compares a VCF's SVs (``--comp``) with a
truth VCF's (``--base``) against a reference, optionally within a BED's
regions (``--includebed``). As upstream, Truvari's outputs are renamed
``<prefix>.<name>``:

- ``<prefix>.tp-base.vcf.gz`` / ``<prefix>.tp-comp.vcf.gz``: the truth's and
  the calls' matched SVs;
- ``<prefix>.fn.vcf.gz`` / ``<prefix>.fp.vcf.gz``: the truth's SVs the calls
  miss, and the calls not in the truth;
- each with its ``.tbi``; ``<prefix>.summary.json`` (precision, recall, F1 …)
  and ``<prefix>.log.txt``.

The VCFs, their indexes and the reference with its ``.fai`` are linked side by
side. Matching options (``--refdist``, ``--pctseq`` …) go in ``args``.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

TRUVARI_IMAGE = "community.wave.seqera.io/library/truvari:5.4.0--a783dff809870f2f"
# Upstream's process_single.
DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory="6Gi")
OUTPUTS = ["fn.vcf.gz", "fp.vcf.gz", "tp-base.vcf.gz", "tp-comp.vcf.gz"]

truvari_bench_cmd = shell.create(
    name="truvari_bench",
    image=TRUVARI_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "vcf": File,
        "tbi": File,
        "truth_vcf": File,
        "truth_tbi": File,
        "bed": File | None,
        "fasta": File,
        "fai": File,
        "prefix": str,
        "args": str,
    },
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        VCF=({inputs.vcf}); TBI=({inputs.tbi}); TRUTH=({inputs.truth_vcf}); TRUTH_TBI=({inputs.truth_tbi})
        FA=({inputs.fasta}); FAI=({inputs.fai})
        shopt -s nullglob; BED=({inputs.bed}); shopt -u nullglob
        P={inputs.prefix}
        ARGS={inputs.args}
        OUT={outputs.results}
        W=$(mktemp -d); cd "$W"
        for f in "${VCF[0]}" "${TBI[0]}" "${TRUTH[0]}" "${TRUTH_TBI[0]}" "${FA[0]}" "${FAI[0]}" "${BED[@]}"; do
            ln -s "$f" "$(basename "$f")"
        done
        REGIONS=()
        if [ ${#BED[@]} -gt 0 ]; then REGIONS=(--includebed "$(basename "${BED[0]}")"); fi
        truvari bench \
            --base "$(basename "${TRUTH[0]}")" \
            --comp "$(basename "${VCF[0]}")" \
            --reference "$(basename "${FA[0]}")" \
            --output "$P" \
            "${REGIONS[@]}" \
            $ARGS || exit 1
        for f in fn.vcf.gz fp.vcf.gz tp-base.vcf.gz tp-comp.vcf.gz; do
            mv "$P/$f" "$OUT/$P.$f"
            mv "$P/$f.tbi" "$OUT/$P.$f.tbi"
        done
        mv "$P/summary.json" "$OUT/$P.summary.json"
        mv "$P/log.txt" "$OUT/$P.log.txt"
    """,
)

env = flyte.TaskEnvironment.from_task("truvari", truvari_bench_cmd.as_task())


@dataclass
class TruvariBenchResult:
    fn_vcf: File  # <prefix>.fn.vcf.gz (+ fn_tbi)
    fn_tbi: File
    fp_vcf: File  # <prefix>.fp.vcf.gz (+ fp_tbi)
    fp_tbi: File
    tp_base_vcf: File  # <prefix>.tp-base.vcf.gz (+ tp_base_tbi)
    tp_base_tbi: File
    tp_comp_vcf: File  # <prefix>.tp-comp.vcf.gz (+ tp_comp_tbi)
    tp_comp_tbi: File
    summary: File  # <prefix>.summary.json
    log: File  # <prefix>.log.txt


async def truvari_bench(
    vcf: File,
    tbi: File,
    truth_vcf: File,
    truth_tbi: File,
    fasta: File,
    fai: File,
    prefix: str,
    bed: File | None = None,
    args: str = "",
) -> TruvariBenchResult:
    """Compare ``vcf``'s SVs with ``truth_vcf``'s against ``fasta``, within ``bed`` when given."""
    results = await truvari_bench_cmd(
        vcf=vcf,
        tbi=tbi,
        truth_vcf=truth_vcf,
        truth_tbi=truth_tbi,
        bed=bed,
        fasta=fasta,
        fai=fai,
        prefix=prefix,
        args=args,
    )

    async def pick(name: str) -> File:
        found = await results.get_file(f"{prefix}.{name}")
        if found is None:
            raise FileNotFoundError(f"truvari bench wrote no {prefix}.{name}")
        return found

    files = {}
    for name in OUTPUTS:
        key = name.removesuffix(".vcf.gz").replace("-", "_")
        files[f"{key}_vcf"] = await pick(name)
        files[f"{key}_tbi"] = await pick(name + ".tbi")
    return TruvariBenchResult(**files, summary=await pick("summary.json"), log=await pick("log.txt"))
