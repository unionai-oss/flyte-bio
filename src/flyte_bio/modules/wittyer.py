"""wittyer — witty.er, Illumina's structural variant and CNV benchmarking.

Exposes :func:`wittyer`, which compares a (plain) VCF's SVs and CNVs with a
truth VCF's, by type and size bin, optionally within a BED's regions (``-b``)
or with a JSON config of per-type matching settings (``-c``; not both),
writing as upstream:

- ``<prefix>.json``: the statistics (``Wittyer.Stats.json``);
- ``<prefix>.vcf.gz`` and its ``.tbi``: the truth's and the calls' variants,
  annotated with their match status.

The VCFs, BED and config are linked into the working directory and named
relative to it, as upstream. witty.er is a .NET tool, run from the image as
upstream (``dotnet /opt/Wittyer/Wittyer.dll``), with ``HOME`` in the working
directory.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

WITTYER_IMAGE = "quay.io/nf-core/wittyer:0.5.2"
# Upstream's process_single.
DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory="6Gi")

wittyer_cmd = shell.create(
    name="wittyer",
    image=WITTYER_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "query_vcf": File,
        "truth_vcf": File,
        "bed": File | None,
        "config": File | None,
        "prefix": str,
        "args": str,
    },
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        QUERY=({inputs.query_vcf}); TRUTH=({inputs.truth_vcf})
        shopt -s nullglob; BED=({inputs.bed}); CONFIG=({inputs.config}); shopt -u nullglob
        P={inputs.prefix}
        ARGS={inputs.args}
        OUT={outputs.results}
        W=$(mktemp -d); cd "$W"
        export HOME="$W"
        for f in "${QUERY[0]}" "${TRUTH[0]}" "${BED[@]}" "${CONFIG[@]}"; do ln -s "$f" "$(basename "$f")"; done
        OPTIONAL=()
        if [ ${#BED[@]} -gt 0 ]; then OPTIONAL+=(-b "$(basename "${BED[0]}")"); fi
        if [ ${#CONFIG[@]} -gt 0 ]; then OPTIONAL+=(-c "$(basename "${CONFIG[0]}")"); fi
        mkdir bench
        dotnet /opt/Wittyer/Wittyer.dll \
            -t "$(basename "${TRUTH[0]}")" \
            -i "$(basename "${QUERY[0]}")" \
            -o bench \
            "${OPTIONAL[@]}" \
            $ARGS || exit 1
        mv bench/Wittyer.Stats.json "$OUT/$P.json"
        mv bench/*.vcf.gz "$OUT/$P.vcf.gz"
        mv bench/*.vcf.gz.tbi "$OUT/$P.vcf.gz.tbi"
    """,
)

env = flyte.TaskEnvironment.from_task("wittyer", wittyer_cmd.as_task())


@dataclass
class WittyerResult:
    report: File  # <prefix>.json
    bench_vcf: File  # <prefix>.vcf.gz
    bench_vcf_tbi: File  # <prefix>.vcf.gz.tbi


async def wittyer(
    query_vcf: File,
    truth_vcf: File,
    prefix: str,
    bed: File | None = None,
    config: File | None = None,
    args: str = "",
) -> WittyerResult:
    """Compare ``query_vcf``'s SVs with ``truth_vcf``'s, within ``bed`` or with ``config`` when given."""
    if bed is not None and config is not None:
        raise ValueError("witty.er takes a BED or a config, not both")
    if f"{prefix}.vcf" in (query_vcf.name, truth_vcf.name):
        raise ValueError("input and output names are the same; choose another prefix")
    results = await wittyer_cmd(
        query_vcf=query_vcf, truth_vcf=truth_vcf, bed=bed, config=config, prefix=prefix, args=args
    )

    async def pick(name: str) -> File:
        found = await results.get_file(name)
        if found is None:
            raise FileNotFoundError(f"witty.er wrote no {name}")
        return found

    return WittyerResult(
        report=await pick(f"{prefix}.json"),
        bench_vcf=await pick(f"{prefix}.vcf.gz"),
        bench_vcf_tbi=await pick(f"{prefix}.vcf.gz.tbi"),
    )
