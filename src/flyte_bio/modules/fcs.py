"""fcs — NCBI Foreign Contamination Screen: adaptor and vector screening (FCS-adaptor).

Exposes :func:`fcs_fcsadaptor`, which screens an assembly (FASTA, may be
gzipped) for adaptor and vector contamination with ``av_screen_x``. ``args``
picks the screen: ``--prok`` (the default, as upstream) or ``--euk``.

As upstream, the sequence names are first cut at their first whitespace (a
work-around for https://github.com/ncbi/fcs/issues/42). Outputs keep
upstream's names:

- ``<prefix>.fcs_adaptor_report.txt``: the contamination found (a header
  line only when there is none);
- ``<prefix>.cleaned_sequences.fa.gz``: the assembly with contamination
  trimmed or masked, only when some was found;
- ``<prefix>.fcs_adaptor.log``, ``<prefix>.pipeline_args.yaml`` and
  ``<prefix>.skipped_trims.jsonl``: the run's log, settings and the trims it
  skipped.

FCS is a US government work, in the public domain.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

# NCBI's own image (upstream's); the tool reports no version.
FCS_ADAPTOR_IMAGE = "docker.io/ncbi/fcs-adaptor:0.5.0"
# Upstream's process_low.
DEFAULT_RESOURCES = flyte.Resources(cpu=2, memory="12Gi")

fcs_fcsadaptor_cmd = shell.create(
    name="fcs_fcsadaptor",
    image=FCS_ADAPTOR_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"assembly": File, "prefix": str, "args": str},
    defaults={"args": "--prok"},
    outputs={"results": Dir},
    script=r"""
        ASM=({inputs.assembly})
        P={inputs.prefix}
        ARGS={inputs.args}
        OUT={outputs.results}
        W=$(mktemp -d); cd "$W"

        # As upstream: keep only each sequence's ID (https://github.com/ncbi/fcs/issues/42).
        gzip -cdf "${ASM[0]}" | sed -E 's/^(>[^[:space:]]+).*/\1/' > "$P.id.only.fasta" || exit 1

        av_screen_x -o output/ $ARGS "$P.id.only.fasta" || exit 1

        if [ "$(wc -l < output/fcs_adaptor_report.txt)" -gt 1 ]; then
            find output/cleaned_sequences/ -type f ! -name "*.gz" -exec gzip {} \;
            cp output/cleaned_sequences/* "$OUT/$P.cleaned_sequences.fa.gz"
        fi
        cp output/fcs_adaptor_report.txt "$OUT/$P.fcs_adaptor_report.txt"
        cp output/fcs_adaptor.log "$OUT/$P.fcs_adaptor.log"
        cp output/pipeline_args.yaml "$OUT/$P.pipeline_args.yaml"
        cp output/skipped_trims.jsonl "$OUT/$P.skipped_trims.jsonl"
    """,
)

env = flyte.TaskEnvironment.from_task("fcs", fcs_fcsadaptor_cmd.as_task())


@dataclass
class FcsAdaptorResult:
    adaptor_report: File  # <prefix>.fcs_adaptor_report.txt
    cleaned_assembly: File | None  # <prefix>.cleaned_sequences.fa.gz, when contamination was found
    log: File  # <prefix>.fcs_adaptor.log
    pipeline_args: File  # <prefix>.pipeline_args.yaml
    skipped_trims: File  # <prefix>.skipped_trims.jsonl


async def fcs_fcsadaptor(assembly: File, prefix: str = "sample", args: str = "--prok") -> FcsAdaptorResult:
    """Screen ``assembly`` for adaptor and vector contamination (``--prok`` or ``--euk`` in ``args``)."""
    results = await fcs_fcsadaptor_cmd(assembly=assembly, prefix=prefix, args=args)

    async def pick(name: str) -> File:
        found = await results.get_file(name)
        if found is None:
            raise FileNotFoundError(f"FCS-adaptor wrote no {name}")
        return found

    return FcsAdaptorResult(
        adaptor_report=await pick(f"{prefix}.fcs_adaptor_report.txt"),
        cleaned_assembly=await results.get_file(f"{prefix}.cleaned_sequences.fa.gz"),
        log=await pick(f"{prefix}.fcs_adaptor.log"),
        pipeline_args=await pick(f"{prefix}.pipeline_args.yaml"),
        skipped_trims=await pick(f"{prefix}.skipped_trims.jsonl"),
    )
