"""freyja — lineage abundances from mixed viral samples (e.g. wastewater).

Exposes:

- :func:`freyja_update` — downloads the current lineage barcodes and
  lineage metadata into a database directory (``--pathogen`` in ``args``
  picks a pathogen other than SARS-CoV-2). The data is live, so its
  contents change over time. Returns the directory and picks out the
  barcodes (``*barcodes.*``), the lineage topology (``*lineages.yml``), the
  pathogen config (``*pathogen_config.yml``) and, when written, the curated
  lineage metadata (``*curated_lineages.json``).
- :func:`freyja_variants` — variant calls and per-site depths of a sorted
  BAM against the reference: ``<prefix>.variants.tsv`` and
  ``<prefix>.depth.tsv``.
- :func:`freyja_demix` — the lineage abundances of one sample, from its
  variants and depths and the barcodes: ``<prefix>.tsv``.
- :func:`freyja_boot` — bootstrapped confidence intervals for those
  abundances, ``repeats`` times: ``<prefix>_lineages.csv`` and
  ``<prefix>_summarized.csv``. Bootstrap threads follow the task's CPUs.

demix and boot take the lineage metadata and topology optionally, as
upstream.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

FREYJA_IMAGE = "community.wave.seqera.io/library/freyja:2.0.6--be3cbd6f14f9d542"
# Upstream: update process_single, variants process_medium, demix process_low, boot process_high.
UPDATE_RESOURCES = flyte.Resources(cpu=1, memory="8Gi")
VARIANTS_RESOURCES = flyte.Resources(cpu=1, memory="4Gi")
DEMIX_RESOURCES = flyte.Resources(cpu=1, memory="8Gi")
BOOT_RESOURCES = flyte.Resources(cpu=4, memory="16Gi")

freyja_update_cmd = shell.create(
    name="freyja_update",
    image=FREYJA_IMAGE,
    resources=UPDATE_RESOURCES,
    inputs={"args": str},
    defaults={"args": ""},
    outputs={"db": Dir},
    script=r"""
        ARGS={inputs.args}
        freyja \
            update \
            $ARGS \
            --outdir {outputs.db} || exit 1
    """,
)

freyja_variants_cmd = shell.create(
    name="freyja_variants",
    image=FREYJA_IMAGE,
    resources=VARIANTS_RESOURCES,
    inputs={"bam": File, "fasta": File, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        BAM=({inputs.bam}); FA=({inputs.fasta})
        P={inputs.prefix}
        ARGS={inputs.args}
        OUT={outputs.results}
        W=$(mktemp -d); cd "$W"
        ln -s "${BAM[0]}" "$(basename "${BAM[0]}")"
        ln -s "${FA[0]}" "$(basename "${FA[0]}")"
        freyja \
            variants \
            $ARGS \
            --ref "$(basename "${FA[0]}")" \
            --variants "$P.variants.tsv" \
            --depths "$P.depth.tsv" \
            "$(basename "${BAM[0]}")" || exit 1
        mv "$P.variants.tsv" "$P.depth.tsv" "$OUT"/
    """,
)

# The optional lineage metadata and topology, as freyja's options.
LINEAGE_OPTIONS = r"""
        shopt -s nullglob
        LINEAGES_META=({inputs.lineages_meta}); TOPOLOGY=({inputs.lineages_topology})
        shopt -u nullglob
        LINEAGE_OPTIONS=()
        if [ ${#LINEAGES_META[@]} -gt 0 ]; then LINEAGE_OPTIONS+=(--meta "${LINEAGES_META[0]}"); fi
        if [ ${#TOPOLOGY[@]} -gt 0 ]; then LINEAGE_OPTIONS+=(--lineageyml "${TOPOLOGY[0]}"); fi
"""

freyja_demix_cmd = shell.create(
    name="freyja_demix",
    image=FREYJA_IMAGE,
    resources=DEMIX_RESOURCES,
    inputs={
        "variants": File,
        "depths": File,
        "barcodes": File,
        "lineages_meta": File | None,
        "lineages_topology": File | None,
        "prefix": str,
        "args": str,
    },
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        VARIANTS=({inputs.variants}); DEPTHS=({inputs.depths}); BARCODES=({inputs.barcodes})
        """
    + LINEAGE_OPTIONS
    + r"""
        P={inputs.prefix}
        ARGS={inputs.args}
        OUT={outputs.results}
        W=$(mktemp -d); cd "$W"
        freyja \
            demix \
            $ARGS \
            --output "$P.tsv" \
            --barcodes "${BARCODES[0]}" \
            "${LINEAGE_OPTIONS[@]}" \
            "${VARIANTS[0]}" \
            "${DEPTHS[0]}" || exit 1
        mv "$P.tsv" "$OUT"/
    """,
)

freyja_boot_cmd = shell.create(
    name="freyja_boot",
    image=FREYJA_IMAGE,
    resources=BOOT_RESOURCES,
    inputs={
        "variants": File,
        "depths": File,
        "repeats": int,
        "barcodes": File,
        "lineages_meta": File | None,
        "lineages_topology": File | None,
        "prefix": str,
        "args": str,
    },
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        VARIANTS=({{inputs.variants}}); DEPTHS=({{inputs.depths}}); BARCODES=({{inputs.barcodes}})
        """
    + LINEAGE_OPTIONS
    + r"""
        P={inputs.prefix}
        ARGS={inputs.args}
        OUT={outputs.results}
        W=$(mktemp -d); cd "$W"
        freyja \
            boot \
            $ARGS \
            --nt $TASK_CPUS \
            --nb {inputs.repeats} \
            --output_base "$P" \
            --barcodes "${BARCODES[0]}" \
            "${LINEAGE_OPTIONS[@]}" \
            "${VARIANTS[0]}" \
            "${DEPTHS[0]}" || exit 1
        mv "${P}_lineages.csv" "${P}_summarized.csv" "$OUT"/
    """,
)

env = flyte.TaskEnvironment.from_task(
    "freyja",
    freyja_update_cmd.as_task(),
    freyja_variants_cmd.as_task(),
    freyja_demix_cmd.as_task(),
    freyja_boot_cmd.as_task(),
)


async def pick_file(results: Dir, name: str, tool: str) -> File:
    found = await results.get_file(name)
    if found is None:
        raise FileNotFoundError(f"{tool} wrote no {name}")
    return found


@dataclass
class FreyjaUpdateResult:
    barcodes: File  # *barcodes.* (the first by name)
    lineages_topology: File  # *lineages.yml
    config: File  # *pathogen_config.yml
    lineages_meta: File | None  # *curated_lineages.json, when written
    db: Dir  # everything freyja update wrote


async def freyja_update(args: str = "") -> FreyjaUpdateResult:
    """The current barcodes and lineage data (for ``--pathogen`` in ``args``, SARS-CoV-2 by default)."""
    db = await freyja_update_cmd(args=args)
    files = sorted(await db.list_files(), key=lambda f: f.name)

    def first(suffix: str) -> File | None:
        return next((f for f in files if f.name.endswith(suffix)), None)

    barcodes = next((f for f in files if "barcodes." in f.name), None)
    topology, config = first("lineages.yml"), first("pathogen_config.yml")
    for name, found in (("*barcodes.*", barcodes), ("*lineages.yml", topology), ("*pathogen_config.yml", config)):
        if found is None:
            raise FileNotFoundError(f"freyja update wrote no {name}")
    return FreyjaUpdateResult(
        barcodes=barcodes,
        lineages_topology=topology,
        config=config,
        lineages_meta=first("curated_lineages.json"),
        db=db,
    )


@dataclass
class FreyjaVariantsResult:
    variants: File  # <prefix>.variants.tsv
    depths: File  # <prefix>.depth.tsv


async def freyja_variants(bam: File, fasta: File, prefix: str = "sample", args: str = "") -> FreyjaVariantsResult:
    """Variant calls and depths of the sorted ``bam`` against the reference ``fasta``."""
    results = await freyja_variants_cmd(bam=bam, fasta=fasta, prefix=prefix, args=args)
    return FreyjaVariantsResult(
        variants=await pick_file(results, f"{prefix}.variants.tsv", "freyja variants"),
        depths=await pick_file(results, f"{prefix}.depth.tsv", "freyja variants"),
    )


async def freyja_demix(
    variants: File,
    depths: File,
    barcodes: File,
    lineages_meta: File | None = None,
    lineages_topology: File | None = None,
    prefix: str = "sample",
    args: str = "",
) -> File:
    """The sample's lineage abundances: ``<prefix>.tsv``."""
    results = await freyja_demix_cmd(
        variants=variants,
        depths=depths,
        barcodes=barcodes,
        lineages_meta=lineages_meta,
        lineages_topology=lineages_topology,
        prefix=prefix,
        args=args,
    )
    return await pick_file(results, f"{prefix}.tsv", "freyja demix")


@dataclass
class FreyjaBootResult:
    lineages: File  # <prefix>_lineages.csv
    summarized: File  # <prefix>_summarized.csv


async def freyja_boot(
    variants: File,
    depths: File,
    repeats: int,
    barcodes: File,
    lineages_meta: File | None = None,
    lineages_topology: File | None = None,
    prefix: str = "sample",
    args: str = "",
) -> FreyjaBootResult:
    """Bootstrapped confidence intervals of the sample's lineage abundances."""
    results = await freyja_boot_cmd(
        variants=variants,
        depths=depths,
        repeats=repeats,
        barcodes=barcodes,
        lineages_meta=lineages_meta,
        lineages_topology=lineages_topology,
        prefix=prefix,
        args=args,
    )
    return FreyjaBootResult(
        lineages=await pick_file(results, f"{prefix}_lineages.csv", "freyja boot"),
        summarized=await pick_file(results, f"{prefix}_summarized.csv", "freyja boot"),
    )
