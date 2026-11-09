"""vireo — assign droplets to donors from their allele counts with Vireo.

Exposes :func:`vireo`, which demultiplexes ``n_donor`` donors from per-cell
allele counts (cellsnp-lite's cell VCF, or VarTrix's matrices) against the
donors' genotypes. Outputs keep upstream's names:

- ``<prefix>_summary.tsv``: droplets per donor;
- ``<prefix>_donor_ids.tsv``: each droplet's donor (or ``doublet`` /
  ``unassigned``);
- ``<prefix>_prob_singlet.tsv.gz``, ``<prefix>_prob_doublet.tsv.gz``: the
  assignment probabilities;
- when Vireo writes donor genotypes (e.g. ``--forceLearnGT``),
  ``<prefix>_GT_donors.vireo.vcf.gz`` and, from ``GTbarcode``,
  ``<prefix>_filtered_variants.tsv``.

Vireo gets the task's CPUs (``-p``); ``args`` go to it too (e.g.
``--forceLearnGT``). As upstream, ``--randSeed 42`` is added unless ``args``
sets a seed, and ``GTbarcode`` gets the same seed.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

VIREO_IMAGE = "community.wave.seqera.io/library/vireosnp:0.5.8--376a9219dbddd077"
# Upstream: process_medium.
DEFAULT_RESOURCES = flyte.Resources(cpu=4, memory="8Gi")

vireo_cmd = shell.create(
    name="vireo",
    image=VIREO_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "cell_data": File | None,
        "vartrix_data": list[File] | None,
        "vartrix_order": str,
        "n_donor": int,
        "donor_file": File,
        "prefix": str,
        "args": str,
    },
    defaults={"vartrix_order": "", "args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        shopt -s nullglob; CELLS=({{inputs.cell_data}}); VARTRIX=({{inputs.vartrix_data}}); shopt -u nullglob
        DONORS=({{inputs.donor_file}})
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        if [[ ! "$ARGS" =~ --randSeed ]]; then ARGS="$ARGS --randSeed 42"; fi
        SEED=""
        if [[ "$ARGS" =~ (--randSeed[[:space:]]+[0-9]+) ]]; then SEED="${{BASH_REMATCH[1]}}"; fi
        if [ ${{#CELLS[@]}} -gt 0 ]; then
            INPUT=(-c "${{CELLS[0]}}")
        else
            # Staged files glob in name order; put them back in the caller's order.
            VARTRIX_FILES=(); ORDER={{inputs.vartrix_order}}
            for NAME in $ORDER; do
                for F in "${{VARTRIX[@]}}"; do
                    if [ "$(basename "$F")" = "$NAME" ]; then VARTRIX_FILES+=("$F"); fi
                done
            done
            INPUT=(--vartrixData "$(IFS=,; echo "${{VARTRIX_FILES[*]}}")")
        fi
        vireo \
            "${{INPUT[@]}}" \
            -N {{inputs.n_donor}} \
            -d "${{DONORS[0]}}" \
            -p $TASK_CPUS \
            -o . \
            $ARGS || exit 1
        mv summary.tsv "$OUT/${{P}}_summary.tsv"
        mv donor_ids.tsv "$OUT/${{P}}_donor_ids.tsv"
        mv prob_singlet.tsv.gz "$OUT/${{P}}_prob_singlet.tsv.gz"
        mv prob_doublet.tsv.gz "$OUT/${{P}}_prob_doublet.tsv.gz"
        if [[ -f GT_donors.vireo.vcf.gz ]]; then
            mv GT_donors.vireo.vcf.gz "${{P}}_GT_donors.vireo.vcf.gz"
            GTbarcode -i "${{P}}_GT_donors.vireo.vcf.gz" -o "./${{P}}_filtered_variants.tsv" $SEED || exit 1
            mv "${{P}}_GT_donors.vireo.vcf.gz" "${{P}}_filtered_variants.tsv" "$OUT"/
        fi
    """,
)

env = flyte.TaskEnvironment.from_task("vireo", vireo_cmd.as_task())


@dataclass
class VireoResult:
    summary: File  # <prefix>_summary.tsv
    donor_ids: File  # <prefix>_donor_ids.tsv
    prob_singlets: File  # <prefix>_prob_singlet.tsv.gz
    prob_doublets: File  # <prefix>_prob_doublet.tsv.gz
    genotype_vcf: File | None  # <prefix>_GT_donors.vireo.vcf.gz
    filtered_variants: File | None  # <prefix>_filtered_variants.tsv


async def vireo(
    n_donor: int,
    donor_file: File,
    prefix: str,
    cell_data: File | None = None,
    vartrix_data: list[File] | None = None,
    args: str = "",
) -> VireoResult:
    """Vireo donor assignments from ``cell_data`` (a cellsnp-lite cell VCF) or ``vartrix_data``.

    ``vartrix_data`` is VarTrix's alt and ref matrices and barcodes (and optionally the SNP VCF), in
    Vireo's ``--vartrixData`` order; their file names must be distinct.
    """
    names = [f.name for f in vartrix_data or []]
    if len(set(names)) != len(names):
        raise ValueError(f"vireo needs distinct vartrix_data file names, got {names}")
    if (cell_data is None) == (not vartrix_data):
        raise ValueError("vireo needs exactly one of cell_data and vartrix_data")
    results = await vireo_cmd(
        cell_data=cell_data,
        vartrix_data=vartrix_data,
        vartrix_order=" ".join(f.name for f in vartrix_data or []),
        n_donor=n_donor,
        donor_file=donor_file,
        prefix=prefix,
        args=args,
    )

    async def pick(name: str) -> File:
        found = await results.get_file(name)
        if found is None:
            raise FileNotFoundError(f"vireo wrote no {name}")
        return found

    return VireoResult(
        summary=await pick(f"{prefix}_summary.tsv"),
        donor_ids=await pick(f"{prefix}_donor_ids.tsv"),
        prob_singlets=await pick(f"{prefix}_prob_singlet.tsv.gz"),
        prob_doublets=await pick(f"{prefix}_prob_doublet.tsv.gz"),
        genotype_vcf=await results.get_file(f"{prefix}_GT_donors.vireo.vcf.gz"),
        filtered_variants=await results.get_file(f"{prefix}_filtered_variants.tsv"),
    )
