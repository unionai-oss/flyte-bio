"""mtmalign — multiple protein structure alignment with mTM-align.

Exposes :func:`mtmalign_align`, which aligns protein structures (PDB files,
may be gzipped) and returns the sequence alignment, ``<prefix>.aln``, and the
superposed structures, ``<prefix>.pdb``; both gzipped with pigz with
``compress``. As upstream, the ``.pdb`` extensions are removed from the
alignment's sequence names.

mTM-align takes no thread count and can't write its alignment elsewhere than
``mTM_result/result.fasta``; it runs in a scratch directory holding the
structures. Upstream decompresses gzipped structures in place; here they are
decompressed to copies next to them, which is the same input.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

# Upstream's Docker image: mTM-align 20220104 and pigz.
MTMALIGN_IMAGE = (
    "quay.io/biocontainers/mulled-v2-5bcf71dc66dac33d8e003c5e78043b80f5c7f269:"
    "8f0e486d46f7ab38892c1a8f78d2894a549d03b5-0"
)
# Upstream labels this process_medium; mTM-align itself is single-threaded.
DEFAULT_RESOURCES = flyte.Resources(cpu=2, memory="8Gi")

mtmalign_align_cmd = shell.create(
    name="mtmalign_align",
    image=MTMALIGN_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"structures": list[File], "prefix": str, "compress": bool},
    defaults={"compress": False},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        STRUCTURES=({{inputs.structures}})
        P={{inputs.prefix}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"

        # decompress input files if required
        for F in "${{STRUCTURES[@]}}"; do
            NAME=$(basename "$F")
            if [[ "$NAME" == *.pdb.gz ]]; then
                unpigz -c "$F" > "${{NAME%.gz}}"
            else
                ln -s "$F" "$NAME"
            fi
        done

        # construct input file for mtmalign
        ls *.pdb | sed s/\ /\n/ > input_list.txt

        mtm-align -i input_list.txt -o "$P.pdb" || exit 1
        # -o does not affect the fasta naming, so move it to the new name
        mv ./mTM_result/result.fasta "$P.aln"
        mv ./mTM_result/"$P.pdb" "$P.pdb"

        # Remove ".pdb" from the ids in the alignment file
        sed -i 's/\.pdb//g' "$P.aln"

        # compress both output files
        if [ {{inputs.compress}} = true ]; then
            pigz -p $TASK_CPUS "$P.aln" "$P.pdb"
        fi
        mv "$P".aln* "$P".pdb* "$OUT"/
    """,
)

env = flyte.TaskEnvironment.from_task("mtmalign", mtmalign_align_cmd.as_task())


@dataclass
class MtmalignResult:
    alignment: File  # <prefix>.aln(.gz): the sequence alignment
    structure: File  # <prefix>.pdb(.gz): the superposed structures


async def mtmalign_align(structures: list[File], prefix: str = "sample", compress: bool = False) -> MtmalignResult:
    """Align ``structures``: the alignment and the superposed structures."""
    results = await mtmalign_align_cmd(structures=structures, prefix=prefix, compress=compress)
    suffix = ".gz" if compress else ""

    async def pick(name: str) -> File:
        found = await results.get_file(name)
        if found is None:
            raise FileNotFoundError(f"mtm-align wrote no {name}")
        return found

    return MtmalignResult(alignment=await pick(f"{prefix}.aln{suffix}"), structure=await pick(f"{prefix}.pdb{suffix}"))
