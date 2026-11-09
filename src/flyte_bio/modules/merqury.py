"""merqury — reference-free assembly evaluation from read k-mers.

Exposes :func:`merqury`, which evaluates one or two assemblies against a
read k-mer database (from :mod:`flyte_bio.modules.meryl`), optionally with
the parents' databases (trio mode: hap-mers, phasing blob plots). As
upstream, the outputs keep ``merqury.sh``'s names, with ``<prefix>`` as the
output name:

- ``<assembly>_only.bed`` / ``.wig``: assembly k-mers missing from the reads;
- ``<prefix>.completeness.stats``, ``<prefix>.dist_only.hist``;
- copy-number spectra ``*.spectra-cn.hist`` and assembly spectra
  ``<prefix>.spectra-asm.hist``, each with ``.fl.png``, ``.ln.png`` and
  ``.st.png`` plots;
- QV: ``<prefix>.qv`` and per-assembly ``<prefix>.<assembly>.qv``;
- ``<read db>.hist.ploidy`` and, in trio mode, ``<prefix>.hapmers.blob.png``.

Databases are ``<name>.meryl.tar.gz`` tarballs, unpacked under their names
(which the outputs carry). meryl's threads (``OMP_NUM_THREADS``) follow the
task's CPUs, as upstream.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

MERQURY_IMAGE = "quay.io/biocontainers/merqury:1.3--hdfd78af_4"
DEFAULT_RESOURCES = flyte.Resources(cpu=2, memory="12Gi")

merqury_cmd = shell.create(
    name="merqury",
    image=MERQURY_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "meryl_db": File,
        "maternal_db": File | None,
        "paternal_db": File | None,
        "assembly": File,
        "assembly_2": File | None,
        "prefix": str,
    },
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        DB=({{inputs.meryl_db}}); ASM=({{inputs.assembly}})
        shopt -s nullglob; MAT=({{inputs.maternal_db}}); PAT=({{inputs.paternal_db}}); ASM2=({{inputs.assembly_2}})
        shopt -u nullglob
        P={{inputs.prefix}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"

        # The container's entrypoint sets merqury's environment; source it, as upstream.
        if [ -f "/usr/local/env-activate.sh" ]; then
            set +u
            . "/usr/local/env-activate.sh"
        fi
        export OMP_NUM_THREADS=$TASK_CPUS

        DBS=()
        for T in "${{DB[@]}}" "${{MAT[@]}}" "${{PAT[@]}}"; do
            tar -xzf "$T"
            DBS+=("$(basename "$T" .tar.gz)")
        done
        ASMS=()
        for A in "${{ASM[@]}}" "${{ASM2[@]}}"; do
            ln -s "$A" "$(basename "$A")"
            ASMS+=("$(basename "$A")")
        done

        merqury.sh "${{DBS[@]}}" "${{ASMS[@]}}" "$P" > merqury.log 2>&1 || {{ tail -n 30 merqury.log >&2; exit 1; }}
        cat merqury.log >&2

        shopt -s nullglob
        for F in ./*_only.bed ./*_only.wig ./*.completeness.stats ./*.dist_only.hist ./*.spectra-cn.* \
                ./*.spectra-asm.* "$P.qv" "$P".*.qv ./*.hist.ploidy ./*.hapmers.blob.png; do
            cp -L "$F" "$OUT"/
        done
    """,
)

env = flyte.TaskEnvironment.from_task("merqury", merqury_cmd.as_task())


@dataclass
class MerquryResult:
    assembly_only_kmers_bed: list[File]  # <assembly>_only.bed, per assembly
    assembly_only_kmers_wig: list[File]  # <assembly>_only.wig, per assembly
    stats: File  # <prefix>.completeness.stats
    dist_hist: File  # <prefix>.dist_only.hist
    spectra_cn_hist: list[File]
    spectra_asm_hist: File  # <prefix>.spectra-asm.hist
    assembly_qv: File  # <prefix>.qv
    scaffold_qv: list[File]  # <prefix>.<assembly>.qv
    read_ploidy: File  # <read db>.hist.ploidy
    plots: list[File]  # every *.png
    results: Dir  # everything above


async def merqury(
    meryl_db: File,
    assembly: File,
    assembly_2: File | None = None,
    maternal_db: File | None = None,
    paternal_db: File | None = None,
    prefix: str = "sample",
) -> MerquryResult:
    """Evaluate ``assembly`` (and ``assembly_2``) against the reads' k-mers; trio mode with both parents."""
    if (maternal_db is None) != (paternal_db is None):
        raise ValueError("merqury trio mode needs both the maternal and the paternal database")
    results = await merqury_cmd(
        meryl_db=meryl_db,
        maternal_db=maternal_db,
        paternal_db=paternal_db,
        assembly=assembly,
        assembly_2=assembly_2,
        prefix=prefix,
    )
    files = {f.name: f async for f in results.walk()}

    def all_ending(suffix: str) -> list[File]:
        return [files[n] for n in sorted(files) if n.endswith(suffix)]

    def named(name: str) -> File:
        if name not in files:
            raise FileNotFoundError(f"merqury wrote no {name}")
        return files[name]

    ploidy = all_ending(".hist.ploidy")
    if len(ploidy) != 1:
        raise FileNotFoundError(f"merqury wrote {len(ploidy)} *.hist.ploidy files, expected one")
    return MerquryResult(
        assembly_only_kmers_bed=all_ending("_only.bed"),
        assembly_only_kmers_wig=all_ending("_only.wig"),
        stats=named(f"{prefix}.completeness.stats"),
        dist_hist=named(f"{prefix}.dist_only.hist"),
        spectra_cn_hist=all_ending(".spectra-cn.hist"),
        spectra_asm_hist=named(f"{prefix}.spectra-asm.hist"),
        assembly_qv=named(f"{prefix}.qv"),
        scaffold_qv=[
            files[n] for n in sorted(files) if n.startswith(f"{prefix}.") and n.endswith(".qv") and n != f"{prefix}.qv"
        ],
        read_ploidy=ploidy[0],
        plots=all_ending(".png"),
        results=results,
    )
