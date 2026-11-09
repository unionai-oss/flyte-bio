"""ragtag — RagTag, reference-guided scaffolding and patching of assemblies.

Exposes:

- :func:`ragtag_scaffold` — order and orient an assembly's contigs against a
  reference (``ragtag.py scaffold``). Outputs, renamed from
  ``ragtag.scaffold.*`` to ``<prefix>.*`` as upstream: the scaffolds
  ``<prefix>.fasta.gz``, ``<prefix>.agp``, ``<prefix>.stats``,
  ``<prefix>.confidence.txt`` and the alignments ``<prefix>.asm.paf.gz``
  (minimap2, with ``<prefix>.asm.paf.log``) or ``<prefix>.asm.delta.gz``
  (nucmer, with ``<prefix>.asm.delta.log``). Contigs to exclude (``-e``),
  skip (``-j``) or hard-skip (``-J``) are given as files of names.
- :func:`ragtag_patch` — fill gaps and join contigs of a target assembly
  with sequences from a query assembly (``ragtag.py patch``). Outputs:
  ``<prefix>.patch.fasta``, ``<prefix>.patch.agp``, ``<prefix>.comps.fasta``
  (the patch components), the target split at gaps ``<prefix>.ctg.fasta`` /
  ``<prefix>.ctg.agp``, when written the renamed query
  ``<prefix>.rename.fasta`` / ``.agp``, the alignments
  ``<prefix>_ragtag.patch.asm.*`` and RagTag's log ``<prefix>.patch.err``.

Inputs may be plain or gzipped FASTA. ``-t`` follows the task's CPUs; pass
other options (``--aligner``, ``-u``, …) in ``args``.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

RAGTAG_SCAFFOLD_IMAGE = "community.wave.seqera.io/library/ragtag_gzip:5de7cc3ac89b3ac4"
RAGTAG_IMAGE = "community.wave.seqera.io/library/ragtag:2.1.0--6125126fe9d553ab"
DEFAULT_RESOURCES = flyte.Resources(cpu=4, memory="16Gi")

# Write <name>.fa from a plain or gzipped FASTA.
LINK_FASTA = r"""
        link_fasta() {
            if [[ "$1" == *.gz ]]; then zcat "$1" > "$2.fa"; else ln -s "$1" "$2.fa"; fi
        }
"""

ragtag_scaffold_cmd = shell.create(
    name="ragtag_scaffold",
    image=RAGTAG_SCAFFOLD_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "assembly": File,
        "reference": File,
        "exclude": File | None,
        "skip": File | None,
        "hard_skip": File | None,
        "prefix": str,
        "args": str,
    },
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        ASM=({{inputs.assembly}}); REF=({{inputs.reference}})
        shopt -s nullglob
        EX=({{inputs.exclude}}); SKIP=({{inputs.skip}}); HARD=({{inputs.hard_skip}})
        shopt -u nullglob
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        {LINK_FASTA}
        link_fasta "${{ASM[0]}}" assembly
        link_fasta "${{REF[0]}}" reference
        EXTRA=()
        if [ ${{#EX[@]}} -gt 0 ]; then EXTRA+=(-e "${{EX[0]}}"); fi
        if [ ${{#SKIP[@]}} -gt 0 ]; then EXTRA+=(-j "${{SKIP[0]}}"); fi
        if [ ${{#HARD[@]}} -gt 0 ]; then EXTRA+=(-J "${{HARD[0]}}"); fi
        ragtag.py scaffold reference.fa assembly.fa -o . -t $TASK_CPUS "${{EXTRA[@]}}" $ARGS \
            2> "$P.stderr.log" > "$P.stdout.log" || {{ tail -n 20 "$P.stderr.log" >&2; exit 1; }}
        cat "$P.stderr.log" >&2
        for F in ragtag.scaffold.*; do mv "$F" "${{F/ragtag.scaffold/$P}}"; done
        gzip ./*.fasta 2>/dev/null || true
        gzip ./*.delta 2>/dev/null || true
        gzip ./*.paf 2>/dev/null || true
        shopt -s nullglob
        mv ./*.fasta.gz ./*.agp ./*.stats ./*.paf.gz ./*.paf.log ./*.delta.gz ./*.delta.log ./*.confidence.txt "$OUT"/
    """,
)

ragtag_patch_cmd = shell.create(
    name="ragtag_patch",
    image=RAGTAG_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"target": File, "query": File, "exclude": File | None, "skip": File | None, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        TGT=({{inputs.target}}); QRY=({{inputs.query}})
        shopt -s nullglob; EX=({{inputs.exclude}}); SKIP=({{inputs.skip}}); shopt -u nullglob
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        {LINK_FASTA}
        link_fasta "${{TGT[0]}}" target
        link_fasta "${{QRY[0]}}" query
        EXTRA=()
        if [ ${{#EX[@]}} -gt 0 ]; then EXTRA+=(-e "${{EX[0]}}"); fi
        if [ ${{#SKIP[@]}} -gt 0 ]; then EXTRA+=(-j "${{SKIP[0]}}"); fi
        ragtag.py patch target.fa query.fa -o "$P" -t $TASK_CPUS "${{EXTRA[@]}}" $ARGS \
            2> "$P.stderr.log" > "$P.stdout.log" || {{
            tail -n 20 "$P.stderr.log" "$P/ragtag.patch.err" >&2; exit 1;
        }}
        cat "$P/ragtag.patch.err" >&2
        mv "$P/ragtag.patch.agp" "$OUT/$P.patch.agp"
        mv "$P/ragtag.patch.fasta" "$OUT/$P.patch.fasta"
        mv "$P/ragtag.patch.comps.fasta" "$OUT/$P.comps.fasta"
        mv "$P/ragtag.patch.ctg.agp" "$OUT/$P.ctg.agp"
        mv "$P/ragtag.patch.ctg.fasta" "$OUT/$P.ctg.fasta"
        if [ -f "$P/ragtag.patch.rename.agp" ]; then mv "$P/ragtag.patch.rename.agp" "$OUT/$P.rename.agp"; fi
        if [ -f "$P/ragtag.patch.rename.fasta" ]; then mv "$P/ragtag.patch.rename.fasta" "$OUT/$P.rename.fasta"; fi
        mv "$P/ragtag.patch.err" "$OUT/$P.patch.err"
        shopt -s nullglob
        for F in "$P"/ragtag.patch.asm.*; do mv "$F" "$OUT/${{P}}_$(basename "$F")"; done
    """,
)

ragtag_scaffold_env = flyte.TaskEnvironment.from_task("ragtag_scaffold", ragtag_scaffold_cmd.as_task())
ragtag_patch_env = flyte.TaskEnvironment.from_task("ragtag_patch", ragtag_patch_cmd.as_task())
env = flyte.TaskEnvironment(name="ragtag", depends_on=[ragtag_scaffold_env, ragtag_patch_env])


@dataclass
class RagtagScaffoldResult:
    corrected_assembly: File  # <prefix>.fasta.gz
    corrected_agp: File  # <prefix>.agp
    corrected_stats: File  # <prefix>.stats
    confidence_txt: File  # <prefix>.confidence.txt
    ragtag_paf: File | None  # <prefix>.asm.paf.gz (minimap2)
    ragtag_paf_log: File | None  # <prefix>.asm.paf.log
    ragtag_delta: File | None  # <prefix>.asm.delta.gz (nucmer)
    ragtag_delta_log: File | None  # <prefix>.asm.delta.log
    results: Dir  # everything above


@dataclass
class RagtagPatchResult:
    patch_fasta: File  # <prefix>.patch.fasta
    patch_agp: File  # <prefix>.patch.agp
    patch_components_fasta: File  # <prefix>.comps.fasta
    target_splits_agp: File  # <prefix>.ctg.agp
    target_splits_fasta: File  # <prefix>.ctg.fasta
    qry_rename_agp: File | None  # <prefix>.rename.agp
    qry_rename_fasta: File | None  # <prefix>.rename.fasta
    assembly_alignments: list[File]  # <prefix>_ragtag.patch.asm.*
    stderr: File  # <prefix>.patch.err
    results: Dir  # everything above


async def named(results: Dir, name: str, tool: str) -> File:
    found = await results.get_file(name)
    if found is None:
        raise FileNotFoundError(f"ragtag {tool} wrote no {name}")
    return found


async def ragtag_scaffold(
    assembly: File,
    reference: File,
    exclude: File | None = None,
    skip: File | None = None,
    hard_skip: File | None = None,
    prefix: str = "sample",
    args: str = "",
) -> RagtagScaffoldResult:
    """Scaffold ``assembly`` against ``reference``."""
    results = await ragtag_scaffold_cmd(
        assembly=assembly,
        reference=reference,
        exclude=exclude,
        skip=skip,
        hard_skip=hard_skip,
        prefix=prefix,
        args=args,
    )
    return RagtagScaffoldResult(
        corrected_assembly=await named(results, f"{prefix}.fasta.gz", "scaffold"),
        corrected_agp=await named(results, f"{prefix}.agp", "scaffold"),
        corrected_stats=await named(results, f"{prefix}.stats", "scaffold"),
        confidence_txt=await named(results, f"{prefix}.confidence.txt", "scaffold"),
        ragtag_paf=await results.get_file(f"{prefix}.asm.paf.gz"),
        ragtag_paf_log=await results.get_file(f"{prefix}.asm.paf.log"),
        ragtag_delta=await results.get_file(f"{prefix}.asm.delta.gz"),
        ragtag_delta_log=await results.get_file(f"{prefix}.asm.delta.log"),
        results=results,
    )


async def ragtag_patch(
    target: File,
    query: File,
    exclude: File | None = None,
    skip: File | None = None,
    prefix: str = "sample",
    args: str = "",
) -> RagtagPatchResult:
    """Patch ``target`` with sequences from ``query``."""
    results = await ragtag_patch_cmd(target=target, query=query, exclude=exclude, skip=skip, prefix=prefix, args=args)
    alignments = sorted(
        [f async for f in results.walk() if f.name.startswith(f"{prefix}_ragtag.patch.asm.")], key=lambda f: f.name
    )
    return RagtagPatchResult(
        patch_fasta=await named(results, f"{prefix}.patch.fasta", "patch"),
        patch_agp=await named(results, f"{prefix}.patch.agp", "patch"),
        patch_components_fasta=await named(results, f"{prefix}.comps.fasta", "patch"),
        target_splits_agp=await named(results, f"{prefix}.ctg.agp", "patch"),
        target_splits_fasta=await named(results, f"{prefix}.ctg.fasta", "patch"),
        qry_rename_agp=await results.get_file(f"{prefix}.rename.agp"),
        qry_rename_fasta=await results.get_file(f"{prefix}.rename.fasta"),
        assembly_alignments=alignments,
        stderr=await named(results, f"{prefix}.patch.err", "patch"),
        results=results,
    )
