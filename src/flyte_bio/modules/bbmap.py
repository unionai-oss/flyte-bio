"""bbmap — BBSplit: bin reads by which reference they map to best.

Exposes :func:`bbsplit_index`, which builds a BBSplit index from a primary
FASTA plus other references (e.g. contaminant genomes), and
:func:`bbsplit`, which splits one sample's reads against that index and
returns the reads assigned to the primary reference.

The index-reuse handling is upstream's: BBSplit records each reference's
source path and modification time in ``ref/genome/*/summary.txt`` and
rebuilds the index if they don't match. Staging an index changes both, so
before splitting the index is laid out in a writable copy whose summaries
point at the staged files with their current timestamps. Both steps run in a
scratch directory under upstream's names (``bbsplit_build`` →
``bbsplit_index``, ``input_index`` → ``index_writable``) so the path
rewrites behave as upstream intends.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

BBMAP_IMAGE = "community.wave.seqera.io/library/bbmap_pigz:07416fe99b090fa9"

BBSPLIT_CPUS = 8
BBSPLIT_MEMORY_MB = 16384
# As upstream: give the JVM 80% of the task's memory.
BBSPLIT_XMX_MB = int(BBSPLIT_MEMORY_MB * 0.8)
DEFAULT_RESOURCES = flyte.Resources(cpu=BBSPLIT_CPUS, memory=f"{BBSPLIT_MEMORY_MB}Mi")

# Other references are passed as explicit `name=filename` pairs: a list[File]
# is staged under original basenames and globbed alphabetically, so file
# order can't carry the names.
bbsplit_index_cmd = shell.create(
    name="bbsplit_index",
    image=BBMAP_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"primary_ref": File, "other_refs": list[File], "other_ref_names": str, "args": str},
    defaults={"args": ""},
    outputs={"index": Dir},
    script=rf"""
        PRIMARY=({{inputs.primary_ref}})
        O=({{inputs.other_refs}})
        ODIR=$(dirname "${{O[0]}}")
        REFS=()
        IFS=',' read -ra PAIRS <<< {{inputs.other_ref_names}}
        for pair in "${{PAIRS[@]}}"; do REFS+=("ref_${{pair%%=*}}=$ODIR/${{pair#*=}}"); done
        ARGS={{inputs.args}}
        W=$(mktemp -d); cd "$W"
        bbsplit.sh \
            -Xmx{BBSPLIT_XMX_MB}M \
            ref_primary="${{PRIMARY[0]}}" "${{REFS[@]}}" path=bbsplit_build \
            threads={BBSPLIT_CPUS} \
            $ARGS

        # Summary files will have an absolute path that will make the index
        # impossible to use in other processes - fix paths and rename atomically
        if [ -d bbsplit_build/ref/genome ]; then
            find bbsplit_build/ref/genome -name summary.txt | while read -r summary_file; do
                sed "s|^source.*|source\t$(grep '^source' "$summary_file" | cut -f2- -d$'\t' | sed 's|.*/bbsplit_build|bbsplit_index|')|" "$summary_file" > ${{summary_file}}.tmp && mv ${{summary_file}}.tmp ${{summary_file}}
            done
            mv bbsplit_build bbsplit_index
        fi
        cp -r bbsplit_index/. {{outputs.index}}/
    """,
)


bbsplit_cmd = shell.create(
    name="bbsplit",
    image=BBMAP_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"reads_1": File, "reads_2": File | None, "index": Dir, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        R1=({{inputs.reads_1}})
        shopt -s nullglob; R2=({{inputs.reads_2}}); shopt -u nullglob
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        W=$(mktemp -d); cd "$W"
        ln -s {{inputs.index}} input_index

        # If using a pre-built index, create writable structure: symlink all files except
        # summary.txt (which we copy to modify). When we stage in the index files the time
        # stamps get disturbed, which bbsplit doesn't like. Fix the time stamps in summaries.
        find input_index/ref -type f | while read -r f; do
            target="index_writable/${{f#input_index/}}"
            mkdir -p "$(dirname "$target")"
            [[ $(basename "$f") == "summary.txt" ]] && cp "$f" "$target" || ln -s "$(realpath "$f")" "$target"
        done
        find index_writable/ref/genome -name summary.txt | while read -r summary_file; do
            src=$(grep '^source' "$summary_file" | cut -f2- -d$'\t' | sed 's|.*/ref/|index_writable/ref/|')
            mod=$(echo "System.out.println(java.nio.file.Files.getLastModifiedTime(java.nio.file.Paths.get(\"$src\")).toMillis());" | jshell -J-Djdk.lang.Process.launchMechanism=vfork - 2>/dev/null | grep -oE '^[0-9]{{12,14}}$')
            sed -e 's|bbsplit_index/ref|index_writable/ref|' -e "s|^last modified.*|last modified\t$mod|" "$summary_file" > ${{summary_file}}.tmp && mv ${{summary_file}}.tmp ${{summary_file}}
        done

        if [ ${{#R2[@]}} -eq 0 ]; then
            IN=(in="${{R1[0]}}"); OUT="basename=${{P}}_%.fastq.gz"
        else
            IN=(in="${{R1[0]}}" in2="${{R2[0]}}"); OUT="basename=${{P}}_%_#.fastq.gz"
        fi
        bbsplit.sh \
            -Xmx{BBSPLIT_XMX_MB}M \
            path=index_writable \
            threads={BBSPLIT_CPUS} \
            "${{IN[@]}}" \
            "$OUT" \
            refstats="$P.stats.txt" \
            $ARGS 2>| >(tee "$P.log" >&2)
        sleep 1  # let the tee'd log finish writing
        mv ./*.fastq.gz ./*.txt ./*.log {{outputs.results}}/
    """,
)


env = flyte.TaskEnvironment.from_task(
    "bbmap",
    bbsplit_index_cmd.as_task(),
    bbsplit_cmd.as_task(),
)


def basename(f: File) -> str:
    return f.path.rstrip("/").rsplit("/", 1)[-1]


async def bbsplit_index(primary_ref: File, other_refs: dict[str, File], args: str = "") -> Dir:
    """Build a BBSplit index from ``primary_ref`` plus ``other_refs`` (name → FASTA)."""
    if not other_refs:
        raise ValueError("bbsplit_index needs at least one other reference")
    names = {name: basename(f) for name, f in other_refs.items()}
    if len(set(names.values())) != len(names):
        raise ValueError(f"other reference files must have distinct basenames: {names}")
    for name in names:
        if not name.replace("_", "").isalnum():
            raise ValueError(f"reference name {name!r} must be alphanumeric/underscore")
    return await bbsplit_index_cmd(
        primary_ref=primary_ref,
        other_refs=list(other_refs.values()),
        other_ref_names=",".join(f"{name}={fname}" for name, fname in names.items()),
        args=args,
    )


@dataclass
class BBSplitResult:
    reads_1: File  # reads assigned to the primary reference
    reads_2: File | None
    stats: File  # <prefix>.stats.txt
    log: File
    results: Dir  # everything BBSplit wrote, incl. reads binned to the other references


async def bbsplit(
    reads_1: File, reads_2: File | None, index: Dir, prefix: str = "bbsplit", args: str = ""
) -> BBSplitResult:
    """Split one sample's reads against ``index``; ``reads_2`` is None for single-end."""
    results = await bbsplit_cmd(
        reads_1=reads_1,
        reads_2=reads_2,
        index=index,
        prefix=prefix,
        args=args,
    )

    async def pick(name: str) -> File:
        found = await results.get_file(name)
        if found is None:
            raise FileNotFoundError(f"bbsplit wrote no {name}")
        return found

    stats, log = await pick(f"{prefix}.stats.txt"), await pick(f"{prefix}.log")
    if reads_2 is None:
        return BBSplitResult(await pick(f"{prefix}_primary.fastq.gz"), None, stats, log, results)
    return BBSplitResult(
        await pick(f"{prefix}_primary_1.fastq.gz"), await pick(f"{prefix}_primary_2.fastq.gz"), stats, log, results
    )
