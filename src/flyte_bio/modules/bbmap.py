"""bbmap — BBTools: BBSplit binning, BBDuk trimming/filtering, BBNorm normalization and FilterByName.

Exposes :func:`bbsplit_index`, which builds a BBSplit index from a primary
FASTA plus other references (e.g. contaminant genomes), and
:func:`bbsplit`, which splits one sample's reads against that index and
returns the reads assigned to the primary reference.

It also exposes :func:`bbduk`, which trims or filters one sample's reads
(quality trimming, adapter or contaminant k-mer matching against an optional
reference FASTA, as ``args`` choose). As upstream, it writes
``<prefix>.fastq.gz`` (single-end) or ``<prefix>_1.fastq.gz`` /
``<prefix>_2.fastq.gz``, and BBDuk's output in ``<prefix>.bbduk.log``.

And it exposes :func:`bbnorm`, which normalizes one sample's read depth
(``target=``/``min=`` in ``args``). As upstream, it writes
``<prefix>.fastq.gz`` (single-end: one or more FASTQs, read in turn, or an
interleaved one) or ``<prefix>_1.nm.fastq.gz`` / ``<prefix>_2.nm.fastq.gz``,
and BBNorm's output in ``<prefix>.bbnorm.log``.

And it exposes :func:`bbmap_filterbyname`, which drops (or, with
``include=t`` in ``args``, keeps only) the reads named in ``names`` and
converts them to ``output_format`` (``fastq.gz``, ``fasta``, ``fastq.bz2``…,
from the extension). As upstream, it writes ``<prefix>.<format>`` (single-end,
or with ``interleaved_output``) or ``<prefix>_1.<format>`` /
``<prefix>_2.<format>``, and FilterByName's standard output in
``<prefix>.log``. Paired reads come as two files, or interleaved in one.

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

from flyte_bio.limits import JAVA_HEAP_MB, TASK_LIMITS

BBMAP_IMAGE = "community.wave.seqera.io/library/bbmap_pigz:07416fe99b090fa9"

BBSPLIT_CPUS = 8
BBSPLIT_MEMORY_MB = 16384
DEFAULT_RESOURCES = flyte.Resources(cpu=BBSPLIT_CPUS, memory=f"{BBSPLIT_MEMORY_MB}Mi")
# Upstream's process_medium is 6 CPUs.
BBDUK_RESOURCES = flyte.Resources(cpu=6, memory="16Gi")
BBNORM_RESOURCES = flyte.Resources(cpu=6, memory="16Gi")
# Upstream labels FilterByName process_single (1 CPU, 6 GB).
FILTERBYNAME_RESOURCES = flyte.Resources(cpu=1, memory="6Gi")

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
        {TASK_LIMITS}
        PRIMARY=({{inputs.primary_ref}})
        O=({{inputs.other_refs}})
        ODIR=$(dirname "${{O[0]}}")
        REFS=()
        IFS=',' read -ra PAIRS <<< {{inputs.other_ref_names}}
        for pair in "${{PAIRS[@]}}"; do REFS+=("ref_${{pair%%=*}}=$ODIR/${{pair#*=}}"); done
        ARGS={{inputs.args}}
        W=$(mktemp -d); cd "$W"
        bbsplit.sh \
            -Xmx{JAVA_HEAP_MB}M \
            ref_primary="${{PRIMARY[0]}}" "${{REFS[@]}}" path=bbsplit_build \
            threads=$TASK_CPUS \
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
        {TASK_LIMITS}
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
            -Xmx{JAVA_HEAP_MB}M \
            path=index_writable \
            threads=$TASK_CPUS \
            "${{IN[@]}}" \
            "$OUT" \
            refstats="$P.stats.txt" \
            $ARGS 2>| >(tee "$P.log" >&2)
        sleep 1  # let the tee'd log finish writing
        mv ./*.fastq.gz ./*.txt ./*.log {{outputs.results}}/
    """,
)

bbduk_cmd = shell.create(
    name="bbduk",
    image=BBMAP_IMAGE,
    resources=BBDUK_RESOURCES,
    inputs={"reads_1": File, "reads_2": File | None, "contaminants": File | None, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        R1=({{inputs.reads_1}})
        shopt -s nullglob; R2=({{inputs.reads_2}}); CONTAMINANTS=({{inputs.contaminants}}); shopt -u nullglob
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        if [ ${{#R2[@]}} -eq 0 ]; then
            RAW=(in="${{R1[0]}}"); TRIMMED=(out="$OUT/$P.fastq.gz")
        else
            RAW=(in1="${{R1[0]}}" in2="${{R2[0]}}"); TRIMMED=(out1="$OUT/${{P}}_1.fastq.gz" out2="$OUT/${{P}}_2.fastq.gz")
        fi
        REF=(); [ ${{#CONTAMINANTS[@]}} -gt 0 ] && REF=(ref="${{CONTAMINANTS[0]}}")
        W=$(mktemp -d); cd "$W"
        bbduk.sh \
            -Xmx{JAVA_HEAP_MB}M \
            "${{RAW[@]}}" \
            "${{TRIMMED[@]}}" \
            threads=$TASK_CPUS \
            $ARGS \
            "${{REF[@]}}" \
            &> "$OUT/$P.bbduk.log"
    """,
)

# Single-end reads may be several files: as upstream, they're joined into one
# comma-separated `in=` (in staging order, i.e. by name). Multipass temp files
# go under the working directory, as upstream keeps them in the task's.
bbnorm_cmd = shell.create(
    name="bbnorm",
    image=BBMAP_IMAGE,
    resources=BBNORM_RESOURCES,
    inputs={"reads": list[File], "reads_2": File | None, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        READS=({{inputs.reads}})
        shopt -s nullglob; R2=({{inputs.reads_2}}); shopt -u nullglob
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        if [ ${{#R2[@]}} -eq 0 ]; then
            IN=(in="$(IFS=,; echo "${{READS[*]}}")"); NORMALIZED=(out="$OUT/$P.fastq.gz")
        else
            IN=(in="${{READS[0]}}" in2="${{R2[0]}}")
            NORMALIZED=(out1="$OUT/${{P}}_1.nm.fastq.gz" out2="$OUT/${{P}}_2.nm.fastq.gz")
        fi
        W=$(mktemp -d); cd "$W"
        mkdir -p tmp
        bbnorm.sh \
            "${{IN[@]}}" \
            "${{NORMALIZED[@]}}" \
            $ARGS \
            threads=$TASK_CPUS \
            -Xmx{JAVA_HEAP_MB}M \
            tmpdir=tmp \
            &> "$OUT/$P.bbnorm.log"
        rm -rf tmp
    """,
)


# As upstream: two mates are read with in= and in2=; one file with
# single_end false is read as interleaved pairs. The output is split into
# mates unless single_end or interleaved_output.
bbmap_filterbyname_cmd = shell.create(
    name="bbmap_filterbyname",
    image=BBMAP_IMAGE,
    resources=FILTERBYNAME_RESOURCES,
    inputs={
        "reads_1": File,
        "reads_2": File | None,
        "names": str,
        "output_format": str,
        "single_end": bool,
        "interleaved_output": bool,
        "prefix": str,
        "args": str,
    },
    defaults={"names": "", "interleaved_output": False, "args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        R1=({{inputs.reads_1}})
        shopt -s nullglob; R2=({{inputs.reads_2}}); shopt -u nullglob
        P={{inputs.prefix}}; FORMAT={{inputs.output_format}}; NAMES={{inputs.names}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        IN=(in="${{R1[0]}}"); if [ ${{#R2[@]}} -gt 0 ]; then IN+=(in2="${{R2[0]}}"); fi
        if [ {{inputs.single_end}} = true ] || [ {{inputs.interleaved_output}} = true ]; then
            FILTERED=(out="$OUT/$P.$FORMAT")
        else
            FILTERED=(out1="$OUT/${{P}}_1.$FORMAT" out2="$OUT/${{P}}_2.$FORMAT")
        fi
        SELECTED=(); if [ -n "$NAMES" ]; then SELECTED=(names="$NAMES"); fi
        W=$(mktemp -d); cd "$W"
        filterbyname.sh \
            -Xmx{JAVA_HEAP_MB}M \
            "${{IN[@]}}" \
            "${{FILTERED[@]}}" \
            "${{SELECTED[@]}}" \
            $ARGS \
            | tee "$OUT/$P.log"
    """,
)

env = flyte.TaskEnvironment.from_task(
    "bbmap",
    bbsplit_index_cmd.as_task(),
    bbsplit_cmd.as_task(),
    bbduk_cmd.as_task(),
    bbnorm_cmd.as_task(),
    bbmap_filterbyname_cmd.as_task(),
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


@dataclass
class BBDukResult:
    reads: list[File]  # <prefix>.fastq.gz, or <prefix>_1.fastq.gz and <prefix>_2.fastq.gz
    log: File  # <prefix>.bbduk.log


async def bbduk(
    reads_1: File, reads_2: File | None = None, contaminants: File | None = None, prefix: str = "bbduk", args: str = ""
) -> BBDukResult:
    """Trim or filter one sample's reads with BBDuk; ``reads_2`` is None for single-end."""
    results = await bbduk_cmd(
        reads_1=reads_1,
        reads_2=reads_2,
        contaminants=contaminants,
        prefix=prefix,
        args=args,
    )

    async def pick(name: str) -> File:
        found = await results.get_file(name)
        if found is None:
            raise FileNotFoundError(f"bbduk wrote no {name}")
        return found

    names = [f"{prefix}.fastq.gz"] if reads_2 is None else [f"{prefix}_1.fastq.gz", f"{prefix}_2.fastq.gz"]
    return BBDukResult(reads=[await pick(name) for name in names], log=await pick(f"{prefix}.bbduk.log"))


@dataclass
class BBNormResult:
    reads: list[File]  # <prefix>.fastq.gz, or <prefix>_1.nm.fastq.gz and <prefix>_2.nm.fastq.gz
    log: File  # <prefix>.bbnorm.log


async def bbnorm(
    reads_1: File | list[File], reads_2: File | None = None, prefix: str = "bbnorm", args: str = ""
) -> BBNormResult:
    """Normalize one sample's read depth with BBNorm.

    Single-end (``reads_2`` None): ``reads_1`` is one FASTQ (interleaved pairs are detected) or a list
    of FASTQs read in turn, by name. Paired-end: ``reads_1`` and ``reads_2`` are one FASTQ each.
    """
    reads = reads_1 if isinstance(reads_1, list) else [reads_1]
    if not reads:
        raise ValueError("bbnorm needs at least one FASTQ")
    if reads_2 is not None and len(reads) != 1:
        raise ValueError("paired-end bbnorm takes one reads_1 FASTQ")
    results = await bbnorm_cmd(reads=reads, reads_2=reads_2, prefix=prefix, args=args)

    async def pick(name: str) -> File:
        found = await results.get_file(name)
        if found is None:
            raise FileNotFoundError(f"bbnorm wrote no {name}")
        return found

    names = [f"{prefix}.fastq.gz"] if reads_2 is None else [f"{prefix}_1.nm.fastq.gz", f"{prefix}_2.nm.fastq.gz"]
    return BBNormResult(reads=[await pick(name) for name in names], log=await pick(f"{prefix}.bbnorm.log"))


@dataclass
class FilterByNameResult:
    reads: list[File]  # <prefix>.<format>, or <prefix>_1.<format> and <prefix>_2.<format>
    log: File  # <prefix>.log


async def bbmap_filterbyname(
    reads_1: File,
    reads_2: File | None = None,
    names: str = "",
    output_format: str = "fastq.gz",
    single_end: bool | None = None,
    interleaved_output: bool = False,
    prefix: str = "filterbyname",
    args: str = "",
) -> FilterByNameResult:
    """Filter one sample's reads by name with FilterByName.

    ``names`` is a comma-separated list of read names (or a name file's path), dropped unless
    ``args`` has ``include=t``. ``single_end`` defaults to ``reads_2`` being None; with
    ``single_end=False`` and no ``reads_2``, ``reads_1`` holds interleaved pairs. Paired reads
    are written as two mates unless ``interleaved_output``.
    """
    if single_end is None:
        single_end = reads_2 is None
    if single_end and reads_2 is not None:
        raise ValueError("single-end filterbyname takes no reads_2")
    if not output_format or "/" in output_format or output_format.startswith("."):
        raise ValueError(f"output_format must be an extension such as fastq.gz, not {output_format!r}")
    results = await bbmap_filterbyname_cmd(
        reads_1=reads_1,
        reads_2=reads_2,
        names=names,
        output_format=output_format,
        single_end=single_end,
        interleaved_output=interleaved_output,
        prefix=prefix,
        args=args,
    )

    async def pick(name: str) -> File:
        found = await results.get_file(name)
        if found is None:
            raise FileNotFoundError(f"filterbyname wrote no {name}")
        return found

    if single_end or interleaved_output:
        names_out = [f"{prefix}.{output_format}"]
    else:
        names_out = [f"{prefix}_1.{output_format}", f"{prefix}_2.{output_format}"]
    return FilterByNameResult(reads=[await pick(name) for name in names_out], log=await pick(f"{prefix}.log"))
