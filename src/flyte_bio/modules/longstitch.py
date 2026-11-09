"""longstitch — LongStitch, correcting and scaffolding a draft assembly with long reads.

Exposes :func:`longstitch`, which runs one of LongStitch's pipelines on a
draft assembly (FASTA, plain or gzipped; refolded to 120 columns) with long
reads (FASTA/FASTQ, plain or gzipped):

- ``run`` or ``tigmint-ntLink``: misassembly correction (Tigmint-long), then
  scaffolding (ntLink);
- ``tigmint-ntLink-arks``: the same, then further scaffolding (ARKS-long
  with LINKS);
- ``ntLink-arks``: scaffolding only.

``span`` or ``genomesize`` (``G``) must be given; ``longmap`` is ``ont``
(default), ``pb`` or ``hifi``; ``t`` follows the task's CPUs. The thread
count changes Tigmint's corrections, and so everything after them.

As upstream, intermediate outputs keep LongStitch's names: those of the
correction step start with ``assembly.cut250.tigmint.fa.``, and the others'
``assembly.fa.`` becomes ``<prefix>.``. They are the ``*.n*.scaffold.dot``
graph, ntLink's ``*.stitch.path``, ``*.stitch.abyss-scaffold.fa``,
``*.trimmed_scafs.{agp,fa,path,tsv}``, ``*.verbose_mapping.tsv``,
``*.ntLink.scaffolds.fa``, the minimizer ``*.k*.w???.tsv``, and with ARKS
its ``*_main.tsv``, ``*_original.gv``, ``*.dist.gv``,
``*.tigpair_checkpoint.tsv`` and ``*.ntLink.scaffolds.renamed.fa``, and
LINKS' ``*_l*.{scaffolds,scaffolds.fa,assembly_correspondence.tsv,gv,log}``.

The final scaffolds are renamed to ``<prefix>.tigmint-ntLink.fa``,
``<prefix>.tigmint-ntLink-arks.fa`` or ``<prefix>.ntLink-arks.fa``, with
their headers cut to the scaffold name. Upstream intends this, but its
shell tests for ``run``, ``tigmint-ntLink`` and ``ntLink-arks`` are
malformed, so there only the ``tigmint-ntLink-arks`` scaffolds get renamed;
here each pipeline's final scaffolds are.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

LONGSTITCH_IMAGE = "quay.io/biocontainers/longstitch:1.0.5--hdfd78af_0"
DEFAULT_RESOURCES = flyte.Resources(cpu=4, memory="16Gi")
LONGSTITCH_COMMANDS = ("run", "tigmint-ntLink-arks", "tigmint-ntLink", "ntLink-arks")
LONGSTITCH_LONGMAPS = ("ont", "pb", "hifi")

# The output patterns upstream collects.
OUTPUT_GLOBS = " ".join(
    f"./*k*.w*.z*.{suffix}"
    for suffix in [
        "n*.scaffold.dot",
        "ntLink.scaffolds_*.dist.gv",
        "ntLink.scaffolds_*_l*.assembly_correspondence.tsv",
        "ntLink.scaffolds_*_l*.gv",
        "ntLink.scaffolds_*_l*.log",
        "ntLink.scaffolds_*_l*.scaffolds",
        "ntLink.scaffolds_*_l*.scaffolds.fa",
        "ntLink.scaffolds_*_main.tsv",
        "ntLink.scaffolds_*_original.gv",
        "ntLink.scaffolds_*.tigpair_checkpoint.tsv",
        "ntLink.scaffolds.fa",
        "ntLink.scaffolds.renamed.fa",
        "stitch.abyss-scaffold.fa",
        "stitch.path",
        "trimmed_scafs.agp",
        "trimmed_scafs.fa",
        "trimmed_scafs.path",
        "trimmed_scafs.tsv",
        "verbose_mapping.tsv",
    ]
)
OUTPUT_GLOBS += " ./*.k*.w???.tsv"

longstitch_cmd = shell.create(
    name="longstitch",
    image=LONGSTITCH_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "assembly": File,
        "reads": File,
        "command": str,
        "span": str,
        "genomesize": str,
        "longmap": str,
        "prefix": str,
        "args": str,
    },
    defaults={"span": "", "genomesize": "", "longmap": "ont", "args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        ASM=({{inputs.assembly}}); READS=({{inputs.reads}})
        CMD={{inputs.command}}
        SPAN={{inputs.span}}
        GENOMESIZE={{inputs.genomesize}}
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"

        if [[ "${{ASM[0]}}" == *.gz ]]; then
            zcat "${{ASM[0]}}" | fold -w 120 > assembly.fa
        else
            fold -w 120 "${{ASM[0]}}" > assembly.fa
        fi
        if [[ "${{READS[0]}}" == *.gz ]]; then
            ln -s "${{READS[0]}}" reads.fq.gz
        else
            gzip -c "${{READS[0]}}" > reads.fq.gz
        fi
        OPTIONAL=()
        if [ -n "$SPAN" ]; then OPTIONAL+=("span=$SPAN"); fi
        if [ -n "$GENOMESIZE" ]; then OPTIONAL+=("G=$GENOMESIZE"); fi

        longstitch "$CMD" draft=assembly reads=reads longmap={{inputs.longmap}} "${{OPTIONAL[@]}}" \
            t=$TASK_CPUS out_prefix="$P" $ARGS > longstitch.log 2>&1 || {{ tail -n 30 longstitch.log >&2; exit 1; }}
        cat longstitch.log >&2

        shopt -s nullglob
        for PIPELINE in tigmint-ntLink tigmint-ntLink-arks ntLink-arks; do
            for F in ./*."$PIPELINE".longstitch-scaffolds.fa; do
                sed 's/\(scaffold[0-9]*\),.*/\1/' "$F" > "$OUT/$P.$PIPELINE.fa"
            done
        done
        for F in assembly.fa.*; do mv "$F" "${{F/assembly.fa./$P.}}"; done
        # As upstream, a link left dangling by the renaming isn't an output.
        for F in {OUTPUT_GLOBS}; do
            if [ -e "$F" ] && [ ! -e "$OUT/$(basename "$F")" ]; then cp -L "$F" "$OUT"/; fi
        done
    """,
)

env = flyte.TaskEnvironment.from_task("longstitch", longstitch_cmd.as_task())


@dataclass
class LongstitchResult:
    scaffolds: File  # <prefix>.<pipeline>.fa, the final scaffolds
    outputs: dict[str, File]  # every output, by file name
    results: Dir  # everything above


async def longstitch(
    assembly: File,
    reads: File,
    command: str = "run",
    span: str = "",
    genomesize: str = "",
    longmap: str = "ont",
    prefix: str = "sample",
    args: str = "",
) -> LongstitchResult:
    """Correct and/or scaffold ``assembly`` with ``reads`` using LongStitch's ``command`` pipeline."""
    if longmap not in LONGSTITCH_LONGMAPS:
        raise ValueError(f"Unrecognised longmap option {longmap!r}. Options: {', '.join(LONGSTITCH_LONGMAPS)}")
    if not span and not genomesize:
        raise ValueError("longstitch requires either span or genomesize")
    if command not in LONGSTITCH_COMMANDS:
        raise ValueError(f"Unrecognised command {command!r}. Options: {', '.join(LONGSTITCH_COMMANDS)}")
    results = await longstitch_cmd(
        assembly=assembly,
        reads=reads,
        command=command,
        span=span,
        genomesize=genomesize,
        longmap=longmap,
        prefix=prefix,
        args=args,
    )
    return await read_longstitch_results(results, command, prefix)


async def read_longstitch_results(results: Dir, command: str, prefix: str) -> LongstitchResult:
    """Pick the outputs out of ``longstitch_cmd``'s results Dir (e.g. from an overridden task)."""
    outputs = {f.name: f async for f in results.walk()}
    pipeline = "tigmint-ntLink" if command == "run" else command
    final = f"{prefix}.{pipeline}.fa"
    if final not in outputs:
        raise FileNotFoundError(f"longstitch wrote no {final}")
    return LongstitchResult(scaffolds=outputs[final], outputs=outputs, results=results)
