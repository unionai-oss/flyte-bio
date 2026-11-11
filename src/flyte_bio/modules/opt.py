"""opt — off-target probe tracker: aligns probe sequences to transcripts to find off-target probe activity.

opt (JEFworks Lab, GPL-3.0) checks the probe panel of a targeted spatial
transcriptomics experiment (e.g. Xenium) against a reference annotation.
Exposes:

- :func:`opt_flip` — reverse-complement probes that align to the opposite
  strand of their intended target genes (``opt flip``):
  ``<prefix>/fwd_oriented.fa``.
- :func:`opt_track` — align forward-oriented probes to the reference
  transcripts and record every probe-target hit (``opt track``):
  ``<prefix>/probe2targets.tsv``.
- :func:`opt_stat` — summarise those hits per probe and gene, collapsing gene
  synonyms when a synonym table is given (``opt stat``):
  ``<prefix>/collapsed_summary.tsv``.

As upstream, ``flip`` and ``track`` get the task's CPUs (``-p``), and ``args``
are appended after the subcommand's options.
"""

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

OPT_IMAGE = "community.wave.seqera.io/library/opt:0.0.1--f0b1e63f50e38ab1"
# Upstream's process_high; the test panel runs in a few CPUs and GiB.
OPT_RESOURCES = flyte.Resources(cpu=4, memory="8Gi")

opt_flip_cmd = shell.create(
    name="opt_flip",
    image=OPT_IMAGE,
    resources=OPT_RESOURCES,
    inputs={"probes_fasta": File, "annotation": File, "fasta": File, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        PROBES=({{inputs.probes_fasta}})
        ANNOT=({{inputs.annotation}})
        FA=({{inputs.fasta}})
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        opt \
            -o "$P" \
            -p $TASK_CPUS \
            flip \
            -i "${{PROBES[0]}}" \
            -a "${{ANNOT[0]}}" \
            -f "${{FA[0]}}" \
            $ARGS || exit 1
        mv "$P/fwd_oriented.fa" "$OUT"/
    """,
)

opt_track_cmd = shell.create(
    name="opt_track",
    image=OPT_IMAGE,
    resources=OPT_RESOURCES,
    inputs={"fwd_oriented_fasta": File, "annotation": File, "fasta": File, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        PROBES=({{inputs.fwd_oriented_fasta}})
        ANNOT=({{inputs.annotation}})
        FA=({{inputs.fasta}})
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        opt \
            -o "$P" \
            -p $TASK_CPUS \
            track \
            -q "${{PROBES[0]}}" \
            -a "${{ANNOT[0]}}" \
            -t "${{FA[0]}}" \
            $ARGS || exit 1
        mv "$P/probe2targets.tsv" "$OUT"/
    """,
)

opt_stat_cmd = shell.create(
    name="opt_stat",
    image=OPT_IMAGE,
    resources=OPT_RESOURCES,
    inputs={
        "probe_targets": File,
        "fwd_oriented_fasta": File,
        "gene_synonyms": File | None,
        "prefix": str,
        "args": str,
    },
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        TARGETS=({{inputs.probe_targets}})
        PROBES=({{inputs.fwd_oriented_fasta}})
        shopt -s nullglob; SYNONYMS=({{inputs.gene_synonyms}}); shopt -u nullglob
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        EXTRA=()
        if [ ${{#SYNONYMS[@]}} -gt 0 ]; then EXTRA+=(-s "${{SYNONYMS[0]}}"); fi
        opt \
            -o "$P" \
            stat \
            -i "${{TARGETS[0]}}" \
            -q "${{PROBES[0]}}" \
            "${{EXTRA[@]}}" \
            $ARGS || exit 1
        mv "$P/collapsed_summary.tsv" "$OUT"/
    """,
)

env = flyte.TaskEnvironment.from_task("opt", opt_flip_cmd.as_task(), opt_track_cmd.as_task(), opt_stat_cmd.as_task())


async def opt_output(results: Dir, name: str, tool: str) -> File:
    found = await results.get_file(name)
    if found is None:
        raise FileNotFoundError(f"{tool} wrote no {name}")
    return found


async def opt_flip(probes_fasta: File, annotation: File, fasta: File, prefix: str = "opt", args: str = "") -> File:
    """``probes_fasta`` with probes on the wrong strand of their targets reverse-complemented: ``fwd_oriented.fa``.

    ``annotation`` is the reference annotation (GFF or GTF) and ``fasta`` its genome.
    """
    results = await opt_flip_cmd(
        probes_fasta=probes_fasta, annotation=annotation, fasta=fasta, prefix=prefix, args=args
    )
    return await opt_output(results, "fwd_oriented.fa", "opt flip")


async def opt_track(
    fwd_oriented_fasta: File, annotation: File, fasta: File, prefix: str = "opt", args: str = ""
) -> File:
    """Every probe-to-target alignment of :func:`opt_flip`'s probes: ``probe2targets.tsv``."""
    results = await opt_track_cmd(
        fwd_oriented_fasta=fwd_oriented_fasta, annotation=annotation, fasta=fasta, prefix=prefix, args=args
    )
    return await opt_output(results, "probe2targets.tsv", "opt track")


async def opt_stat(
    probe_targets: File,
    fwd_oriented_fasta: File,
    gene_synonyms: File | None = None,
    prefix: str = "opt",
    args: str = "",
) -> File:
    """A per-probe summary of :func:`opt_track`'s hits: ``collapsed_summary.tsv``."""
    results = await opt_stat_cmd(
        probe_targets=probe_targets,
        fwd_oriented_fasta=fwd_oriented_fasta,
        gene_synonyms=gene_synonyms,
        prefix=prefix,
        args=args,
    )
    return await opt_output(results, "collapsed_summary.tsv", "opt stat")
