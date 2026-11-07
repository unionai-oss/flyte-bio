"""vsearch — clustering sequences with VSEARCH.

Exposes :func:`vsearch_cluster`, which clusters the sequences in ``fasta``
(FASTA or FASTQ, may be gzipped) with one of VSEARCH's clustering commands
(``clustering``: ``--cluster_fast``, ``--cluster_size``, ``--cluster_smallmem``
or ``--cluster_unoise``) and writes one output, chosen by ``output``. Its name
follows upstream's: ``<prefix>.<ext>.gz`` (gzipped with ``gzip -n``), where
``ext`` is

- ``aln`` for ``--alnout``, ``biom`` for ``--biomout``, ``blast.tsv`` for
  ``--blast6out``, ``centroids.fasta`` for ``--centroids``, ``mothur.tsv`` for
  ``--mothur_shared_out``, ``msa.fasta`` for ``--msaout``, ``otu.tsv`` for
  ``--otutabout``, ``profile.txt`` for ``--profile``, ``uc.tsv`` for ``--uc``
  and ``out.tsv`` for ``--userout``;
- for ``--clusters``, one FASTA per cluster, ``<prefix>.clusters.fasta<n>.gz``;
- for ``--samout``, a BAM, ``<prefix>.bam``, converted with samtools against
  the input sequences.

Options such as ``--id 0.8`` go in ``args``.

It also exposes :func:`vsearch_sort`, which sorts the sequences in ``fasta``
(FASTA or FASTQ, may be gzipped) by abundance (``--sortbysize``) or length
(``--sortbylength``) and writes them as FASTA to ``<prefix>.fasta``; a prefix
that would overwrite the input is an error. It runs in upstream's
VSEARCH-only image.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

VSEARCH_IMAGE = "community.wave.seqera.io/library/samtools_vsearch:36d56b7ebe9ea248"
# Upstream's process_low.
DEFAULT_RESOURCES = flyte.Resources(cpu=2, memory="4Gi")
CLUSTERINGS = ("--cluster_fast", "--cluster_size", "--cluster_smallmem", "--cluster_unoise")
OUTPUT_EXTENSIONS = {
    "--alnout": "aln",
    "--biomout": "biom",
    "--blast6out": "blast.tsv",
    "--centroids": "centroids.fasta",
    "--clusters": "clusters.fasta",
    "--mothur_shared_out": "mothur.tsv",
    "--msaout": "msa.fasta",
    "--otutabout": "otu.tsv",
    "--profile": "profile.txt",
    "--samout": "sam",
    "--uc": "uc.tsv",
    "--userout": "out.tsv",
}

vsearch_cluster_cmd = shell.create(
    name="vsearch_cluster",
    image=VSEARCH_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"fasta": File, "clustering": str, "output": str, "extension": str, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        FA=({{inputs.fasta}})
        CLUSTERING={{inputs.clustering}}
        OUTPUT={{inputs.output}}
        EXT={{inputs.extension}}
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        # The SAM header records the command line, so the input keeps its name.
        NAME=$(basename "${{FA[0]}}")
        ln -s "${{FA[0]}}" "$NAME"

        vsearch \
            $CLUSTERING "$NAME" \
            $OUTPUT "$P.$EXT" \
            --threads $TASK_CPUS \
            $ARGS || exit 1

        if [ "$OUTPUT" = --clusters ]; then
            find . -type f -name "$P.$EXT*[0-9]" | xargs gzip -n
            mv "$P.$EXT"*.gz "$OUT"/
        elif [ "$OUTPUT" != --samout ]; then
            gzip -n "$P.$EXT" && mv "$P.$EXT.gz" "$OUT"/
        else
            # samtools needs an uncompressed (or bgzipped) reference to build @SQ lines
            gzip -cdf "$NAME" > "$P.reference"
            samtools view -T "$P.reference" -S -b "$P.$EXT" > "$OUT/$P.bam"
        fi
    """,
)

VSEARCH_SORT_IMAGE = "community.wave.seqera.io/library/vsearch:2.31.0--80e3b2da3ae2f3d8"
SORTS = ("--sortbysize", "--sortbylength")

vsearch_sort_cmd = shell.create(
    name="vsearch_sort",
    image=VSEARCH_SORT_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"fasta": File, "sort": str, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        FA=({{inputs.fasta}})
        SORT={{inputs.sort}}
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        vsearch \
            $SORT "${{FA[0]}}" \
            --threads $TASK_CPUS \
            --output "$OUT/$P.fasta" \
            $ARGS || exit 1
    """,
)

vsearch_cluster_env = flyte.TaskEnvironment.from_task("vsearch_samtools", vsearch_cluster_cmd.as_task())
vsearch_sort_env = flyte.TaskEnvironment.from_task("vsearch_only", vsearch_sort_cmd.as_task())
env = flyte.TaskEnvironment(name="vsearch", depends_on=[vsearch_cluster_env, vsearch_sort_env])


@dataclass
class VsearchClusterResult:
    output: File | None  # <prefix>.<ext>.gz, or <prefix>.bam for --samout; None for --clusters
    clusters: list[File]  # --clusters: <prefix>.clusters.fasta<n>.gz, in cluster order
    results: Dir  # everything above


async def vsearch_cluster(
    fasta: File, clustering: str, output: str, prefix: str = "sample", args: str = ""
) -> VsearchClusterResult:
    """Cluster ``fasta`` with ``clustering`` (``--cluster_fast``, ...), writing the ``output`` option's file."""
    if clustering not in CLUSTERINGS:
        raise ValueError(f"Unknown clustering option provided ({clustering}); expected one of {CLUSTERINGS}")
    if output not in OUTPUT_EXTENSIONS:
        raise ValueError(f"Unknown output file format provided ({output}); expected one of {list(OUTPUT_EXTENSIONS)}")
    extension = OUTPUT_EXTENSIONS[output]
    results = await vsearch_cluster_cmd(
        fasta=fasta, clustering=clustering, output=output, extension=extension, prefix=prefix, args=args
    )
    if output == "--clusters":
        stem = f"{prefix}.{extension}"
        found = [f async for f in results.walk() if f.name.startswith(stem) and f.name.endswith(".gz")]
        number = {f.name: int(f.name[len(stem) : -len(".gz")]) for f in found}
        return VsearchClusterResult(output=None, clusters=sorted(found, key=lambda f: number[f.name]), results=results)
    name = f"{prefix}.bam" if output == "--samout" else f"{prefix}.{extension}.gz"
    single = await results.get_file(name)
    if single is None:
        raise FileNotFoundError(f"vsearch wrote no {name}")
    return VsearchClusterResult(output=single, clusters=[], results=results)


async def vsearch_sort(fasta: File, sort: str, prefix: str = "sample", args: str = "") -> File:
    """``fasta``'s sequences sorted by ``sort`` (``--sortbysize`` or ``--sortbylength``): ``<prefix>.fasta``."""
    if sort not in SORTS:
        raise ValueError(f"Unknown sort option provided ({sort}); expected one of {SORTS}")
    filename = f"{prefix}.fasta"
    if fasta.path.rstrip("/").rsplit("/", 1)[-1] == filename:
        raise ValueError("Input and output names are the same, set prefix in module configuration to disambiguate!")
    results = await vsearch_sort_cmd(fasta=fasta, sort=sort, prefix=prefix, args=args)
    found = await results.get_file(filename)
    if found is None:
        raise FileNotFoundError(f"vsearch wrote no {filename}")
    return found
