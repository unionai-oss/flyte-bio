"""trnascanse — tRNAscan-SE, detection of transfer RNA genes.

Exposes :func:`trnascanse`, which scans a FASTA (plain or gzipped) for tRNA
genes. Outputs keep upstream's names:

- ``<prefix>.tsv``: the tRNAs found;
- ``<prefix>.log``: the run log (timestamped);
- ``<prefix>.stats``: summary statistics (timestamped), with ``write_stats``;
- ``<prefix>.fasta.gz``: the tRNA sequences, with ``write_fasta``;
- ``<prefix>.gff`` and ``<prefix>.bed``: the tRNAs as GFF3 and BED, with
  ``write_gff`` and ``write_bed``.

As upstream, tRNAscan-SE's configuration is copied into the working directory
with its temporary directory (``/tmp``) pointed at the working directory, so
large genomes don't fill the container's ``/tmp``. ``--thread`` follows the
task's CPUs; ``args`` (for example ``-B`` for bacteria, ``-A`` for archaea)
come before the outputs.

Upstream gzips the FASTA output but then looks for ``*.fasta``, so its
``fasta`` output is always empty; here the gzipped FASTA is returned.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

TRNASCANSE_IMAGE = "community.wave.seqera.io/library/trnascan-se_gzip_findutils:15c0475b13b49f03"
# Upstream: process_medium. The 4 CPUs match upstream's test CPUs.
DEFAULT_RESOURCES = flyte.Resources(cpu=4, memory="8Gi")

trnascanse_cmd = shell.create(
    name="trnascanse",
    image=TRNASCANSE_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "fasta": File,
        "write_stats": bool,
        "write_fasta": bool,
        "write_gff": bool,
        "write_bed": bool,
        "prefix": str,
        "args": str,
    },
    defaults={"write_stats": True, "write_fasta": False, "write_gff": False, "write_bed": False, "args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        FA=({{inputs.fasta}})
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        NAME=$(basename "${{FA[0]}}")
        if [[ "$NAME" == *.gz ]]; then
            INPUT=${{NAME%.gz}}
            gunzip -c "${{FA[0]}}" > "$INPUT"
        else
            INPUT=$NAME
            ln -s "${{FA[0]}}" "$INPUT"
        fi
        OUTPUTS=()
        if [ {{inputs.write_stats}} = true ]; then OUTPUTS+=(--stats "$P".stats); fi
        if [ {{inputs.write_fasta}} = true ]; then OUTPUTS+=(--fasta "$P".fasta); fi
        if [ {{inputs.write_gff}} = true ]; then OUTPUTS+=(--gff "$P".gff); fi
        if [ {{inputs.write_bed}} = true ]; then OUTPUTS+=(--bed "$P".bed); fi

        ## larger genomes can fill up the limited temp space in the singularity container
        ## expected location of the default config file is with the exectuable
        ## copy this and modify to use the working dir as the temp directory
        CONF=$(which tRNAscan-SE).conf
        cp "$CONF" trnascan.conf
        sed -i s#/tmp#.#g trnascan.conf

        tRNAscan-SE \
            --thread $TASK_CPUS \
            -c trnascan.conf \
            $ARGS \
            -o "$P".tsv \
            -l "$P".log \
            "${{OUTPUTS[@]}}" \
            "$INPUT" || exit 1

        if [ -f "$P".fasta ]; then gzip "$P".fasta; fi
        shopt -s nullglob
        mv "$P".tsv "$P".log "$P".stats "$P".fasta.gz "$P".gff "$P".bed "$OUT"/
    """,
)

env = flyte.TaskEnvironment.from_task("trnascanse", trnascanse_cmd.as_task())


@dataclass
class TrnascanseResult:
    tsv: File
    log: File
    stats: File | None  # with write_stats
    fasta: File | None  # gzipped, with write_fasta
    gff: File | None  # with write_gff
    bed: File | None  # with write_bed
    results: Dir  # everything above


async def trnascanse(
    fasta: File,
    write_stats: bool = True,
    write_fasta: bool = False,
    write_gff: bool = False,
    write_bed: bool = False,
    prefix: str = "sample",
    args: str = "",
) -> TrnascanseResult:
    """Scan ``fasta`` for tRNA genes."""
    results = await trnascanse_cmd(
        fasta=fasta,
        write_stats=write_stats,
        write_fasta=write_fasta,
        write_gff=write_gff,
        write_bed=write_bed,
        prefix=prefix,
        args=args,
    )

    async def pick(name: str) -> File:
        found = await results.get_file(name)
        if found is None:
            raise FileNotFoundError(f"tRNAscan-SE wrote no {name}")
        return found

    return TrnascanseResult(
        tsv=await pick(f"{prefix}.tsv"),
        log=await pick(f"{prefix}.log"),
        stats=await pick(f"{prefix}.stats") if write_stats else None,
        fasta=await pick(f"{prefix}.fasta.gz") if write_fasta else None,
        gff=await pick(f"{prefix}.gff") if write_gff else None,
        bed=await pick(f"{prefix}.bed") if write_bed else None,
        results=results,
    )
