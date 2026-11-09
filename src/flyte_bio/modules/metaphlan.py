"""metaphlan — MetaPhlAn 4, taxonomic profiling of metagenomes with clade-specific marker genes.

Exposes :func:`metaphlan_metaphlan`, which profiles one sample's reads (FASTQ,
single- or paired-end), contigs (FASTA), or an earlier run's bowtie2 mapping
(``*.bowtie2out.txt``) or SAM. The input type is guessed from the first file's
name, as upstream: ``.fastq``/``.fq`` is FASTQ, ``.fasta``/``.fna``/``.fa`` is
FASTA, ``.bowtie2out.txt`` is a bowtie2 mapping, anything else is SAM. Paired
FASTQ files are passed as ``R1,R2``.

The database is a tar archive (any compression) holding a MetaPhlAn bowtie2
index: it's extracted, and the index is found by its ``*.rev.1.bt2*`` file, as
upstream finds it in its database directory. Outputs keep upstream's names:

- ``<prefix>_profile.txt``: the taxonomic profile;
- ``<prefix>.biom``: the profile as a BIOM table;
- ``<prefix>.bowtie2out.txt``: the reads' marker mapping (FASTQ and FASTA
  input only), reusable as input;
- with ``save_samfile``: ``<prefix>.sam``, the bowtie2 alignments.

MetaPhlAn runs on the task's CPUs (``--nproc``).
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

METAPHLAN_IMAGE = "quay.io/biocontainers/metaphlan:4.1.1--pyhdfd78af_0"
# Upstream labels it process_medium.
DEFAULT_RESOURCES = flyte.Resources(cpu=4, memory="16Gi")

metaphlan_metaphlan_cmd = shell.create(
    name="metaphlan_metaphlan",
    image=METAPHLAN_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"input": list[File], "database": File, "save_samfile": bool, "prefix": str, "args": str},
    defaults={"save_samfile": False, "args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        INPUT=({{inputs.input}})
        DATABASE=({{inputs.database}})
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        mkdir database
        python -c 'import sys, tarfile; tarfile.open(sys.argv[1]).extractall(sys.argv[2])' \
            "${{DATABASE[0]}}" database || exit 1
        REV=$(find -L database -name '*.rev.1.bt2*' | head -n 1)
        if [ -z "$REV" ]; then echo "no bowtie2 index (*.rev.1.bt2*) in the database" >&2; exit 1; fi
        BT2_DB=$(dirname "$REV")
        BT2_DB_INDEX=$(basename "$REV" | sed 's/\.rev\.1\.bt2.*$//')

        NAME=$(basename "${{INPUT[0]}}")
        if [[ "$NAME" =~ \.(fastq|fq) ]]; then
            TYPE=fastq
        elif [[ "$NAME" =~ \.(fasta|fna|fa) ]]; then
            TYPE=fasta
        elif [[ "$NAME" == *.bowtie2out.txt ]]; then
            TYPE=bowtie2out
        else
            TYPE=sam
        fi
        DATA="${{INPUT[0]}}"
        if [ $TYPE = fastq ] && [ ${{#INPUT[@]}} -gt 1 ]; then DATA="${{INPUT[0]}},${{INPUT[1]}}"; fi
        EXTRA=()
        if [ $TYPE = fastq ] || [ $TYPE = fasta ]; then EXTRA+=(--bowtie2out "$P.bowtie2out.txt"); fi
        if [ {{inputs.save_samfile}} = true ]; then EXTRA+=(-s "$P.sam"); fi

        metaphlan \
            --nproc $TASK_CPUS \
            --input_type $TYPE \
            "$DATA" \
            $ARGS \
            "${{EXTRA[@]}}" \
            --bowtie2db "$BT2_DB" \
            --index "$BT2_DB_INDEX" \
            --biom "$P.biom" \
            --output_file "${{P}}_profile.txt" || exit 1
        shopt -s nullglob
        mv ./*_profile.txt ./*.biom ./*.bowtie2out.txt ./*.sam "$OUT"/
    """,
)

env = flyte.TaskEnvironment.from_task("metaphlan", metaphlan_metaphlan_cmd.as_task())


@dataclass
class MetaphlanResult:
    profile: File  # <prefix>_profile.txt
    biom: File  # <prefix>.biom
    bowtie2out: File | None  # <prefix>.bowtie2out.txt, for FASTQ and FASTA input
    sam: File | None  # <prefix>.sam, with save_samfile
    results: Dir  # everything MetaPhlAn wrote


async def metaphlan_metaphlan(
    input: list[File], database: File, prefix: str = "sample", save_samfile: bool = False, args: str = ""
) -> MetaphlanResult:
    """Profile one sample's ``input`` (one or two FASTQs, a FASTA, a bowtie2 mapping or a SAM) against ``database``."""
    results = await metaphlan_metaphlan_cmd(
        input=input, database=database, save_samfile=save_samfile, prefix=prefix, args=args
    )
    return await metaphlan_result(results, prefix)


async def metaphlan_result(results: Dir, prefix: str) -> MetaphlanResult:
    """Pick MetaPhlAn's outputs out of ``results``."""
    profile = await results.get_file(f"{prefix}_profile.txt")
    biom = await results.get_file(f"{prefix}.biom")
    if profile is None or biom is None:
        raise FileNotFoundError(f"MetaPhlAn wrote no {prefix}_profile.txt or {prefix}.biom")
    return MetaphlanResult(
        profile=profile,
        biom=biom,
        bowtie2out=await results.get_file(f"{prefix}.bowtie2out.txt"),
        sam=await results.get_file(f"{prefix}.sam"),
        results=results,
    )
