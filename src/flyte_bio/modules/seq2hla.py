"""seq2hla — HLA typing and expression from paired-end RNA-seq reads.

Exposes :func:`seq2hla`, which maps one sample's reads against seq2HLA's HLA
allele references with bowtie and reports, for the classical class I genes
(``ClassI-class``: HLA-A, -B, -C), the non-classical class I genes
(``ClassI-nonclass``) and the class II genes (``ClassII``):

- ``<prefix>-<class>.HLAgenotype2digits`` / ``.HLAgenotype4digits``: the
  most likely alleles at 2- and 4-digit resolution, with confidence;
- ``<prefix>-<class>.expression``: each gene's expression (RPKM);
- ``<prefix>-<class>.bowtielog``: bowtie's mapping log;
- ``<prefix>.ambiguity``: written only when a typing is ambiguous.

As upstream, the reads are first linked to ``<prefix>_1.fastq.gz`` and
``<prefix>_2.fastq.gz``, and bowtie uses the task's CPUs.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

SEQ2HLA_IMAGE = "community.wave.seqera.io/library/seq2hla:2.3--26344ebe67dd0e8e"
DEFAULT_RESOURCES = flyte.Resources(cpu=4, memory="8Gi")

HLA_CLASSES = ("ClassI-class", "ClassI-nonclass", "ClassII")
CLASS_SUFFIXES = ("bowtielog", "expression", "HLAgenotype2digits", "HLAgenotype4digits")

# seq2HLA is a Python script; a crash can exit 255, which would be retried,
# so any failure exits 1. Intermediate files stay in the working directory.
seq2hla_cmd = shell.create(
    name="seq2hla",
    image=SEQ2HLA_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"reads_1": File, "reads_2": File, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        R1=({{inputs.reads_1}}); R2=({{inputs.reads_2}})
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        W=$(mktemp -d); cd "$W"
        ln -s "${{R1[0]}}" "${{P}}_1.fastq.gz"
        ln -s "${{R2[0]}}" "${{P}}_2.fastq.gz"
        seq2HLA -1 "${{P}}_1.fastq.gz" -2 "${{P}}_2.fastq.gz" -r "$P" -p $TASK_CPUS $ARGS || exit 1
        shopt -s nullglob
        mv ./*.bowtielog ./*.expression ./*.HLAgenotype2digits ./*.HLAgenotype4digits ./*.ambiguity \
            {{outputs.results}}/
    """,
)

env = flyte.TaskEnvironment.from_task("seq2hla", seq2hla_cmd.as_task())


@dataclass
class Seq2HlaClassResult:
    bowtielog: File
    expression: File
    genotype_2d: File  # HLAgenotype2digits
    genotype_4d: File  # HLAgenotype4digits


@dataclass
class Seq2HlaResult:
    class1: Seq2HlaClassResult  # ClassI-class: HLA-A, -B, -C
    class1_nonclass: Seq2HlaClassResult  # ClassI-nonclass: HLA-E, -F, -G, ...
    class2: Seq2HlaClassResult  # ClassII
    ambiguity: File | None  # <prefix>.ambiguity, only when a typing is ambiguous


async def seq2hla(reads_1: File, reads_2: File, prefix: str, args: str = "") -> Seq2HlaResult:
    """HLA genotypes and expression of one sample from its paired-end reads, named ``<prefix>-<class>.<kind>``."""
    results = await seq2hla_cmd(reads_1=reads_1, reads_2=reads_2, prefix=prefix, args=args)
    per_class = []
    for hla_class in HLA_CLASSES:
        files = []
        for suffix in CLASS_SUFFIXES:
            name = f"{prefix}-{hla_class}.{suffix}"
            found = await results.get_file(name)
            if found is None:
                raise FileNotFoundError(f"seq2HLA wrote no {name}")
            files.append(found)
        per_class.append(Seq2HlaClassResult(*files))
    ambiguity = await results.get_file(f"{prefix}.ambiguity")
    return Seq2HlaResult(class1=per_class[0], class1_nonclass=per_class[1], class2=per_class[2], ambiguity=ambiguity)
