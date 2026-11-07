"""rsemmergecounts — merge per-sample RSEM results into count and TPM matrices.

Exposes :func:`rsem_merge_counts`, which takes the samples'
``<sample>.genes.results`` and ``<sample>.isoforms.results`` and writes, as
upstream (samples in file-name order, named by file name):

- ``<prefix>.gene_counts.tsv`` / ``<prefix>.gene_tpm.tsv``: gene and
  transcript ids, then one expected-count (TPM) column per sample;
- ``<prefix>.transcript_counts.tsv`` / ``<prefix>.transcript_tpm.tsv``: the
  same per transcript;
- ``<prefix>.genes_long.tsv`` / ``<prefix>.isoforms_long.tsv``: every
  sample's results stacked, with a leading ``sample_name`` column.

Shell only (cut, paste, awk), on upstream's sed image.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

RSEMMERGECOUNTS_IMAGE = "community.wave.seqera.io/library/sed:4.7--2af80cd00892ec9f"
DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory="2Gi")

# The results are linked into genes/ and isoforms/ under their own names, as
# upstream stages them; samples are then taken in `ls` order.
rsem_merge_counts_cmd = shell.create(
    name="rsem_merge_counts",
    image=RSEMMERGECOUNTS_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"genes": list[File], "isoforms": list[File], "prefix": str},
    outputs={"results": Dir},
    script=r"""
        GENES=({inputs.genes}); ISOFORMS=({inputs.isoforms})
        P={inputs.prefix}
        OUT={outputs.results}
        W=$(mktemp -d); cd "$W"
        mkdir genes isoforms
        for F in "${GENES[@]}"; do ln -s "$F" genes/"$(basename "$F")"; done
        for F in "${ISOFORMS[@]}"; do ln -s "$F" isoforms/"$(basename "$F")"; done

        mkdir -p tmp/genes
        cut -f 1,2 `ls ./genes/* | head -n 1` > gene_ids.txt
        for fileid in `ls ./genes/*`; do
            samplename=`basename $fileid | sed s/\\.genes.results$//g`
            echo $samplename > tmp/genes/${samplename}.counts.txt
            cut -f 5 ${fileid} | tail -n+2 >> tmp/genes/${samplename}.counts.txt
            echo $samplename > tmp/genes/${samplename}.tpm.txt
            cut -f 6 ${fileid} | tail -n+2 >> tmp/genes/${samplename}.tpm.txt
        done

        mkdir -p tmp/isoforms
        cut -f 1,2 `ls ./isoforms/* | head -n 1` > transcript_ids.txt
        for fileid in `ls ./isoforms/*`; do
            samplename=`basename $fileid | sed s/\\.isoforms.results$//g`
            echo $samplename > tmp/isoforms/${samplename}.counts.txt
            cut -f 5 ${fileid} | tail -n+2 >> tmp/isoforms/${samplename}.counts.txt
            echo $samplename > tmp/isoforms/${samplename}.tpm.txt
            cut -f 6 ${fileid} | tail -n+2 >> tmp/isoforms/${samplename}.tpm.txt
        done

        paste gene_ids.txt tmp/genes/*.counts.txt > $P.gene_counts.tsv
        paste gene_ids.txt tmp/genes/*.tpm.txt > $P.gene_tpm.tsv
        paste transcript_ids.txt tmp/isoforms/*.counts.txt > $P.transcript_counts.tsv
        paste transcript_ids.txt tmp/isoforms/*.tpm.txt > $P.transcript_tpm.tsv

        GENES_HEADER="sample_name\tgene_id\ttranscript_id(s)\tlength\teffective_length\texpected_count\tTPM\tFPKM"
        echo -e "$GENES_HEADER" > $P.genes_long.tsv
        for fileid in `ls ./genes/*`; do
            samplename=`basename $fileid | sed s/\\.genes.results$//g`
            tail -n+2 $fileid \
                | awk -v sample=$samplename 'BEGIN{OFS="\t"}{print sample,$1,$2,$3,$4,$5,$6,$7}' >> $P.genes_long.tsv
        done

        ISOFORMS_HEADER="sample_name\ttranscript_id\tgene_id\tlength\teffective_length\texpected_count\tTPM\tFPKM\tIsoPct"
        echo -e "$ISOFORMS_HEADER" > $P.isoforms_long.tsv
        for fileid in `ls ./isoforms/*`; do
            samplename=`basename $fileid | sed s/\\.isoforms.results$//g`
            tail -n+2 $fileid \
                | awk -v sample=$samplename 'BEGIN{OFS="\t"}{print sample,$1,$2,$3,$4,$5,$6,$7,$8}' \
                >> $P.isoforms_long.tsv
        done
        mv $P.*.tsv "$OUT"/
    """,
)

env = flyte.TaskEnvironment.from_task("rsemmergecounts", rsem_merge_counts_cmd.as_task())


@dataclass
class RsemMergedCounts:
    counts_gene: File  # <prefix>.gene_counts.tsv
    tpm_gene: File  # <prefix>.gene_tpm.tsv
    counts_transcript: File  # <prefix>.transcript_counts.tsv
    tpm_transcript: File  # <prefix>.transcript_tpm.tsv
    genes_long: File  # <prefix>.genes_long.tsv
    isoforms_long: File  # <prefix>.isoforms_long.tsv


async def rsem_merge_counts(genes: list[File], isoforms: list[File], prefix: str = "rsem") -> RsemMergedCounts:
    """Merge samples' ``*.genes.results`` and ``*.isoforms.results`` into matrices and long tables."""
    if not genes or not isoforms:
        raise ValueError("rsem_merge_counts needs gene and isoform results")
    results = await rsem_merge_counts_cmd(genes=genes, isoforms=isoforms, prefix=prefix)

    async def pick(name: str) -> File:
        found = await results.get_file(f"{prefix}.{name}")
        if found is None:
            raise FileNotFoundError(f"rsem_merge_counts wrote no {prefix}.{name}")
        return found

    return RsemMergedCounts(
        counts_gene=await pick("gene_counts.tsv"),
        tpm_gene=await pick("gene_tpm.tsv"),
        counts_transcript=await pick("transcript_counts.tsv"),
        tpm_transcript=await pick("transcript_tpm.tsv"),
        genes_long=await pick("genes_long.tsv"),
        isoforms_long=await pick("isoforms_long.tsv"),
    )
