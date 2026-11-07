"""crisprcleanr — CRISPRcleanR, correction of copy-number bias in CRISPR screens.

Exposes :func:`crisprcleanr_normalize`, which normalises a screen's sgRNA
counts, corrects their fold changes for gene-independent (copy-number)
effects along the genome and writes the corrected counts as
``<prefix>_norm_table.tsv``, with a ``versions.yml``.

The R code is upstream's inline script, unchanged but for taking its values
(the files, ``sample`` as the experiment name, ``min_reads``,
``min_targeted_genes`` and ``prefix``) as command-line arguments. The library
annotation is a tab-separated table keyed by its ``CODE`` column.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

CRISPRCLEANR_IMAGE = "quay.io/biocontainers/r-crisprcleanr:3.0.0--r42hdfd78af_1"
DEFAULT_RESOURCES = flyte.Resources(cpu=2, memory="8Gi")

NORMALIZE_R = r"""
args <- commandArgs(trailingOnly = TRUE)
count_file <- args[1]
library_file <- args[2]
min_reads <- as.numeric(args[3])
min_targeted_genes <- as.numeric(args[4])
sample <- args[5]
prefix <- args[6]
suppressPackageStartupMessages(library(CRISPRcleanR))
library <- read.delim(library_file, header=T,sep="\t")
row.names(library) <- library[["CODE"]]
normANDfcs <- ccr.NormfoldChanges(count_file,saveToFig = FALSE,min_reads=min_reads,EXPname=sample,
    libraryAnnotation=library,display=FALSE)
gwSortedFCs <- ccr.logFCs2chromPos(normANDfcs[["logFCs"]],library)
correctedFCs <- ccr.GWclean(gwSortedFCs,display=FALSE,label=sample)
correctedCounts <- ccr.correctCounts(sample,
                        normANDfcs[["norm_counts"]],
                        correctedFCs,
                        library,
                        minTargetedGenes=min_targeted_genes,
                        OutDir='./')

write.table(correctedCounts, file=paste0(prefix,"_norm_table.tsv"),row.names=FALSE,quote=FALSE,sep="\t")

version_crisprcleanr <- as.character(packageVersion("CRISPRcleanR"))
version_rbase <- paste(R.version[['major']],R.version[['minor']], sep = ".")
writeLines(c(
        '"CRISPRCLEANR_NORMALIZE":',
        paste('    r-base:', version_rbase),
        paste('    crisprcleanr:', version_crisprcleanr)
), 'versions.yml')
""".strip()

crisprcleanr_normalize_cmd = shell.create(
    name="crisprcleanr_normalize",
    image=CRISPRCLEANR_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "count_file": File,
        "library_file": File,
        "min_reads": int,
        "min_targeted_genes": int,
        "sample": str,
        "prefix": str,
    },
    outputs={"results": Dir},
    script=rf"""
        COUNTS=({{inputs.count_file}})
        LIB=({{inputs.library_file}})
        SAMPLE={{inputs.sample}}
        P={{inputs.prefix}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        cat > normalize.R <<'END_OF_SCRIPT'
{NORMALIZE_R}
END_OF_SCRIPT
        Rscript normalize.R "${{COUNTS[0]}}" "${{LIB[0]}}" {{inputs.min_reads}} {{inputs.min_targeted_genes}} \
            "$SAMPLE" "$P" || exit 1
        mv "${{P}}_norm_table.tsv" versions.yml "$OUT"/
    """,
)

env = flyte.TaskEnvironment.from_task("crisprcleanr", crisprcleanr_normalize_cmd.as_task())


@dataclass
class CrisprcleanrResult:
    norm_count_file: File  # <prefix>_norm_table.tsv
    versions: File  # versions.yml
    results: Dir


async def crisprcleanr_normalize(
    count_file: File,
    library_file: File,
    min_reads: int,
    min_targeted_genes: int,
    sample: str = "sample",
    prefix: str = "",
) -> CrisprcleanrResult:
    """Normalise and copy-number-correct ``count_file``'s counts; ``prefix`` defaults to ``sample``."""
    prefix = prefix or sample
    results = await crisprcleanr_normalize_cmd(
        count_file=count_file,
        library_file=library_file,
        min_reads=min_reads,
        min_targeted_genes=min_targeted_genes,
        sample=sample,
        prefix=prefix,
    )

    async def pick(name: str) -> File:
        found = await results.get_file(name)
        if found is None:
            raise FileNotFoundError(f"crisprcleanr wrote no {name}")
        return found

    return CrisprcleanrResult(
        norm_count_file=await pick(f"{prefix}_norm_table.tsv"), versions=await pick("versions.yml"), results=results
    )
