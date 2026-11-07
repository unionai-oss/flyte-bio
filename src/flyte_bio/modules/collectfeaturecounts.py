"""collectfeaturecounts — featureCounts tables collected into one long table with TPMs.

Exposes :func:`collect_featurecounts`, which reads one or more featureCounts
tables (one sample each: the table's last column) and writes
``<prefix>.counts.tsv.gz``: one row per feature with a nonzero count, columns
``orf chr start end strand length sample count tpm``. The sample name is the
count column's name without a trailing ``.sorted.bam``; TPMs are per sample,
rounded to 6 decimals. Tables are read one at a time, in file-name order, and
a sample in more than one table is an error (R's message is the task's
error message). A ``versions.yml`` is written alongside.

The R code is upstream's template, unchanged but for taking its thread count
(the task's CPUs) and ``prefix`` as command-line arguments. As upstream, the
tables are staged under ``input/`` with their own names.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

COLLECTFEATURECOUNTS_IMAGE = (
    "community.wave.seqera.io/library/r-base_r-data.table_r-dplyr_r-dtplyr_pruned:53b2279b83bc6eef"
)
# Upstream: process_low.
DEFAULT_RESOURCES = flyte.Resources(cpu=2, memory="8Gi")

COLLECTFEATURECOUNTS_R = r"""
args <- commandArgs(trailingOnly = TRUE)

library(data.table)
library(dtplyr)
library(readr)
library(dplyr)
library(stringr)

setDTthreads(as.integer(args[1]))

prefix <- args[2]

files <- Sys.glob('input/*')
if ( length(files) == 0 ) stop('No input files')

# One file at a time: memory use is set by the largest file, not by the number of samples.
# TPM is computed per file, so each sample must be in exactly one file.
out <- gzfile(paste0(prefix, '.counts.tsv.gz'), 'wb')
first <- TRUE
seen <- character()

for ( file in files ) {
    sample <- str_remove(tail(names(fread(file, sep = '\t', skip = 1, nrows = 0)), 1), '[.]sorted[.]bam$')
    if ( sample %in% seen ) stop(paste0('Sample ', sample, ' is in more than one input file'))
    seen <- c(seen, sample)
    fread(file, sep = '\t', skip = 1) %>%
        melt(measure.vars = c(ncol(.)), variable.name = 'sample', value.name = 'count') %>%
        lazy_dt() %>%
        filter(count > 0) %>%
        mutate(
            sample = str_remove(sample, '[.]sorted[.]bam$'),
            r = count/Length
        ) %>%
        rename(orf = Geneid, chr = Chr, start = Start, end = End, strand = Strand, length = Length) %>%
        group_by(sample) %>%
        # Rounded so tables independently computed elsewhere from the same underlying
        # counts (e.g. after a caller-consolidation step using a different R backend)
        # stay byte-comparable -- data.table/dtplyr and plain dplyr can otherwise round
        # the last significant digit of the same mathematical value differently.
        mutate(tpm = round(r/sum(r) * 1e6, 6)) %>% ungroup() %>%
        select(-r) %>%
        as_tibble() %>%
        write_tsv(out, col_names = first)
    first <- FALSE
}

close(out)

writeLines(
    c(
        "\"CUSTOM_COLLECTFEATURECOUNTS\":",
        paste0("    R: ", paste0(R.Version()[c("major","minor")], collapse = ".")),
        paste0("    dplyr: ", packageVersion('dplyr')),
        paste0("    readr: ", packageVersion('readr')),
        paste0("    stringr: ", packageVersion('stringr')),
        paste0("    dtplyr: ", packageVersion('dtplyr')),
        paste0("    data.table: ", packageVersion('data.table'))
    ),
    "versions.yml"
)
""".strip()

collect_featurecounts_cmd = shell.create(
    name="collect_featurecounts",
    image=COLLECTFEATURECOUNTS_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"tables": list[File], "prefix": str},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        TABLES=({{inputs.tables}})
        P={{inputs.prefix}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        mkdir input
        for F in "${{TABLES[@]}}"; do ln -s "$F" input/; done
        cat > collectfeaturecounts.R <<'END_OF_SCRIPT'
{COLLECTFEATURECOUNTS_R}
END_OF_SCRIPT
        if ! Rscript collectfeaturecounts.R "$TASK_CPUS" "$P" 2> collectfeaturecounts.err; then
            cat collectfeaturecounts.err >&2
            {{ tail -c 2048 collectfeaturecounts.err > /dev/termination-log; }} 2> /dev/null
            exit 1
        fi
        cat collectfeaturecounts.err >&2
        mv "$P.counts.tsv.gz" versions.yml "$OUT"/
    """,
)

env = flyte.TaskEnvironment.from_task("collectfeaturecounts", collect_featurecounts_cmd.as_task())


@dataclass
class CollectedFeatureCounts:
    counts: File  # <prefix>.counts.tsv.gz
    versions: File  # versions.yml


async def collect_featurecounts(tables: list[File], prefix: str = "sample") -> CollectedFeatureCounts:
    """Collect featureCounts ``tables`` (one sample each) into ``<prefix>.counts.tsv.gz`` with per-sample TPMs."""
    results = await collect_featurecounts_cmd(tables=tables, prefix=prefix)

    async def pick(name: str) -> File:
        found = await results.get_file(name)
        if found is None:
            raise FileNotFoundError(f"collect_featurecounts wrote no {name}")
        return found

    return CollectedFeatureCounts(counts=await pick(f"{prefix}.counts.tsv.gz"), versions=await pick("versions.yml"))
