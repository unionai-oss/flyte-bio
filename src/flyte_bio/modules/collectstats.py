"""collectstats — per-sample read statistics collected into one table.

Exposes :func:`collect_stats`, which joins, per sample, the read counts of a
metagenome/metatranscriptome run into ``<prefix>.overall_stats.tsv.gz``
(sorted by sample), with a ``versions.yml``:

- ``n_post_trimming``: reads written by Trim Galore, from
  ``<sample>.*_trimming_report.txt`` (single-end) or
  ``<sample>_1.*_trimming_report.txt`` (paired-end, doubled); ``NA`` without
  trimming reports;
- ``n_non_contaminated``: the ``Result:`` reads of ``<sample>.bbduk.log``
  (the column is left out without BBDuk logs);
- ``idxs_n_mapped``, ``idxs_n_unmapped``: totals of ``<sample>.idxstats``,
  without the ``*`` line;
- one column per feature-count table ``<name>.<feature>.<rest>`` (columns
  ``sample`` and ``count``): the sample's summed counts of that feature;
- the columns of ``mergetab`` (a table keyed by ``sample``), if given.

The read counts are all expected in reads, both mates of a pair counted.
``samples`` lists every sample's id with whether it is single-end.

The R code is upstream's inline script, unchanged but for taking the samples
(a table of ``id`` and ``single_end``, written by the wrapper), the
feature-count file names and whether trimming reports and ``mergetab`` were
given as command-line arguments. As upstream, all files are staged side by
side under their own names.
"""

import tempfile
from dataclasses import dataclass
from pathlib import Path

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

COLLECTSTATS_IMAGE = "community.wave.seqera.io/library/r-base_r-dplyr_r-purrr_r-readr_pruned:24026ac6303954f2"
# Upstream: process_medium; the script is single-threaded.
DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory="8Gi")

COLLECTSTATS_R = r"""
args <- commandArgs(trailingOnly = TRUE)
samples_file <- args[1]
outfile <- args[2]
have_trimlogs <- args[3] == 'true'
mergetab_file <- args[4]
fcs <- readLines(args[5])

# All read counts joined into this table (n_post_trimming, n_non_contaminated,
# idxs_n_mapped/idxs_n_unmapped, and the featureCounts columns) are expected to use
# the same unit -- individual reads, with both mates of a pair counted separately,
# never fragments/pairs as a single unit -- so the columns stay directly comparable.
# This depends on how the upstream trimming/decontamination/alignment/counting steps
# are invoked (e.g. whether a paired-end-aware counting tool is asked to count reads
# or fragments), not on anything this module itself controls.

library(dplyr)
library(readr)
library(purrr)
library(tidyr)
library(stringr)

samples_meta <- read_tsv(samples_file, col_types = 'cl')

start    <- tibble(sample = samples_meta$id)

if ( have_trimlogs ) {
    trimming <- bind_rows(
        samples_meta %>% filter(single_end) %>%
            transmute(sample = id, tlogglob = paste0(id, '.*_trimming_report.txt'), mult = 1),
        samples_meta %>% filter(! single_end) %>%
            transmute(sample = id, tlogglob = paste0(id, '_1.*_trimming_report.txt'), mult = 2)
    ) %>%
        mutate(
            d = map(
                tlogglob,
                function(g) {
                    # Deliberately using Sys.glob() + readLines() rather than shelling
                    # out to grep/sed via pipe(sprintf(...)): the glob and file content
                    # are only ever passed as R function arguments here, never
                    # interpolated into a shell command string.
                    f <- Sys.glob(g)[1]
                    line <- grep('Reads written \\(passing filters\\)', readLines(f), value = TRUE)
                    n <- as.integer(gsub(',', '', sub(' .*', '', sub('.*: *', '', line))))
                    tibble(n_post_trimming = n)
                }
            )
        ) %>%
        unnest(d) %>%
        transmute(sample, n_post_trimming = n_post_trimming * mult)
} else {
    trimming <- tibble(sample = character(), n_post_trimming = integer())
}

# Reading each idxstats file directly instead of shelling out to grep, to avoid
# passing a sample-derived filename through a shell command string.
idxs <- tibble(fname = Sys.glob('*.idxstats')) %>%
    mutate(
        sample = str_remove(basename(fname), '.idxstats'),
        d = map(
            fname,
            \(f) {
                read_tsv(f, col_names = c('chr', 'length', 'idxs_n_mapped', 'idxs_n_unmapped'), col_types = 'ciii') %>%
                    filter(chr != '*') %>%
                    summarise(idxs_n_mapped = sum(idxs_n_mapped), idxs_n_unmapped = sum(idxs_n_unmapped))
            }
        )
    ) %>%
    unnest(d) %>%
    select(-fname)

# Column types given by name rather than position: this only depends on the
# feature-counts input's columns existing and being correctly typed, not also on
# their order, so a future column reorder upstream won't silently misparse here.
# Only sample and count are read, one file at a time, so memory does not grow with
# the number of files, and the files need not have the same other columns.
counts <- map(
    fcs,
    function(f) {
        read_tsv(
            f,
            col_types = cols_only(sample = col_character(), count = col_integer())
        ) %>%
            group_by(sample) %>%
            summarise(count = sum(count)) %>%
            mutate(feature = str_replace(f, '^[^.]+\\.([^.]+)\\..*', '\\1'))
    }
) %>%
    list_rbind() %>%
    group_by(sample, feature) %>%
    summarise(count = sum(count), .groups = 'drop') %>%
    pivot_wider(names_from = feature, values_from = count, values_fill = 0)

# Reading each bbduk log directly (instead of shelling out to grep/sed) avoids
# passing a sample-derived filename through a shell command string.
bbduk <- tibble(sample = character(), n_non_contaminated = integer())
for ( f in Sys.glob('*.bbduk.log') ) {
    s <- str_remove(f, '.bbduk.log')
    line <- grep('Result:', readLines(f), value = TRUE)
    n <- as.integer(sub(' reads.*', '', sub('Result:[[:space:]]*', '', line)))
    bbduk <- bbduk %>% union(tibble(n_non_contaminated = n, sample = s))
}
if ( nrow(bbduk) == 0 ) bbduk <- bbduk %>% select(sample)

# Add in stats from taxonomy and function, if provided
if ( mergetab_file != '' ) {
    mergetab <- read_tsv(mergetab_file, show_col_types = FALSE)
} else {
    mergetab <- tibble(sample = character())
}

# Write output
start %>%
    left_join(trimming, by = join_by(sample)) %>%
    left_join(bbduk, by = join_by(sample)) %>%
    left_join(idxs, by = join_by(sample)) %>%
    left_join(counts, by = join_by(sample)) %>%
    left_join(mergetab, by = join_by(sample)) %>%
    arrange(sample) %>%
    write_tsv(outfile)

writeLines(
    c(
        "\"CUSTOM_COLLECTSTATS\":",
        paste0("    R: ", paste0(R.Version()[c("major","minor")], collapse = ".")),
        paste0("    dplyr: ", packageVersion('dplyr')),
        paste0("    readr: ", packageVersion('readr')),
        paste0("    purrr: ", packageVersion('purrr')),
        paste0("    tidyr: ", packageVersion('tidyr')),
        paste0("    stringr: ", packageVersion('stringr'))
    ),
    "versions.yml"
)
""".strip()

collect_stats_cmd = shell.create(
    name="collect_stats",
    image=COLLECTSTATS_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "samples": File,
        "trimlogs": list[File],
        "bblogs": list[File],
        "idxstats": list[File],
        "featurecounts": list[File],
        "mergetab": File | None,
        "prefix": str,
    },
    outputs={"results": Dir},
    script=rf"""
        SAMPLES=({{inputs.samples}})
        shopt -s nullglob
        TRIMLOGS=({{inputs.trimlogs}}); BBLOGS=({{inputs.bblogs}}); IDXSTATS=({{inputs.idxstats}})
        FCS=({{inputs.featurecounts}}); MERGETAB=({{inputs.mergetab}})
        shopt -u nullglob
        P={{inputs.prefix}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        for F in "${{TRIMLOGS[@]}}" "${{BBLOGS[@]}}" "${{IDXSTATS[@]}}" "${{FCS[@]}}" "${{MERGETAB[@]}}"; do
            ln -s "$F" .
        done
        : > featurecounts.list
        for F in "${{FCS[@]}}"; do basename "$F" >> featurecounts.list; done
        HAVE_TRIMLOGS=false
        if [ ${{#TRIMLOGS[@]}} -gt 0 ]; then HAVE_TRIMLOGS=true; fi
        MERGETAB_NAME=""
        if [ ${{#MERGETAB[@]}} -gt 0 ]; then MERGETAB_NAME=$(basename "${{MERGETAB[0]}}"); fi
        cat > collectstats.R <<'END_OF_SCRIPT'
{COLLECTSTATS_R}
END_OF_SCRIPT
        Rscript collectstats.R "${{SAMPLES[0]}}" "$P.overall_stats.tsv.gz" $HAVE_TRIMLOGS "$MERGETAB_NAME" \
            featurecounts.list || exit 1
        mv "$P.overall_stats.tsv.gz" versions.yml "$OUT"/
    """,
)

env = flyte.TaskEnvironment.from_task("collectstats", collect_stats_cmd.as_task())


@dataclass
class CollectedStats:
    overall_stats: File  # <prefix>.overall_stats.tsv.gz
    versions: File  # versions.yml


async def collect_stats(
    samples: list[tuple[str, bool]],
    idxstats: list[File],
    featurecounts: list[File],
    trimlogs: list[File] | None = None,
    bblogs: list[File] | None = None,
    mergetab: File | None = None,
    prefix: str = "run",
) -> CollectedStats:
    """Join ``samples``' (id, single-end) read statistics into ``<prefix>.overall_stats.tsv.gz``."""
    table = Path(tempfile.mkdtemp()) / "samples.tsv"
    rows = [f"{sample_id}\t{'TRUE' if single_end else 'FALSE'}" for sample_id, single_end in samples]
    table.write_text("\n".join(["id\tsingle_end", *rows]) + "\n")
    results = await collect_stats_cmd(
        samples=await File.from_local(str(table)),
        trimlogs=trimlogs or [],
        bblogs=bblogs or [],
        idxstats=idxstats,
        featurecounts=featurecounts,
        mergetab=mergetab,
        prefix=prefix,
    )

    async def pick(name: str) -> File:
        found = await results.get_file(name)
        if found is None:
            raise FileNotFoundError(f"collect_stats wrote no {name}")
        return found

    return CollectedStats(
        overall_stats=await pick(f"{prefix}.overall_stats.tsv.gz"), versions=await pick("versions.yml")
    )
