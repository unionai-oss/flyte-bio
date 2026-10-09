#!/usr/bin/env python
"""Vendored helper: filter a GTF down to a genome's sequences.

Runs as a stock-python-biocontainer shell task (the script is staged into
the container as a File input) — see :mod:`flyte_bio.modules.gtffilter`. So
no tool pod needs ``flyte_bio`` installed.

Filtering logic preserved verbatim from its upstream origin (MIT-licensed;
originally by Olga Botvinnik, reworked by Jonathan Manning and Nico
Trummer). Only the entry point changed: the original Nextflow ``template``
with ``${...}`` interpolation is replaced by an argparse CLI.
"""



import argparse
import gzip
import logging
import os
import re
import statistics
from typing import Optional, Set

logging.basicConfig(format="%(name)s - %(asctime)s %(levelname)s: %(message)s")
logger = logging.getLogger("gtf_filter")
logger.setLevel(logging.INFO)


def extract_fasta_seq_names(fasta_name: str) -> Set[str]:
    """Extract the sequence names from a FASTA file."""
    is_gz = fasta_name.endswith(".gz")
    open_fn = gzip.open if is_gz else open

    with open_fn(fasta_name) as fasta:
        sequences = set()
        for line in fasta:
            line = line.decode("utf-8") if is_gz else line
            if line.startswith(">"):
                sequences.add(line[1:].split(None, 1)[0])

        return sequences


def tab_delimited(file: str) -> float:
    """Check a file is tab-delimited; return the median tab count per line."""
    with open(file) as f:
        data = f.read(102400)
        return statistics.median(line.count("\t") for line in data.split("\n"))


def filter_gtf(fasta: Optional[str], gtf_in: str, filtered_gtf_out: str, skip_transcript_id_check: bool) -> None:
    """Filter a GTF file based on FASTA sequence names."""
    if tab_delimited(gtf_in) != 8:
        raise ValueError("Invalid GTF file: Expected 9 tab-separated columns.")

    seq_names_in_genome = None
    if fasta and os.path.isfile(fasta):
        seq_names_in_genome = extract_fasta_seq_names(fasta)
        logger.info(f"Extracted chromosome sequence names from {fasta}")

    seq_names_in_gtf = set()
    is_gz = gtf_in.endswith(".gz")
    open_fn = gzip.open if is_gz else open
    with open_fn(gtf_in) as gtf, open_fn(filtered_gtf_out, "wb" if is_gz else "w") as out:
        line_count = 0
        for line in gtf:
            line = line.decode("utf-8") if is_gz else line
            seq_name = line.split("\t")[0]
            seq_names_in_gtf.add(seq_name)

            if seq_names_in_genome is None or seq_name in seq_names_in_genome:
                if skip_transcript_id_check or re.search(r'transcript_id "([^"]+)"', line):
                    out.write(line.encode() if is_gz else line)
                    line_count += 1

        if line_count == 0:
            raise ValueError("All GTF lines removed by filters")

    logger.info(f"Extracted {line_count} matching sequences from {gtf_in} into {filtered_gtf_out}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--gtf", required=True)
    parser.add_argument("--fasta", default="")
    parser.add_argument("--output", required=True)
    parser.add_argument("--skip-transcript-id-check", action="store_true", default=False)
    args = parser.parse_args()
    filter_gtf(args.fasta or None, args.gtf, args.output, args.skip_transcript_id_check)
