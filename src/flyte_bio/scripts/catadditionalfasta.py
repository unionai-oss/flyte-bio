#!/usr/bin/env python
"""Vendored helper: append an extra FASTA (e.g. spike-ins) to a genome.

Generates a minimal GTF describing each record of the additional FASTA,
then concatenates the additional FASTA onto the genome FASTA and the
generated GTF onto the genome GTF. Runs as a stock-python-biocontainer
shell task (staged in as a File input) — see
:mod:`flyte_bio.modules.catadditionalfasta` — so no tool pod needs
``flyte_bio`` installed.

The FASTA parsing and GTF-line generation are preserved verbatim from the
upstream origin (MIT-licensed; originally by Pranathi Vemuri, modified by
Jonathan Manning). Only the entry point changed: the Nextflow ``template``
is replaced by an argparse CLI.
"""



import argparse
import os
from itertools import groupby
from typing import Iterator, Tuple


def parse_fasta(fasta_file: str) -> Iterator[Tuple[str, str]]:
    """Parse a FASTA file, yielding ``(header, sequence)`` tuples."""
    with open(fasta_file) as file_handle:
        fasta_iter = (x[1] for x in groupby(file_handle, lambda line: line[0] == ">"))
        for header in fasta_iter:
            header_str = next(header)[1:].strip()
            sequence = "".join(s.strip() for s in next(fasta_iter))
            yield (header_str, sequence)


def generate_gtf_line(name: str, length: int, biotype: str) -> str:
    """Generate one GTF line for a sequence of the given name and length."""
    biotype_attr = f' {biotype} "transgene";' if biotype else ""
    attributes = f'exon_id "{name}.1"; exon_number "1";{biotype_attr} gene_id "{name}_gene"; gene_name "{name}_gene"; gene_source "custom"; transcript_id "{name}_gene"; transcript_name "{name}_gene";\n'
    return f"{name}\ttransgene\texon\t1\t{length}\t.\t+\t.\t{attributes}"


def fasta_to_gtf(fasta: str, output_file: str, biotype: str) -> None:
    """Write a GTF describing every record in ``fasta``."""
    lines = []
    for header, sequence in parse_fasta(fasta):
        seq_name = header.split()[0].replace(" ", "_")
        lines.append(generate_gtf_line(seq_name, len(sequence), biotype))
    with open(output_file, "w") as file_handle:
        file_handle.writelines(lines)


def concat_bytes(parts: list[str], out_path: str) -> None:
    """Byte-for-byte concatenate ``parts`` into ``out_path`` (like ``cat``)."""
    with open(out_path, "wb") as out:
        for part in parts:
            with open(part, "rb") as src:
                out.write(src.read())


def cat_additional_fasta(fasta: str, gtf: str, add_fasta: str, biotype: str, out_fasta: str, out_gtf: str) -> None:
    """Append ``add_fasta`` (and a generated GTF for it) onto a genome pair."""
    add_name = os.path.splitext(os.path.basename(add_fasta))[0]
    add_gtf = f"{add_name}.gtf"
    fasta_to_gtf(add_fasta, add_gtf, biotype)
    concat_bytes([fasta, add_fasta], out_fasta)
    concat_bytes([gtf, add_gtf], out_gtf)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--fasta", required=True)
    parser.add_argument("--gtf", required=True)
    parser.add_argument("--add-fasta", required=True)
    parser.add_argument("--biotype", default="")
    parser.add_argument("--out-fasta", required=True)
    parser.add_argument("--out-gtf", required=True)
    args = parser.parse_args()
    cat_additional_fasta(args.fasta, args.gtf, args.add_fasta, args.biotype, args.out_fasta, args.out_gtf)
