"""hifiasm — haplotype-resolved de novo assembly of PacBio HiFi reads.

Exposes :func:`hifiasm`, which assembles one sample's HiFi reads (one or more
files), optionally with ultra-long ONT reads (``--ul``), and phases the
assembly with parental k-mer dumps (yak, trio binning: ``-1``/``-2``) or Hi-C
reads (``--h1``/``--h2``). Earlier runs' ``.bin`` files can be passed back to
skip read correction and overlapping, as hifiasm reuses the bins it finds
under its output prefix.

As upstream, the outputs keep hifiasm's names under ``-o <prefix>``; which
ones exist depends on the mode (``bp``, ``dip`` for trio, ``hic`` for Hi-C):

- ``<prefix>.<mode>.r_utg.gfa`` (raw unitigs), ``<prefix>.<mode>.p_utg.gfa``
  (processed unitigs);
- primary contigs ``<prefix>.{p_ctg,bp.p_ctg,hic.p_ctg}.gfa`` and alternate
  contigs ``<prefix>.{a_ctg,hic.a_ctg}.gfa``;
- haplotype contigs ``<prefix>.<mode>.hap1.p_ctg.gfa`` and ``.hap2.p_ctg.gfa``;
- the ``*.bin`` files, and with ``--write-ec`` / ``--write-paf`` in ``args``
  the corrected reads ``<prefix>.ec.fa.gz`` and overlaps
  ``<prefix>.ovlp.paf.gz``;
- hifiasm's stderr, ``<prefix>.stderr.log``.

Long reads (and ultra-long reads) are passed sorted by file name, as
upstream. ``-t`` follows the task's CPUs; options in ``args`` come first.
Trio binning and Hi-C phasing can't be combined, and each needs both of its
inputs: upstream logs an error for these, here the wrapper raises.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

HIFIASM_IMAGE = "quay.io/biocontainers/hifiasm:0.25.0--h5ca1c30_0"

# Upstream's label asks for 12 CPUs and 72 GB. These defaults fit the test
# data and small genomes; raise them with hifiasm_cmd.as_task().override(...)
# for real assemblies. The 4 CPUs match upstream's test CPUs.
DEFAULT_RESOURCES = flyte.Resources(cpu=4, memory="16Gi")

hifiasm_cmd = shell.create(
    name="hifiasm",
    image=HIFIASM_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "long_reads": list[File],
        "ul_reads": list[File],
        "paternal_kmer_dump": File | None,
        "maternal_kmer_dump": File | None,
        "hic_read1": File | None,
        "hic_read2": File | None,
        "bin_files": list[File],
        "prefix": str,
        "args": str,
    },
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        LONG=({{inputs.long_reads}})
        shopt -s nullglob
        UL=({{inputs.ul_reads}}); PAT=({{inputs.paternal_kmer_dump}}); MAT=({{inputs.maternal_kmer_dump}})
        H1=({{inputs.hic_read1}}); H2=({{inputs.hic_read2}}); BINS=({{inputs.bin_files}})
        shopt -u nullglob
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"

        # hifiasm reuses the bins it finds under its output prefix.
        for B in "${{BINS[@]}}"; do ln -s "$B" "$(basename "$B")"; done
        EXTRA=()
        if [ ${{#PAT[@]}} -gt 0 ]; then EXTRA+=(-1 "${{PAT[0]}}" -2 "${{MAT[0]}}"); fi
        if [ ${{#H1[@]}} -gt 0 ]; then EXTRA+=(--h1 "${{H1[0]}}" --h2 "${{H2[0]}}"); fi
        if [ ${{#UL[@]}} -gt 0 ]; then EXTRA+=(--ul "${{UL[@]}}"); fi

        hifiasm $ARGS -t $TASK_CPUS "${{EXTRA[@]}}" -o "$P" "${{LONG[@]}}" 2> "$P.stderr.log" || {{
            tail -n 20 "$P.stderr.log" >&2; exit 1;
        }}
        cat "$P.stderr.log" >&2

        if [ -f "$P.ec.fa" ]; then gzip "$P.ec.fa"; fi
        if [ -f "$P.ovlp.paf" ]; then gzip "$P.ovlp.paf"; fi
        shopt -s nullglob
        for F in ./*.r_utg.gfa ./*.p_utg.gfa ./*.p_ctg.gfa ./*.a_ctg.gfa ./*.bin ./*.ec.fa.gz ./*.ovlp.paf.gz \
                "$P.stderr.log"; do
            cp -L "$F" "$OUT"/
        done
    """,
)

env = flyte.TaskEnvironment.from_task("hifiasm", hifiasm_cmd.as_task())


@dataclass
class HifiasmResult:
    raw_unitigs: File | None  # *.r_utg.gfa
    processed_unitigs: File | None  # *.p_utg.gfa
    primary_contigs: File | None  # <prefix>.{p_ctg,bp.p_ctg,hic.p_ctg}.gfa
    alternate_contigs: File | None  # <prefix>.{a_ctg,hic.a_ctg}.gfa
    hap1_contigs: File | None  # <prefix>.<mode>.hap1.p_ctg.gfa
    hap2_contigs: File | None  # <prefix>.<mode>.hap2.p_ctg.gfa
    corrected_reads: File | None  # <prefix>.ec.fa.gz, with --write-ec
    read_overlaps: File | None  # <prefix>.ovlp.paf.gz, with --write-paf
    bin_files: list[File]  # *.bin, to pass to a later run
    log: File  # <prefix>.stderr.log
    results: Dir  # everything above


async def hifiasm(
    long_reads: list[File],
    ul_reads: list[File] | None = None,
    paternal_kmer_dump: File | None = None,
    maternal_kmer_dump: File | None = None,
    hic_read1: File | None = None,
    hic_read2: File | None = None,
    bin_files: list[File] | None = None,
    prefix: str = "sample",
    args: str = "",
) -> HifiasmResult:
    """Assemble one sample's HiFi reads, optionally phased by parental k-mers or Hi-C reads."""
    if not long_reads:
        raise ValueError("hifiasm needs long reads")
    trio = [paternal_kmer_dump, maternal_kmer_dump]
    hic = [hic_read1, hic_read2]
    if any(f is not None for f in trio) and any(f is not None for f in hic):
        raise ValueError("hifiasm trio binning mode and Hi-C phasing can not be used at the same time")
    if any(f is not None for f in trio) and not all(f is not None for f in trio):
        raise ValueError("hifiasm trio binning needs both the paternal and the maternal k-mer dump")
    if any(f is not None for f in hic) and not all(f is not None for f in hic):
        raise ValueError("hifiasm Hi-C phasing needs both the forward and the reverse Hi-C reads")
    results = await hifiasm_cmd(
        long_reads=long_reads,
        ul_reads=ul_reads or [],
        paternal_kmer_dump=paternal_kmer_dump,
        maternal_kmer_dump=maternal_kmer_dump,
        hic_read1=hic_read1,
        hic_read2=hic_read2,
        bin_files=bin_files or [],
        prefix=prefix,
        args=args,
    )
    return await read_hifiasm_results(results, prefix)


async def read_hifiasm_results(results: Dir, prefix: str) -> HifiasmResult:
    """Pick the outputs out of ``hifiasm_cmd``'s results Dir (e.g. from an overridden task)."""
    files = {f.name: f async for f in results.walk()}

    def one(names: list[str]) -> File | None:
        found = [files[n] for n in sorted(files) if n in names]
        if len(found) > 1:
            raise FileNotFoundError(f"hifiasm wrote several of {names}")
        return found[0] if found else None

    def ending(suffix: str, start: str = "") -> File | None:
        return one([n for n in files if n.endswith(suffix) and n.startswith(start)])

    if f"{prefix}.stderr.log" not in files:
        raise FileNotFoundError(f"hifiasm wrote no {prefix}.stderr.log")
    return HifiasmResult(
        raw_unitigs=ending(".r_utg.gfa"),
        processed_unitigs=ending(".p_utg.gfa"),
        primary_contigs=one([f"{prefix}.{m}.gfa" for m in ("p_ctg", "bp.p_ctg", "hic.p_ctg")]),
        alternate_contigs=one([f"{prefix}.{m}.gfa" for m in ("a_ctg", "hic.a_ctg")]),
        hap1_contigs=ending(".hap1.p_ctg.gfa", f"{prefix}."),
        hap2_contigs=ending(".hap2.p_ctg.gfa", f"{prefix}."),
        corrected_reads=ending(".ec.fa.gz"),
        read_overlaps=ending(".ovlp.paf.gz"),
        bin_files=[files[n] for n in sorted(files) if n.endswith(".bin")],
        log=files[f"{prefix}.stderr.log"],
        results=results,
    )
