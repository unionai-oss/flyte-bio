"""trust4 — TCR and BCR repertoire reconstruction from bulk or single-cell RNA-seq with TRUST4.

Exposes :func:`trust4`, which runs ``run-trust4`` on one sample: an aligned
BAM (``bam``), or reads, single-end (``single_end``, one FASTQ) or paired-end
(two). ``fasta`` holds the V/D/J/C gene coordinates (TRUST4's ``-f``) and
``vdj_reference`` their full IMGT sequences for annotation (``--ref``,
optional for reads). For single-cell data, ``cell_barcode_read`` and
``umi_read`` (``R1`` or ``R2``) name the read holding the cell barcode and
the UMI, ``read_format`` their layout (``--readFormat``) and
``barcode_whitelist`` the barcodes to keep.

Outputs keep upstream's names, under ``<prefix>``:

- ``<prefix>_airr.tsv``, ``<prefix>_airr_align.tsv`` and ``<prefix>_report.tsv``:
  the rearrangements (AIRR format), their alignments and the CDR3 report;
- ``<prefix>_annot.fa`` and ``<prefix>_assembled_reads.fa``: the annotated
  contigs and assembled reads;
- ``<prefix>_cdr3.out``, ``<prefix>_raw.out`` and ``<prefix>_final.out``;
- ``<prefix>_toassemble.fq`` (or ``_toassemble_1.fq`` / ``_toassemble_2.fq``):
  the candidate receptor reads;
- for single-cell data, ``<prefix>_barcode_report.tsv`` and friends.

TRUST4's threads are the task's CPUs.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

TRUST4_IMAGE = "quay.io/biocontainers/trust4:1.1.5--h5ca1c30_0"
# Upstream labels TRUST4 process_medium.
DEFAULT_RESOURCES = flyte.Resources(cpu=4, memory="8Gi")
READS = ("R1", "R2")

trust4_cmd = shell.create(
    name="trust4",
    image=TRUST4_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "bam": File | None,
        "reads_1": File | None,
        "reads_2": File | None,
        "single_end": bool,
        "fasta": File,
        "vdj_reference": File | None,
        "barcode_whitelist": File | None,
        "cell_barcode_read": str,
        "umi_read": str,
        "read_format": str,
        "prefix": str,
        "args": str,
    },
    defaults={"cell_barcode_read": "", "umi_read": "", "read_format": "", "args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        FA=({{inputs.fasta}})
        shopt -s nullglob
        BAM=({{inputs.bam}}); R1=({{inputs.reads_1}}); R2=({{inputs.reads_2}})
        REF=({{inputs.vdj_reference}}); WHITELIST=({{inputs.barcode_whitelist}})
        shopt -u nullglob
        BARCODE_READ={{inputs.cell_barcode_read}}
        UMI_READ={{inputs.umi_read}}
        READ_FORMAT={{inputs.read_format}}
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        cd "$OUT"

        EXTRA=()
        if [ ${{#BAM[@]}} -gt 0 ]; then EXTRA+=(-b "${{BAM[0]}}"); fi
        if [ ${{#R1[@]}} -gt 0 ]; then
            if [ {{inputs.single_end}} = true ]; then
                EXTRA+=(-u "${{R1[0]}}")
            else
                EXTRA+=(-1 "${{R1[0]}}" -2 "${{R2[0]}}")
            fi
        fi
        if [ "$BARCODE_READ" = R1 ]; then EXTRA+=(--barcode "${{R1[0]}}"); fi
        if [ "$BARCODE_READ" = R2 ]; then EXTRA+=(--barcode "${{R2[0]}}"); fi
        if [ -n "$READ_FORMAT" ]; then EXTRA+=(--readFormat "$READ_FORMAT"); fi
        if [ "$UMI_READ" = R1 ]; then EXTRA+=(--UMI "${{R1[0]}}"); fi
        if [ "$UMI_READ" = R2 ]; then EXTRA+=(--UMI "${{R2[0]}}"); fi
        if [ ${{#REF[@]}} -gt 0 ]; then EXTRA+=(--ref "${{REF[0]}}"); fi
        if [ ${{#WHITELIST[@]}} -gt 0 ]; then EXTRA+=(--barcodeWhitelist "${{WHITELIST[0]}}"); fi

        run-trust4 \
            "${{EXTRA[@]}}" \
            -t $TASK_CPUS \
            -f "${{FA[0]}}" \
            -o "$P" \
            $ARGS || exit 1
    """,
)

env = flyte.TaskEnvironment.from_task("trust4", trust4_cmd.as_task())


@dataclass
class Trust4Result:
    airr: File  # <prefix>_airr.tsv
    report: File  # <prefix>_report.tsv
    tsv: list[File]  # *.tsv: <prefix>_airr.tsv, <prefix>_airr_align.tsv, <prefix>_report.tsv, ...
    fasta: list[File]  # *.fa: <prefix>_annot.fa, <prefix>_assembled_reads.fa
    out: list[File]  # *.out: <prefix>_cdr3.out, <prefix>_final.out, <prefix>_raw.out
    fq: list[File]  # *.fq: <prefix>_toassemble.fq, or _toassemble_1.fq and _toassemble_2.fq
    results: Dir  # everything


def file_name(f: File) -> str:
    return f.path.rstrip("/").rsplit("/", 1)[-1]


async def trust4(
    fasta: File,
    prefix: str,
    reads: list[File] | None = None,
    bam: File | None = None,
    single_end: bool = False,
    vdj_reference: File | None = None,
    barcode_whitelist: File | None = None,
    cell_barcode_read: str = "",
    umi_read: str = "",
    read_format: str = "",
    args: str = "",
) -> Trust4Result:
    """TRUST4's TCR/BCR reconstruction from ``bam`` or ``reads`` (one FASTQ if ``single_end``, else R1 and R2)."""
    reads = reads or []
    if bam is None and not reads:
        raise ValueError("trust4 needs a bam or reads")
    if reads and len(reads) != (1 if single_end else 2):
        raise ValueError("trust4 needs one FASTQ single-end, or two (R1, R2) paired-end")
    for option, value in (("cell_barcode_read", cell_barcode_read), ("umi_read", umi_read)):
        if value and value not in READS:
            raise ValueError(f"{option} must be R1 or R2, got {value!r}")
        if value == "R2" and single_end:
            raise ValueError(f"{option} R2 needs paired-end reads")
    results = await trust4_cmd(
        bam=bam,
        reads_1=reads[0] if reads else None,
        reads_2=reads[1] if len(reads) > 1 else None,
        single_end=single_end,
        fasta=fasta,
        vdj_reference=vdj_reference,
        barcode_whitelist=barcode_whitelist,
        cell_barcode_read=cell_barcode_read,
        umi_read=umi_read,
        read_format=read_format,
        prefix=prefix,
        args=args,
    )
    files = sorted(await results.list_files(), key=file_name)

    def ending(suffix: str) -> list[File]:
        return [f for f in files if file_name(f).endswith(suffix)]

    by_name = {file_name(f): f for f in files}
    for name in (f"{prefix}_airr.tsv", f"{prefix}_report.tsv"):
        if name not in by_name:
            raise FileNotFoundError(f"trust4 wrote no {name}")
    return Trust4Result(
        airr=by_name[f"{prefix}_airr.tsv"],
        report=by_name[f"{prefix}_report.tsv"],
        tsv=ending(".tsv"),
        fasta=ending(".fa"),
        out=ending(".out"),
        fq=ending(".fq"),
        results=results,
    )
