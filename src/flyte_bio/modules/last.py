"""last — LAST, pairwise alignment of genomes and reads with trained scoring parameters.

Exposes:

- :func:`last_lastdb` — index a FASTA or FASTQ (plain or gzipped): a Dir of
  ``<prefix>.{bck,des,prj,sds,ssp,suf,tis}``.
- :func:`last_train` — fit substitution and gap parameters for aligning
  ``fastx`` to an index: ``<prefix>.train`` and its MultiQC summary
  ``<prefix>.train.tsv``.
- :func:`last_lastal` — align ``fastx`` to an index, optionally with trained
  parameters: ``<prefix>.maf.gz`` and the MultiQC table ``<prefix>.tsv``.
- :func:`last_split` — keep each query's best alignments (``last-split``):
  ``<prefix>.maf.gz`` and ``<prefix>.tsv``.
- :func:`last_mafconvert` — convert a MAF to another format: ``sam``
  (``<prefix>.sam.gz``), ``bam``, ``cram``, ``bcf`` (variants called with
  ``bcftools``, plus ``<prefix>.stats``), ``gff`` or any other ``maf-convert``
  format (``<prefix>.<format>.gz``; ``blasttab+`` gives
  ``<prefix>.blasttabplus.gz``).
- :func:`last_dotplot` — a dot plot of a MAF: ``<prefix>.<format>`` (``png``
  or ``gif``), optionally after keeping only linked alignments
  (``maf-linked``).

As upstream:

- An index is passed between tasks as a Dir; the tools find its prefix from
  its ``.des`` file. ``lastal`` records the index path in its MAF header, so
  the Dir is linked into the working directory under ``index_dir_name``:
  ``lastdb`` by default, the directory upstream's lastdb writes.
- ``lastdb``, ``last-train`` and ``lastal`` get the task's CPUs as threads
  (``-P``); options in ``args`` come later and win. The thread count can
  change ``lastdb``'s suffix array for many sequences, and ``last-train``
  records it in its output.
- The ``.tsv`` tables (alignment length, percent identity, genome sizes and
  sequence counts) are computed while the MAF streams through, under the
  sample name ``sample_id`` (upstream's ``meta.id``; it defaults to
  ``prefix``).
- SAM, BAM, CRAM and BCF output take a read group ``ID:<sample_id>
  SM:<sample_id>``, and their sequence header from ``dict`` (else from the
  MAF, ``-d``).

``args`` go to the main tool. For :func:`last_mafconvert`, ``args2`` go to
``samtools sort``, ``args3`` to ``bcftools mpileup`` and ``args4`` to
``bcftools call``; for :func:`last_dotplot`, ``args2`` go to ``maf-linked``.
"""

import re
from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

LAST_IMAGE = "community.wave.seqera.io/library/last_gzip:40e65bc1d7d8624e"
# maf-convert's image adds samtools and bcftools for SAM/BAM/CRAM/BCF output.
LAST_MAFCONVERT_IMAGE = "community.wave.seqera.io/library/bcftools_last_samtools_gzip:6baa7d64a57919ca"
# last-dotplot's image adds the font it draws labels with.
LAST_DOTPLOT_IMAGE = "community.wave.seqera.io/library/last_gzip_open-fonts:92dd7f8fc3f0c4fd"

# Upstream labels lastdb and last-split process_medium, last-train, lastal and
# maf-convert process_high, and last-dotplot process_low.
MEDIUM_RESOURCES = flyte.Resources(cpu=4, memory="8Gi")
HIGH_RESOURCES = flyte.Resources(cpu=4, memory="16Gi")
LOW_RESOURCES = flyte.Resources(cpu=1, memory="4Gi")

# Links the index Dir as ./$INDEX_DIR_NAME and sets INDEX to its prefix path,
# $INDEX_DIR_NAME/<name>, found from its .des file.
LINK_INDEX = r"""
        ln -s "$INDEX_DIR" "$INDEX_DIR_NAME"
        DES=("$INDEX_DIR_NAME"/*.des)
        if [ ! -e "${DES[0]}" ]; then echo "no .des file in the index" >&2; exit 1; fi
        INDEX="$INDEX_DIR_NAME/$(basename "${DES[0]}" .des)"
"""

# Shell functions that summarise a MAF stream for MultiQC, as upstream:
# get_genome_stats reads LAST's header and footer counts; calculate_psl_metrics
# reads the PSL conversion. The sample name comes from $SAMPLE.
MAF_STATS_FUNCTIONS = r"""
        get_genome_stats() { awk '
            BEGIN {
                OFS = "\t"
                in_params = 1
            }

            function extract(line, regex,   s) {
                if (match(line, regex)) {
                    s = substr(line, RSTART, RLENGTH)
                    sub(/^[^=]*=/, "", s)
                    return s
                }
                return ""
            }

            /^# Reference sequences=/ {
                in_params = 0
                ref_seq     = extract($0, "^# Reference sequences=[0-9]+")
                ref_letters = extract($0,                "letters=[0-9]+")
            }

            /^#/ {
                if (in_params && gsub(/=/, "&") >= 5) {
                    if (match($0, /s=([0-9]+)/)) {
                        strand_mode = substr($0, RSTART+2, RLENGTH-2)
                    }
                }
            }

            /^# Query sequences=/ {
                qry_seq     = extract($0, "^# Query sequences=[0-9]+")
                qry_letters = extract($0,            "letters=[0-9]+")
            }

            END {
            if (strand_mode == 1) {
                ref_seq     /= 2
                ref_letters /= 2
            }

                print "TargetSequences", "TargetLength", "QuerySequences", "QueryLength"
                print ref_seq+0,         ref_letters+0,  qry_seq+0,         qry_letters+0
            }'
        }

        calculate_psl_metrics() {
            awk -v sample="$SAMPLE" 'BEGIN {
                OFS="\t";
                FS="\t";
                totalMatches = 0;
                totalAlignmentLength = 0;
                totalAlignedBases = 0;
            }
            {
                totalMatches         += $1 +       $3            ;
                totalAlignmentLength += $1 + $2 + $3 + $6 + $8;
                totalAlignedBases    += $1 + $2 + $3            ;
            }
            END {
                percentIdentity       = (totalAlignmentLength > 0) ? (totalMatches / totalAlignmentLength * 100) : 0;
                percentIdentityNoGaps = (totalAlignmentLength > 0) ? (totalMatches / totalAlignedBases    * 100) : 0;
                print "Sample",  "TotalAlignmentLength", "PercentIdentity", "PercentIdentityNoGaps";
                print sample, totalAlignmentLength,   percentIdentity,   percentIdentityNoGaps;
            }'
        }
"""

# After the alignments stream through calculate_psl_metrics into $P.maf: the
# genome counts, the combined MultiQC table $P.tsv, and $P.maf.gz. Upstream
# does the same with process substitutions; reading the written MAF instead
# makes sure each part is complete before the table is pasted together.
SUMMARISE_MAF = r"""
        get_genome_stats < "$P.maf" > "$P.genomestats.txt"
        paste "$P.alignmentstats.txt" "$P.genomestats.txt" > "$P.tsv"
        gzip --no-name "$P.maf"
"""

last_lastdb_cmd = shell.create(
    name="last_lastdb",
    image=LAST_IMAGE,
    resources=MEDIUM_RESOURCES,
    inputs={"fastx": File, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"index": Dir},
    script=rf"""
        {TASK_LIMITS}
        set -e
        FX=({{inputs.fastx}})
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        W=$(mktemp -d); cd "$W"
        ln -s "${{FX[0]}}" "$(basename "${{FX[0]}}")"
        mkdir lastdb
        lastdb $ARGS -P $TASK_CPUS "lastdb/$P" "$(basename "${{FX[0]}}")" || exit 1
        mv lastdb/* {{outputs.index}}/
    """,
)

last_train_cmd = shell.create(
    name="last_train",
    image=LAST_IMAGE,
    resources=HIGH_RESOURCES,
    inputs={"fastx": File, "index": Dir, "index_dir_name": str, "prefix": str, "args": str},
    defaults={"index_dir_name": "lastdb", "args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        set -e
        FX=({{inputs.fastx}})
        INDEX_DIR={{inputs.index}}
        INDEX_DIR_NAME={{inputs.index_dir_name}}
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        W=$(mktemp -d); cd "$W"
        FASTX=$(basename "${{FX[0]}}")
        ln -s "${{FX[0]}}" "$FASTX"
        {LINK_INDEX}
        last-train $ARGS -P $TASK_CPUS "$INDEX" "$FASTX" > "$P.train" || exit 1

        T="$P.train"
        printf 'id\tsubstitution_percent_identity\tlast -t\tlast -a\tlast -A\tlast -b\tlast -B\tlast -S\n' > "$T.tsv"
        printf "$(basename "$T" .target.train)\t" >> "$T.tsv"
        # As upstream, not the last 'substitution percent identity': it comes from a matrix rounded to integers.
        grep 'substitution percent identity' "$T" | tail -n 2 | head -n 1 \
            | awk '{{print $5}}' | tr '\n' '\t' >> "$T.tsv"
        grep 'last -t' "$T" | tail -n 1 | awk '{{print $2}}' | sed -e 's/-t//' | tr '\n' '\t' >> "$T.tsv"
        grep 'last -a' "$T" | tail -n 1 | awk '{{print $3}}' | tr '\n' '\t' >> "$T.tsv"
        grep 'last -A' "$T" | tail -n 1 | awk '{{print $3}}' | tr '\n' '\t' >> "$T.tsv"
        grep 'last -b' "$T" | tail -n 1 | awk '{{print $3}}' | tr '\n' '\t' >> "$T.tsv"
        grep 'last -B' "$T" | tail -n 1 | awk '{{print $3}}' | tr '\n' '\t' >> "$T.tsv"
        grep 'last -S' "$T" | tail -n 1 | awk '{{print $3}}' >> "$T.tsv"
        mv "$T" "$T.tsv" {{outputs.results}}/
    """,
)

last_lastal_cmd = shell.create(
    name="last_lastal",
    image=LAST_IMAGE,
    resources=HIGH_RESOURCES,
    inputs={
        "fastx": File,
        "param_file": File | None,
        "index": Dir,
        "index_dir_name": str,
        "prefix": str,
        "sample_id": str,
        "args": str,
    },
    defaults={"index_dir_name": "lastdb", "args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        set -e
        set -o pipefail
        FX=({{inputs.fastx}})
        shopt -s nullglob; PAR=({{inputs.param_file}}); shopt -u nullglob
        INDEX_DIR={{inputs.index}}
        INDEX_DIR_NAME={{inputs.index_dir_name}}
        P={{inputs.prefix}}
        SAMPLE={{inputs.sample_id}}
        ARGS={{inputs.args}}
        W=$(mktemp -d); cd "$W"
        FASTX=$(basename "${{FX[0]}}")
        ln -s "${{FX[0]}}" "$FASTX"
        TRAINED=()
        if [ ${{#PAR[@]}} -gt 0 ]; then
            ln -s "${{PAR[0]}}" "$(basename "${{PAR[0]}}")"
            TRAINED=(-p "$(basename "${{PAR[0]}}")")
        fi
        {LINK_INDEX}
        {MAF_STATS_FUNCTIONS}
        lastal -P $TASK_CPUS "${{TRAINED[@]}}" $ARGS "$INDEX" "$FASTX" \
            | tee "$P.maf" \
            | maf-convert psl \
            | calculate_psl_metrics > "$P.alignmentstats.txt" || exit 1
        {SUMMARISE_MAF}
        mv "$P.maf.gz" "$P.tsv" {{outputs.results}}/
    """,
)

last_split_cmd = shell.create(
    name="last_split",
    image=LAST_IMAGE,
    resources=MEDIUM_RESOURCES,
    inputs={"maf": File, "prefix": str, "sample_id": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        set -e
        set -o pipefail
        MAF=({{inputs.maf}})
        P={{inputs.prefix}}
        SAMPLE={{inputs.sample_id}}
        ARGS={{inputs.args}}
        W=$(mktemp -d); cd "$W"
        {MAF_STATS_FUNCTIONS}
        zcat < "${{MAF[0]}}" \
            | last-split $ARGS \
            | tee "$P.maf" \
            | maf-convert psl \
            | calculate_psl_metrics > "$P.alignmentstats.txt" || exit 1
        {SUMMARISE_MAF}
        mv "$P.maf.gz" "$P.tsv" {{outputs.results}}/
    """,
)

last_mafconvert_cmd = shell.create(
    name="last_mafconvert",
    image=LAST_MAFCONVERT_IMAGE,
    resources=HIGH_RESOURCES,
    inputs={
        "maf": File,
        "format": str,
        "fasta": File | None,
        "fai": File | None,
        "gzi": File | None,
        "sizes": File | None,
        "dict": File | None,
        "prefix": str,
        "sample_id": str,
        "args": str,
        "args2": str,
        "args3": str,
        "args4": str,
    },
    defaults={"args": "", "args2": "", "args3": "", "args4": ""},
    outputs={"results": Dir},
    script=r"""
        set -e
        set -o pipefail
        MAF=({inputs.maf})
        shopt -s nullglob
        FA=({inputs.fasta}); FAI=({inputs.fai}); GZI=({inputs.gzi}); SIZES=({inputs.sizes}); DICT=({inputs.dict})
        shopt -u nullglob
        FORMAT={inputs.format}
        P={inputs.prefix}
        SAMPLE={inputs.sample_id}
        ARGS={inputs.args}
        ARGS2={inputs.args2}
        ARGS3={inputs.args3}
        ARGS4={inputs.args4}
        W=$(mktemp -d); cd "$W"
        # The FASTA and its indexes side by side, under their own names.
        for f in "${FA[@]}" "${FAI[@]}" "${GZI[@]}" "${SIZES[@]}" "${DICT[@]}" "${MAF[0]}"; do
            ln -s "$f" "$(basename "$f")"
        done
        M=$(basename "${MAF[0]}")
        FASTA=""; if [ ${#FA[@]} -gt 0 ]; then FASTA=$(basename "${FA[0]}"); fi
        REF_ARGS=()
        if [ ${#DICT[@]} -gt 0 ]; then
            D=$(basename "${DICT[0]}")
            DICT_ARGS=(-f "$D")
            if [ "$FORMAT" = cram ]; then
                REF_CRAM=$(grep '^@SQ' "$D" | sed -n 's/.*UR:\([^ \t]*\).*/\1/p' | uniq)
                if [ -r $REF_CRAM ]; then REF_ARGS=(); else REF_ARGS=(--reference "$FASTA"); fi
            fi
        else
            DICT_ARGS=(-d)
        fi
        RG="ID:$SAMPLE SM:$SAMPLE"

        case "$FORMAT" in
            gff)
                {
                    echo "##gff-version 3"
                    if [ ${#SIZES[@]} -gt 0 ]; then
                        awk '{ printf "##sequence-region %s 1 %s\n", $1, $2 }' "$(basename "${SIZES[0]}")"
                    fi
                    maf-convert $ARGS -n gff "$M"
                } | gzip --no-name > "$P.gff.gz" || exit 1
                ;;
            sam)
                maf-convert $ARGS "${DICT_ARGS[@]}" sam "$M" -r "$RG" \
                    | samtools sort -O sam \
                    | gzip --no-name > "$P.sam.gz" || exit 1
                ;;
            bam)
                maf-convert $ARGS "${DICT_ARGS[@]}" sam "$M" -r "$RG" \
                    | samtools sort $ARGS2 -O bam -o "$P.bam" || exit 1
                ;;
            cram)
                # As upstream, CRAM needs a reference that is plain or bgzipped.
                maf-convert $ARGS "${DICT_ARGS[@]}" sam "$M" -r "$RG" \
                    | samtools sort $ARGS2 -O cram "${REF_ARGS[@]}" -o "$P.cram" || exit 1
                ;;
            bcf)
                maf-convert $ARGS "${DICT_ARGS[@]}" sam "$M" -r "$RG" \
                    | samtools sort $ARGS2 -u \
                    | bcftools mpileup $ARGS3 --fasta-ref "$FASTA" -Ou - \
                    | bcftools call $ARGS4 -Ob -o "$P.bcf" || exit 1
                bcftools stats "$P.bcf" > "$P.stats" || exit 1
                ;;
            blasttab+)
                maf-convert $ARGS "$FORMAT" "$M" | gzip --no-name > "$P.blasttabplus.gz" || exit 1
                ;;
            *)
                maf-convert $ARGS "$FORMAT" "$M" | gzip --no-name > "$P.$FORMAT.gz" || exit 1
                ;;
        esac
        rm -f "$M"
        for f in "$P".*; do
            if [ -f "$f" ] && [ ! -L "$f" ]; then mv "$f" {outputs.results}/; fi
        done
    """,
)

last_dotplot_cmd = shell.create(
    name="last_dotplot",
    image=LAST_DOTPLOT_IMAGE,
    resources=LOW_RESOURCES,
    inputs={
        "maf": File,
        "annot_a": File | None,
        "annot_b": File | None,
        "format": str,
        "filter": bool,
        "prefix": str,
        "args": str,
        "args2": str,
    },
    defaults={"args": "", "args2": ""},
    outputs={"results": Dir},
    script=r"""
        set -e
        set -o pipefail
        MAF=({inputs.maf})
        shopt -s nullglob; ANNOT_A=({inputs.annot_a}); ANNOT_B=({inputs.annot_b}); shopt -u nullglob
        FORMAT={inputs.format}
        P={inputs.prefix}
        ARGS={inputs.args}
        ARGS2={inputs.args2}
        W=$(mktemp -d); cd "$W"
        ANNOT=()
        if [ ${#ANNOT_A[@]} -gt 0 ]; then
            ln -s "${ANNOT_A[0]}" "$(basename "${ANNOT_A[0]}")"; ANNOT+=(-a "$(basename "${ANNOT_A[0]}")")
        fi
        if [ ${#ANNOT_B[@]} -gt 0 ]; then
            ln -s "${ANNOT_B[0]}" "$(basename "${ANNOT_B[0]}")"; ANNOT+=(-b "$(basename "${ANNOT_B[0]}")")
        fi
        TTF=/home/runner/conda_pkgs_dir/open-fonts-0.7.0-1/fonts/open-fonts/DejaVuSansMono-Regular.ttf
        [ -e "$TTF" ] || TTF="/opt/conda/fonts/open-fonts/DejaVuSansMono-Regular.ttf"
        if [ {inputs.filter} = true ]; then
            READ=(maf-linked $ARGS2)
        else
            READ=(zcat -f)
        fi
        "${READ[@]}" "${MAF[0]}" \
            | last-dotplot -f "$TTF" $ARGS "${ANNOT[@]}" - "$P.$FORMAT" || exit 1
        mv "$P.$FORMAT" {outputs.results}/
    """,
)

# One image per environment: maf-convert's adds samtools and bcftools, and
# last-dotplot's a font, so `env` aggregates three.
last_tools_env = flyte.TaskEnvironment.from_task(
    "last_tools",
    last_lastdb_cmd.as_task(),
    last_train_cmd.as_task(),
    last_lastal_cmd.as_task(),
    last_split_cmd.as_task(),
)
last_mafconvert_env = flyte.TaskEnvironment.from_task("last_mafconvert", last_mafconvert_cmd.as_task())
last_dotplot_env = flyte.TaskEnvironment.from_task("last_dotplot", last_dotplot_cmd.as_task())
env = flyte.TaskEnvironment(name="last", depends_on=[last_tools_env, last_mafconvert_env, last_dotplot_env])


@dataclass
class LastTrainResult:
    param_file: File  # <prefix>.train, for lastal -p
    multiqc: File  # <prefix>.train.tsv
    results: Dir


@dataclass
class LastAlignmentResult:
    maf: File  # <prefix>.maf.gz
    multiqc: File  # <prefix>.tsv
    results: Dir


@dataclass
class LastMafconvertResult:
    alignment: File  # <prefix>.<format>[.gz]
    index: File | None  # <prefix>.<format>.{bai,crai,csi}, when args2 asks samtools sort for one
    stats: File | None  # <prefix>.stats, for bcf
    results: Dir


async def pick(results: Dir, name: str, tool: str) -> File:
    found = await results.get_file(name)
    if found is None:
        raise FileNotFoundError(f"{tool} wrote no {name}")
    return found


async def last_lastdb(fastx: File, prefix: str = "sample", args: str = "") -> Dir:
    """LAST index of ``fastx``: a Dir of ``<prefix>.{bck,des,prj,sds,ssp,suf,tis}``."""
    return await last_lastdb_cmd(fastx=fastx, prefix=prefix, args=args)


async def last_train(
    fastx: File, index: Dir, prefix: str = "sample", args: str = "", index_dir_name: str = "lastdb"
) -> LastTrainResult:
    """Train alignment parameters for ``fastx`` against ``index`` (from :func:`last_lastdb`)."""
    results = await last_train_cmd(fastx=fastx, index=index, index_dir_name=index_dir_name, prefix=prefix, args=args)
    return LastTrainResult(
        param_file=await pick(results, f"{prefix}.train", "last-train"),
        multiqc=await pick(results, f"{prefix}.train.tsv", "last-train"),
        results=results,
    )


async def last_lastal(
    fastx: File,
    index: Dir,
    param_file: File | None = None,
    prefix: str = "sample",
    sample_id: str | None = None,
    args: str = "",
    index_dir_name: str = "lastdb",
) -> LastAlignmentResult:
    """Align ``fastx`` to ``index``, with ``param_file`` (from :func:`last_train`) when given."""
    results = await last_lastal_cmd(
        fastx=fastx,
        param_file=param_file,
        index=index,
        index_dir_name=index_dir_name,
        prefix=prefix,
        sample_id=sample_id or prefix,
        args=args,
    )
    return LastAlignmentResult(
        maf=await pick(results, f"{prefix}.maf.gz", "lastal"),
        multiqc=await pick(results, f"{prefix}.tsv", "lastal"),
        results=results,
    )


async def last_split(
    maf: File, prefix: str = "sample", sample_id: str | None = None, args: str = ""
) -> LastAlignmentResult:
    """Keep each query's best alignments in ``maf`` (gzipped), as ``<prefix>.maf.gz``."""
    name = maf.name or maf.path.rsplit("/", 1)[-1]
    if name == f"{prefix}.maf.gz":
        raise ValueError("last_split: input and output names are the same; pass another prefix")
    results = await last_split_cmd(maf=maf, prefix=prefix, sample_id=sample_id or prefix, args=args)
    return LastAlignmentResult(
        maf=await pick(results, f"{prefix}.maf.gz", "last-split"),
        multiqc=await pick(results, f"{prefix}.tsv", "last-split"),
        results=results,
    )


def mafconvert_output_name(prefix: str, format: str) -> str:
    if format in ("bam", "cram", "bcf"):
        return f"{prefix}.{format}"
    if format == "blasttab+":
        return f"{prefix}.blasttabplus.gz"
    return f"{prefix}.{format}.gz"


async def last_mafconvert(
    maf: File,
    format: str,
    fasta: File | None = None,
    fai: File | None = None,
    gzi: File | None = None,
    sizes: File | None = None,
    sequence_dict: File | None = None,
    prefix: str = "sample",
    sample_id: str | None = None,
    args: str = "",
    args2: str = "",
    args3: str = "",
    args4: str = "",
) -> LastMafconvertResult:
    """Convert ``maf`` to ``format``; ``bcf`` needs ``fasta`` (with ``fai``), ``gff`` uses ``sizes``.

    ``sequence_dict`` is the reference's ``.dict`` (upstream's ``dict``).
    """
    if format == "bcf":
        if fasta is None:
            raise ValueError("last_mafconvert: bcf output needs a fasta")
        # As upstream: samtools sort writes to a pipe here, where it can't write an index.
        args2 = re.sub(r"\s*--write-index\b", "", args2)
    results = await last_mafconvert_cmd(
        maf=maf,
        format=format,
        fasta=fasta,
        fai=fai,
        gzi=gzi,
        sizes=sizes,
        dict=sequence_dict,
        prefix=prefix,
        sample_id=sample_id or prefix,
        args=args,
        args2=args2,
        args3=args3,
        args4=args4,
    )
    name = mafconvert_output_name(prefix, format)
    index = None
    for ext in ("bai", "crai", "csi"):
        index = index or await results.get_file(f"{name}.{ext}")
    return LastMafconvertResult(
        alignment=await pick(results, name, "maf-convert"),
        index=index,
        stats=await results.get_file(f"{prefix}.stats") if format == "bcf" else None,
        results=results,
    )


async def last_dotplot(
    maf: File,
    annot_a: File | None = None,
    annot_b: File | None = None,
    format: str = "png",
    filter: bool = False,
    prefix: str = "sample",
    args: str = "",
    args2: str = "",
) -> File:
    """Dot plot of ``maf`` as ``<prefix>.<format>``; ``filter`` keeps only linked alignments first.

    ``annot_a`` and ``annot_b`` are BED-like annotations of the reference and
    query sequences.
    """
    results = await last_dotplot_cmd(
        maf=maf,
        annot_a=annot_a,
        annot_b=annot_b,
        format=format,
        filter=filter,
        prefix=prefix,
        args=args,
        args2=args2,
    )
    return await pick(results, f"{prefix}.{format}", "last-dotplot")
