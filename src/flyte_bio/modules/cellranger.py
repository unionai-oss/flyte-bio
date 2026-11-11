"""cellranger — 10x Genomics Cell Ranger: references, gene expression, multi and V(D)J.

Exposes:

- :func:`cellranger_mkgtf`, which filters a GTF by attributes
  (``cellranger mkgtf``; pass ``--attribute=gene_biotype:...`` in ``args``)
  into ``<prefix>.gtf``.
- :func:`cellranger_mkref`, which builds a Cell Ranger transcriptome
  reference (``cellranger mkref``: STAR index, ``fasta/``, ``genes/``,
  ``reference.json``) from a genome FASTA and a GTF, returned as
  ``<reference_name>.tar.gz`` (an archive of the ``<reference_name>``
  folder).
- :func:`cellranger_count`, which quantifies one sample's gene expression
  (``cellranger count``) against such a reference, extracting it first. As
  upstream's template does, the FASTQs are taken as R1/R2 pairs in file-name
  order and renamed ``<sample_id>_S1_L00<n>_R{1,2}_001.fastq.gz`` before the
  run. Cell Ranger's ``<prefix>/outs`` folder comes back as
  ``<prefix>.outs.tar.gz``, with its main files alongside
  (``metrics_summary.csv``, ``web_summary.html``, the filtered and raw
  ``*_feature_bc_matrix.h5``, ``molecule_info.h5`` and, with
  ``--create-bam true``, ``possorted_genome_bam.bam``).
- :func:`cellranger_multi`, which runs ``cellranger multi`` on one sample's
  libraries (gene expression, V(D)J, antibody, antigen, multiplexing, CRISPR)
  with their references. As upstream, the multi config CSV is built from what
  is given (:func:`cellranger_multi_config`), the FASTQs of each library are
  renamed ``<prefix>_S1_L00<n>_R{1,2}_001.fastq.gz`` (unless
  ``skip_renaming``), and the feature reference is copied before use. It
  returns the config (``cellranger_multi_config.csv``) and ``<prefix>/outs``
  as ``<prefix>.outs.tar.gz``.
- :func:`cellranger_mkfastq`, which demultiplexes an Illumina run folder,
  given as a ``.tar.gz`` archive, with a simple or Illumina sample sheet
  (``cellranger mkfastq``, which drives bcl2fastq). It returns
  ``fastq.tar.gz`` (every FASTQ under ``outs/fastq_path``, with its path
  there: ``<project>/<sample>/<sample>_S1_L001_R1_001.fastq.gz``, index reads
  ``..._I1_001.fastq.gz``, ``Undetermined_S0_...``), ``Reports.tar.gz``,
  ``Stats.tar.gz`` and the ``outs/interop_path/*.bin`` files.
- :func:`cellranger_mkvdjref`, which builds a Cell Ranger V(D)J reference
  (``cellranger mkvdjref``) from a genome FASTA and a GTF, or from
  ready-made V(D)J segment sequences (``seqs``). It returns the
  ``<reference_name>`` folder (``reference.json``, ``fasta/regions.fa``, ...)
  as ``<reference_name>.tar.gz``, plus its ``regions.fa`` and
  ``reference.json``.
- :func:`cellranger_vdj`, which assembles one sample's V(D)J transcripts and
  calls clonotypes (``cellranger vdj``) against such a reference, extracting
  it first. ``reads`` are the sample's FASTQs, named as ``bcl2fastq`` names
  them (``<sample>_S1_L001_R1_001.fastq.gz``). It returns Cell Ranger's
  ``<prefix>/outs`` folder (``airr_rearrangement.tsv``,
  ``filtered_contig_annotations.csv``, ``clonotypes.csv``,
  ``metrics_summary.csv``, ``web_summary.html``, ``vdj_reference/``, ...),
  and the same as ``<prefix>.outs.tar.gz`` (an archive of ``<prefix>/outs``)
  to pass on whole.

References go between tasks as ``.tar.gz`` archives of the reference folder
(a folder with ``reference.json`` at the top, or a single such folder). The
run commands pass the task's cores and memory as ``--localcores`` and
``--localmem`` (and mkref's ``--nthreads``). Telemetry is disabled
(``TENX_DISABLE_TELEMETRY``) and ``HOME`` is the task's working directory.

The containers are upstream's builds of 10x's Cell Ranger (10.1.0; 9.0.1 for
mkfastq, which also bundles Illumina's bcl2fastq), distributed under 10x
Genomics' licence: the software may be used only in connection with 10x
Genomics products or data generated with them, solely for internal use, and
not redistributed or sublicensed. bcl2fastq is under Illumina's end-user
licence.
"""

from dataclasses import dataclass, field

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

CELLRANGER_IMAGE = "quay.io/nf-core/cellranger:10.1.0"
# Upstream's mkfastq image: Cell Ranger 9.0.1 with Illumina's bcl2fastq.
CELLRANGER_MKFASTQ_IMAGE = "quay.io/nf-core/cellrangermkfastq:9.0.1"
NO_TELEMETRY = {"TENX_DISABLE_TELEMETRY": "1"}
# Upstream's process_high asks for 12 CPUs and 72 GB; its tests run in 6 GB.
# Real samples want 16-64 GiB (10x recommends 64 GB or more for vdj): raise
# with <cmd>.override(resources=...).
MKVDJREF_RESOURCES = flyte.Resources(cpu=4, memory="8Gi", disk="20Gi")
VDJ_RESOURCES = flyte.Resources(cpu=4, memory="16Gi", disk="50Gi")
# Upstream: process_low.
MKGTF_RESOURCES = flyte.Resources(cpu=1, memory="4Gi", disk="10Gi")
# Upstream: process_high. A human genome needs about 32 GiB for STAR's index.
MKREF_RESOURCES = flyte.Resources(cpu=4, memory="16Gi", disk="50Gi")
COUNT_RESOURCES = flyte.Resources(cpu=4, memory="16Gi", disk="50Gi")
MULTI_RESOURCES = flyte.Resources(cpu=4, memory="16Gi", disk="50Gi")
# Upstream: process_medium.
MKFASTQ_RESOURCES = flyte.Resources(cpu=4, memory="16Gi", disk="50Gi")

# A bash function: extract the reference archive $1 and move its reference
# folder (the archive's top, or its single top-level folder) to $2.
EXTRACT_REFERENCE = r"""
        extract_reference() {
            local TMP
            TMP=$(mktemp -d -p .)
            tar -xzf "$1" -C "$TMP" || exit 1
            if [ -f "$TMP/reference.json" ]; then
                mv "$TMP" "$2"
            else
                local INNER=("$TMP"/*/)
                mv "${INNER[0]%/}" "$2"
                rm -rf "$TMP"
            fi
        }
""".strip()

# Upstream's template, unchanged but for taking its templated values (sample
# id, prefix, reference, cores, memory in GB, args) as command-line arguments,
# and Groovy's escaping undone.
CELLRANGER_COUNT_PY = r'''
#!/usr/bin/env python3
"""
Automatically rename staged files for input into cellranger count.

Copyright (c) Gregor Sturm 2023 - MIT License
"""

import re
import shlex
import sys
from pathlib import Path
from subprocess import run
from textwrap import dedent

sample_id, prefix, reference, cpus, memory, args = sys.argv[1:7]


def chunk_iter(seq, size):
    """iterate over `seq` in chunks of `size`"""
    return (seq[pos : pos + size] for pos in range(0, len(seq), size))


# get fastqs, ordered by path. Files are staged into
#   - "fastq_001/{original_name.fastq.gz}"
#   - "fastq_002/{original_name.fastq.gz}"
#   - ...
# Since we require fastq files in the input channel to be ordered such that a R1/R2 pair
# of files follows each other, ordering will get us a sequence of [R1, R2, R1, R2, ...]
fastqs = sorted(Path(".").glob("fastq_*/*"))
assert len(fastqs) % 2 == 0

# target directory in which the renamed fastqs will be placed
fastq_all = Path("./fastq_all")
fastq_all.mkdir(exist_ok=True)

# Match R1 in the filename, but only if it is followed by a non-digit or non-character
# match "file_R1.fastq.gz", "file.R1_000.fastq.gz", etc. but
# do not match "SRR12345", "file_INFIXR12", etc
filename_pattern = r"([^a-zA-Z0-9])R1([^a-zA-Z0-9])"

for i, (r1, r2) in enumerate(chunk_iter(fastqs, 2), start=1):
    if re.sub(filename_pattern, r"\1R2\2", r1.name) != r2.name:
        raise AssertionError(
            dedent(
                f"""\
                We expect R1 and R2 of the same sample to have the same filename except for R1/R2.
                This has been checked by replacing "R1" with "R2" in the first filename and comparing it to the second filename.
                If you believe this check shouldn't have failed on your filenames, please report an issue on GitHub!

                Files involved:
                    - {r1}
                    - {r2}
                """
            )
        )
    r1.rename(fastq_all / f"{sample_id}_S1_L{i:03d}_R1_001.fastq.gz")
    r2.rename(fastq_all / f"{sample_id}_S1_L{i:03d}_R2_001.fastq.gz")

# fmt: off
run(
    [
        "cellranger", "count",
        "--id", prefix,
        "--fastqs", str(fastq_all),
        "--transcriptome", reference,
        "--localcores", cpus,
        "--localmem", memory,
        *shlex.split(args),
    ],
    check=True,
)
# fmt: on

# Output version information
version = run(
    ["cellranger", "-V"],
    text=True,
    check=True,
    capture_output=True,
).stdout.replace("cellranger cellranger-", "")

# alas, no `pyyaml` pre-installed in the cellranger container
with open("versions.yml", "w") as f:
    f.write('"CELLRANGER_COUNT":\n')
    f.write(f'  cellranger: "{version}"\n')
'''.strip()  # noqa: E501

cellranger_mkgtf_cmd = shell.create(
    name="cellranger_mkgtf",
    image=CELLRANGER_IMAGE,
    resources=MKGTF_RESOURCES,
    env_vars=NO_TELEMETRY,
    inputs={"gtf": File, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        GTF=({inputs.gtf})
        P={inputs.prefix}
        ARGS={inputs.args}
        OUT={outputs.results}
        W=$(mktemp -d); cd "$W"
        export HOME="$W"
        cellranger \
            mkgtf \
            "${GTF[0]}" \
            "$P.gtf" \
            $ARGS > run.log 2>&1 || { tail -c 3000 run.log >&2; exit 1; }
        mv "$P.gtf" "$OUT"/
    """,
)

cellranger_mkref_cmd = shell.create(
    name="cellranger_mkref",
    image=CELLRANGER_IMAGE,
    resources=MKREF_RESOURCES,
    env_vars=NO_TELEMETRY,
    inputs={"fasta": File, "gtf": File, "reference_name": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        FA=({{inputs.fasta}}); GTF=({{inputs.gtf}})
        NAME={{inputs.reference_name}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        export HOME="$W"
        cellranger \
            mkref \
            --genome="$NAME" \
            --fasta="${{FA[0]}}" \
            --genes="${{GTF[0]}}" \
            --localcores=$TASK_CPUS \
            --localmem=$(( TASK_MEMORY_MB / 1024 )) \
            --nthreads=$TASK_CPUS \
            $ARGS > run.log 2>&1 || {{ tail -c 3000 run.log >&2; exit 1; }}
        tar -czf "$OUT/$NAME.tar.gz" "$NAME"
        cp "$NAME/reference.json" "$OUT"/
    """,
)

cellranger_count_cmd = shell.create(
    name="cellranger_count",
    image=CELLRANGER_IMAGE,
    resources=COUNT_RESOURCES,
    env_vars=NO_TELEMETRY,
    inputs={"reads": list[File], "reference": File, "sample_id": str, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        {EXTRACT_REFERENCE}
        READS=({{inputs.reads}})
        REF_ARCHIVE=({{inputs.reference}})
        ID={{inputs.sample_id}}
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        export HOME="$W"

        extract_reference "${{REF_ARCHIVE[0]}}" reference
        I=0
        for R in "${{READS[@]}}"; do
            I=$((I + 1))
            D=fastq_$(printf %03d $I)
            mkdir "$D"; ln -s "$R" "$D"/
        done
        cat > cellranger_count.py <<'END_OF_SCRIPT'
{CELLRANGER_COUNT_PY}
END_OF_SCRIPT
        echo "$ARGS"
        python3 cellranger_count.py "$ID" "$P" reference "$TASK_CPUS" "$(( TASK_MEMORY_MB / 1024 ))" "$ARGS" \
            > run.log 2>&1 || {{ tail -c 3000 run.log >&2; exit 1; }}

        tar -czhf "$OUT/$P.outs.tar.gz" "$P/outs"
        for F in metrics_summary.csv web_summary.html filtered_feature_bc_matrix.h5 raw_feature_bc_matrix.h5 \
                molecule_info.h5 possorted_genome_bam.bam possorted_genome_bam.bam.bai; do
            if [ -e "$P/outs/$F" ]; then cp -L "$P/outs/$F" "$OUT"/; fi
        done
        cp versions.yml "$OUT"/
    """,
)

MULTI_MODALITIES = ("gex", "vdj", "vdj_t", "vdj_t_gd", "vdj_b", "ab", "beam", "cmo", "crispr")
# Where upstream stages each single-file reference (as references/.../<its name>).
MULTI_REFERENCE_DIRS = {
    "gex_frna_probeset": "references/gex/probeset",
    "gex_targetpanel": "references/gex/targetpanel",
    "vdj_primer_index": "references/vdj/primers",
    "fb_reference": "references/fb",
    "beam_antigen_panel": "references/beam/panel/antigens",
    "beam_control_panel": "references/beam/panel/controls",
    "cmo_reference": "references/cmo",
    "cmo_barcodes": "references/cmo/barcodes",
    "cmo_barcode_assignment": "references/cmo/sample_barcode_assignment",
    "frna_sampleinfo": "references/frna",
    "ocm_barcodes": "references/ocm/barcodes",
}
MULTI_REFERENCE_STAGING = "\n".join(
    f"            {folder} {{inputs.{kind}}} -" for kind, folder in MULTI_REFERENCE_DIRS.items()
)

cellranger_multi_cmd = shell.create(
    name="cellranger_multi",
    image=CELLRANGER_IMAGE,
    resources=MULTI_RESOURCES,
    env_vars=NO_TELEMETRY,
    inputs={
        **{f"{m}_reads": list[File] | None for m in MULTI_MODALITIES},
        "gex_reference": File | None,
        "vdj_reference": File | None,
        **dict.fromkeys(MULTI_REFERENCE_DIRS, File | None),
        "config": str,
        "prefix": str,
        "skip_renaming": str,
        "args": str,
    },
    defaults={"args": "", "skip_renaming": "false"},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        {EXTRACT_REFERENCE}
        shopt -s nullglob
        GEX=({{inputs.gex_reads}}); VDJ=({{inputs.vdj_reads}}); VDJ_T=({{inputs.vdj_t_reads}})
        VDJ_T_GD=({{inputs.vdj_t_gd_reads}}); VDJ_B=({{inputs.vdj_b_reads}}); AB=({{inputs.ab_reads}})
        BEAM=({{inputs.beam_reads}}); CMO=({{inputs.cmo_reads}}); CRISPR=({{inputs.crispr_reads}})
        GEX_REF=({{inputs.gex_reference}}); VDJ_REF=({{inputs.vdj_reference}})
        REFERENCE_FILES=(
{MULTI_REFERENCE_STAGING}
        )
        FB=({{inputs.fb_reference}})
        shopt -u nullglob
        prefix={{inputs.prefix}}
        CONFIG={{inputs.config}}
        skip_renaming={{inputs.skip_renaming}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        export HOME="$W"

        # Stage each library's FASTQs as fastqs/<modality>/fastq_<nnn>/<name>, in file-name order.
        stage_fastqs() {{
            local MODALITY=$1; shift
            local I=0
            for F in "$@"; do
                I=$((I + 1))
                local D=fastqs/$MODALITY/fastq_$(printf %03d $I)
                mkdir -p "$D"; ln -s "$F" "$D"/
            done
        }}
        stage_fastqs gex "${{GEX[@]}}"; stage_fastqs vdj "${{VDJ[@]}}"; stage_fastqs vdj_t "${{VDJ_T[@]}}"
        stage_fastqs vdj_t_gd "${{VDJ_T_GD[@]}}"; stage_fastqs vdj_b "${{VDJ_B[@]}}"; stage_fastqs ab "${{AB[@]}}"
        stage_fastqs beam "${{BEAM[@]}}"; stage_fastqs cmo "${{CMO[@]}}"; stage_fastqs crispr "${{CRISPR[@]}}"

        # References where upstream stages them: references/{{gex,vdj}}/<archive name without .tar.gz>, and
        # files as references/<folder>/<name>.
        if [ ${{#GEX_REF[@]}} -gt 0 ]; then
            mkdir -p references/gex
            extract_reference "${{GEX_REF[0]}}" "references/gex/$(basename "${{GEX_REF[0]}}" .tar.gz)"
        fi
        if [ ${{#VDJ_REF[@]}} -gt 0 ]; then
            mkdir -p references/vdj
            extract_reference "${{VDJ_REF[0]}}" "references/vdj/$(basename "${{VDJ_REF[0]}}" .tar.gz)"
        fi
        FOLDER=""
        for ITEM in "${{REFERENCE_FILES[@]}}"; do
            if [ "$ITEM" = - ]; then FOLDER=""
            elif [ -z "$FOLDER" ]; then FOLDER=$ITEM
            else mkdir -p "$FOLDER"; ln -s "$ITEM" "$FOLDER"/
            fi
        done

        #
        # Symlink FASTQs into fastq_all/, maintaining R1/R2 lane-pairing order.
        # skip_renaming=false (default): rename to Cell Ranger convention ${{prefix}}_S1_L00N_R[12]_001.fastq.gz
        # skip_renaming=true:           keep original filenames as-is
        #
        mkdir -p fastq_all/{{gex,vdj,vdj_t,vdj_t_gd,vdj_b,ab,beam,cmo,crispr}}

        for modality in gex vdj vdj_t vdj_t_gd vdj_b ab beam cmo crispr; do
            # absent modalities are not staged; skip them so find does not fail under pipefail
            [ -d "fastqs/${{modality}}" ] || continue
            lane=1
            n_fastq_dirs=$(find fastqs/${{modality}} -maxdepth 1 -type d -name "fastq_*" | wc -l)
            if [ $((n_fastq_dirs % 2)) -ne 0 ]; then
                echo "ERROR: Found an odd number (${{n_fastq_dirs}}) of staged FASTQ files for modality '${{modality}}'. Expected R1/R2 pairs." >&2
                exit 1
            fi
            while IFS= read -r -d '' r1_dir && IFS= read -r -d '' r2_dir; do
                if [ "${{skip_renaming}}" = "true" ]; then
                    r1=$(find "${{r1_dir}}" -maxdepth 1 -name "*.fastq.gz" | head -1)
                    r2=$(find "${{r2_dir}}" -maxdepth 1 -name "*.fastq.gz" | head -1)
                    [ -z "${{r1}}" ] || [ -z "${{r2}}" ] && continue
                    ln -sf "$(readlink -f "${{r1}}")" "fastq_all/${{modality}}/$(basename "${{r1}}")"
                    ln -sf "$(readlink -f "${{r2}}")" "fastq_all/${{modality}}/$(basename "${{r2}}")"
                else
                    r1=$(find "${{r1_dir}}" -maxdepth 1 -name "*_R1_*.fastq.gz" | head -1)
                    r2=$(find "${{r2_dir}}" -maxdepth 1 -name "*_R2_*.fastq.gz" | head -1)
                    [ -z "${{r1}}" ] || [ -z "${{r2}}" ] && continue

                    r1_base="$(basename "${{r1}}")"
                    r2_base="$(basename "${{r2}}")"
                    if [ "${{r1_base/_R1_/_R2_}}" != "${{r2_base}}" ]; then
                        echo "ERROR: R1 and R2 basenames do not match except for R1/R2 replacement." >&2
                        echo "       R1: ${{r1_base}}" >&2
                        echo "       R2: ${{r2_base}}" >&2
                        exit 1
                    fi

                    ln -sf "$(readlink -f "${{r1}}")" "fastq_all/${{modality}}/${{prefix}}_S1_L$(printf %03d ${{lane}})_R1_001.fastq.gz"
                    ln -sf "$(readlink -f "${{r2}}")" "fastq_all/${{modality}}/${{prefix}}_S1_L$(printf %03d ${{lane}})_R2_001.fastq.gz"
                fi
                lane=$((lane + 1))
            done < <(find fastqs/${{modality}} -maxdepth 1 -type d -name "fastq_*" | sort | xargs -n1 printf '%s\0')
        done

        #
        # Copy fb_reference to avoid symlink corruption
        # Cell Ranger writes to this file during validation, which corrupts the symlinked original
        #
        if [ ${{#FB[@]}} -gt 0 ]; then
            cp "${{FB[0]}}" "fb_reference_copy.csv"
        fi

        # The config as upstream's unquoted here-document: $PWD and $(cat ...) expand.
        eval "cat <<CONFIG_EOF"$'\n'"$CONFIG"$'\n'"CONFIG_EOF" > cellranger_multi_config.csv

        cellranger multi \
            --id="$prefix" \
            --csv=cellranger_multi_config.csv \
            --localcores=$TASK_CPUS \
            --localmem=$(( TASK_MEMORY_MB / 1024 )) \
            $ARGS > run.log 2>&1 || {{ tail -c 3000 run.log >&2; tail -c 2000 run.log > /dev/termination-log 2> /dev/null; exit 1; }}

        cp cellranger_multi_config.csv "$OUT"/
        tar -czhf "$OUT/$prefix.outs.tar.gz" "$prefix/outs"
    """,  # noqa: E501
)

cellranger_mkfastq_cmd = shell.create(
    name="cellranger_mkfastq",
    image=CELLRANGER_MKFASTQ_IMAGE,
    resources=MKFASTQ_RESOURCES,
    env_vars=NO_TELEMETRY,
    inputs={"samplesheet": File, "run_dir": File, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        SHEET=({{inputs.samplesheet}})
        RUN=({{inputs.run_dir}})
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        export HOME="$W"

        # The run folder: the archive's top, or its single top-level folder.
        mkdir run
        if [[ $(tar -tzf "${{RUN[0]}}" | grep -o -P "^.*?\/" | uniq | wc -l) -eq 1 ]]; then
            tar -C run --strip-components 1 -xzf "${{RUN[0]}}"
        else
            tar -C run -xzf "${{RUN[0]}}"
        fi

        cellranger \
            mkfastq \
            --id="${{P}}_outs" \
            --run="$W/run" \
            --csv="${{SHEET[0]}}" \
            --localcores=$TASK_CPUS \
            --localmem=$(( TASK_MEMORY_MB / 1024 )) \
            $ARGS > run.log 2>&1 || {{ tail -c 3000 run.log >&2; exit 1; }}

        FASTQ_PATH="${{P}}_outs/outs/fastq_path"
        (cd "$FASTQ_PATH" && find . -name '*.fastq.gz' | sort | tar -czf "$OUT/fastq.tar.gz" -T -)
        if [ -d "$FASTQ_PATH/Reports" ]; then tar -C "$FASTQ_PATH" -czf "$OUT/Reports.tar.gz" Reports; fi
        if [ -d "$FASTQ_PATH/Stats" ]; then tar -C "$FASTQ_PATH" -czf "$OUT/Stats.tar.gz" Stats; fi
        mkdir "$OUT/interop"
        cp "${{P}}_outs"/outs/interop_path/*.bin "$OUT/interop"/
    """,
)

cellranger_mkvdjref_cmd = shell.create(
    name="cellranger_mkvdjref",
    image=CELLRANGER_IMAGE,
    resources=MKVDJREF_RESOURCES,
    env_vars=NO_TELEMETRY,
    inputs={
        "fasta": File | None,
        "gtf": File | None,
        "seqs": File | None,
        "reference_name": str,
        "args": str,
    },
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        shopt -s nullglob
        FA=({{inputs.fasta}}); GTF=({{inputs.gtf}}); SEQS=({{inputs.seqs}})
        shopt -u nullglob
        NAME={{inputs.reference_name}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        export HOME="$W"

        EXTRA=()
        if [ ${{#GTF[@]}} -gt 0 ]; then EXTRA+=(--genes "${{GTF[0]}}"); fi
        if [ ${{#FA[@]}} -gt 0 ]; then EXTRA+=(--fasta "${{FA[0]}}"); fi
        if [ ${{#SEQS[@]}} -gt 0 ]; then EXTRA+=(--seqs "${{SEQS[0]}}"); fi

        cellranger \
            mkvdjref \
            --genome="$NAME" \
            "${{EXTRA[@]}}" \
            --localcores=$TASK_CPUS \
            --localmem=$(( TASK_MEMORY_MB / 1024 )) \
            $ARGS || exit 1

        tar -czf "$OUT/$NAME.tar.gz" "$NAME"
        cp "$NAME/fasta/regions.fa" "$NAME/reference.json" "$OUT"/
    """,
)

cellranger_vdj_cmd = shell.create(
    name="cellranger_vdj",
    image=CELLRANGER_IMAGE,
    resources=VDJ_RESOURCES,
    env_vars=NO_TELEMETRY,
    inputs={"reads": list[File], "reference": File, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        READS=({{inputs.reads}})
        REF_ARCHIVE=({{inputs.reference}})
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        export HOME="$W"

        mkdir reference fastqs
        tar -xzf "${{REF_ARCHIVE[0]}}" -C reference
        if [ -f reference/reference.json ]; then REF="$W/reference"; else REF=(reference/*/); REF="$W/${{REF[0]%/}}"; fi
        for R in "${{READS[@]}}"; do ln -s "$R" fastqs/; done

        cellranger \
            vdj \
            --id="$P" \
            --fastqs="$W/fastqs" \
            --reference="$REF" \
            --localcores=$TASK_CPUS \
            --localmem=$(( TASK_MEMORY_MB / 1024 )) \
            $ARGS || exit 1

        tar -czhf "$OUT/$P.outs.tar.gz" "$P/outs"
        mkdir "$OUT/outs"
        cp -rL "$P"/outs/. "$OUT/outs"/
    """,
)

cellranger_main_env = flyte.TaskEnvironment.from_task(
    "cellranger_main",
    cellranger_mkgtf_cmd.as_task(),
    cellranger_mkref_cmd.as_task(),
    cellranger_count_cmd.as_task(),
    cellranger_multi_cmd.as_task(),
    cellranger_mkvdjref_cmd.as_task(),
    cellranger_vdj_cmd.as_task(),
)
cellranger_mkfastq_env = flyte.TaskEnvironment.from_task("cellranger_mkfastq", cellranger_mkfastq_cmd.as_task())
env = flyte.TaskEnvironment(name="cellranger", depends_on=[cellranger_main_env, cellranger_mkfastq_env])


async def pick_file(results: Dir, name: str, tool: str) -> File:
    found = await results.get_file(name)
    if found is None:
        raise FileNotFoundError(f"{tool} wrote no {name}")
    return found


async def cellranger_mkgtf(gtf: File, prefix: str, args: str = "") -> File:
    """``gtf`` filtered by ``cellranger mkgtf`` (attributes in ``args``) into ``<prefix>.gtf``."""
    if gtf.name == f"{prefix}.gtf":
        raise ValueError("Input and output names are the same, use a different prefix to disambiguate!")
    results = await cellranger_mkgtf_cmd(gtf=gtf, prefix=prefix, args=args)
    return await pick_file(results, f"{prefix}.gtf", "cellranger mkgtf")


@dataclass
class CellrangerReference:
    reference: File  # <reference_name>.tar.gz, an archive of the <reference_name> folder
    reference_json: File  # <reference_name>/reference.json


async def cellranger_mkref(fasta: File, gtf: File, reference_name: str, args: str = "") -> CellrangerReference:
    """A Cell Ranger transcriptome reference from ``fasta`` and ``gtf``: ``<reference_name>.tar.gz``."""
    results = await cellranger_mkref_cmd(fasta=fasta, gtf=gtf, reference_name=reference_name, args=args)
    return CellrangerReference(
        reference=await pick_file(results, f"{reference_name}.tar.gz", "cellranger mkref"),
        reference_json=await pick_file(results, "reference.json", "cellranger mkref"),
    )


@dataclass
class CountResult:
    archive: File  # <prefix>.outs.tar.gz, an archive of <prefix>/outs
    metrics_summary: File  # outs/metrics_summary.csv
    web_summary: File  # outs/web_summary.html
    filtered_h5: File  # outs/filtered_feature_bc_matrix.h5
    raw_h5: File  # outs/raw_feature_bc_matrix.h5
    molecule_info: File  # outs/molecule_info.h5
    bam: File | None  # outs/possorted_genome_bam.bam, with --create-bam true
    versions: File  # versions.yml
    results: Dir  # everything above


async def cellranger_count(
    reads: list[File], reference: File, prefix: str, args: str = "", sample_id: str | None = None
) -> CountResult:
    """Cell Ranger gene-expression counts of one sample's ``reads`` (R1/R2 pairs) against ``reference`` (``.tar.gz``).

    The FASTQs are renamed ``<sample_id>_S1_L00<n>_R{1,2}_001.fastq.gz`` (``sample_id`` defaults to ``prefix``).
    """
    if not reads or len(reads) % 2:
        raise ValueError(f"cellranger_count needs R1/R2 pairs of FASTQs, got {len(reads)} files")
    results = await cellranger_count_cmd(
        reads=reads, reference=reference, sample_id=sample_id or prefix, prefix=prefix, args=args
    )
    return await read_count_results(results, prefix)


async def read_count_results(results: Dir, prefix: str) -> CountResult:
    """Pick the outputs out of ``cellranger_count_cmd``'s results Dir (e.g. from an overridden task)."""
    tool = "cellranger count"
    return CountResult(
        archive=await pick_file(results, f"{prefix}.outs.tar.gz", tool),
        metrics_summary=await pick_file(results, "metrics_summary.csv", tool),
        web_summary=await pick_file(results, "web_summary.html", tool),
        filtered_h5=await pick_file(results, "filtered_feature_bc_matrix.h5", tool),
        raw_h5=await pick_file(results, "raw_feature_bc_matrix.h5", tool),
        molecule_info=await pick_file(results, "molecule_info.h5", tool),
        bam=await results.get_file("possorted_genome_bam.bam"),
        versions=await pick_file(results, "versions.yml", tool),
        results=results,
    )


@dataclass
class MultiLibrary:
    """One library of a ``cellranger multi`` sample: its FASTQs (R1/R2 pairs) and options.

    ``options`` are the config keys upstream forwards (for gene expression:
    ``expect-cells``, ``chemistry``, ``create-bam``, ``no-secondary``, ...;
    ``r1-length``/``r2-length`` for the others; ``min-assignment-confidence``
    for multiplexing). ``fastq_id`` is the library's ``fastq_id`` in the config
    and defaults to the sample's prefix, which the renamed FASTQs carry.
    """

    reads: list[File]
    options: dict[str, str | int | bool] = field(default_factory=dict)
    fastq_id: str | None = None


def option_text(value: str | int | bool) -> str:
    """A config value as Groovy prints it: booleans in lower case."""
    if isinstance(value, bool):
        return "true" if value else "false"
    return str(value)


def multi_reference_paths(
    gex_reference: str | None, vdj_reference: str | None, files: dict[str, str]
) -> dict[str, str]:
    """Where :func:`cellranger_multi` stages its references, by kind, from their file names.

    ``gex_reference`` and ``vdj_reference`` are the archives' names (their
    folders are named after them, without ``.tar.gz``); ``files`` maps the
    kinds of :data:`MULTI_REFERENCE_DIRS` to their files' names.
    """
    paths = {kind: f"{MULTI_REFERENCE_DIRS[kind]}/{name}" for kind, name in files.items()}
    if gex_reference is not None:
        paths["gex"] = "references/gex/" + gex_reference.removesuffix(".tar.gz")
    if vdj_reference is not None:
        paths["vdj"] = "references/vdj/" + vdj_reference.removesuffix(".tar.gz")
    return paths


def cellranger_multi_config(
    prefix: str,
    libraries: dict[str, MultiLibrary],
    references: dict[str, str],
) -> str:
    """The ``cellranger multi`` config CSV upstream writes, with ``$PWD``-relative paths left for the shell.

    ``libraries`` maps a modality (``gex``, ``vdj``, ``vdj_t``, ``vdj_t_gd``,
    ``vdj_b``, ``ab``, ``beam``, ``cmo``, ``crispr``) to its library;
    ``references`` maps a reference kind (``gex``, ``vdj`` and the single-file
    kinds of :data:`MULTI_REFERENCE_DIRS`) to its path under the working
    directory (:func:`multi_reference_paths`). ``[samples]`` and ``[antigen-specificity]`` sections are
    ``$(cat ...)`` of their files, as upstream.
    """
    ref = references.get

    def options(modality: str) -> dict[str, str | int | bool]:
        library = libraries.get(modality)
        return library.options if library else {}

    def fastq_id(modality: str) -> str:
        return libraries[modality].fastq_id or prefix

    # Validate mutually exclusive barcode types
    if sum(1 for kind in ("ocm_barcodes", "cmo_barcodes", "frna_sampleinfo") if ref(kind)) >= 2:
        raise ValueError(
            "The ocm barcodes, cmo barcodes, and frna probes are mutually exclusive features. "
            "Please use only one per sample."
        )

    # Determine which library types are present based on FASTQs and references
    has_gex = "gex" in libraries and bool(ref("gex"))
    has_vdj = "vdj" in libraries and bool(ref("vdj"))
    has_ab = "ab" in libraries and bool(ref("fb_reference"))
    has_beam = "beam" in libraries and bool(ref("beam_control_panel"))
    has_cmo = "cmo" in libraries and bool(ref("cmo_barcodes"))
    has_crispr = "crispr" in libraries and bool(ref("fb_reference"))
    has_frna = bool(ref("gex_frna_probeset")) and bool(ref("frna_sampleinfo"))
    has_ocm = bool(ref("ocm_barcodes"))

    has_vdj_t = "vdj_t" in libraries and bool(ref("vdj"))
    has_vdj_t_gd = "vdj_t_gd" in libraries and bool(ref("vdj"))
    has_vdj_b = "vdj_b" in libraries and bool(ref("vdj"))
    has_any_vdj = has_vdj or has_vdj_t or has_vdj_t_gd or has_vdj_b

    if has_vdj_t_gd and not ref("vdj_primer_index"):
        raise ValueError(
            "VDJ-T-GD libraries require inner enrichment primers. Please provide them via the vdj_primer_index input."
        )

    # [gene-expression]: Cell Ranger requires the transcriptome reference for feature barcode libraries,
    # even without GEX FASTQs
    gex_section = []
    if has_gex or (ref("gex") and (has_ab or has_crispr or has_beam)):
        gex_options = options("gex")
        gex_section.append("[gene-expression]")
        gex_section.append(f"reference,$PWD/{ref('gex')}")
        # Flex probe sets are mutually exclusive with CMO or OCM barcodes and V(D)J libraries
        if ref("gex_frna_probeset") and not has_cmo and not has_ocm and not has_any_vdj:
            gex_section.append(f"probe-set,$PWD/{ref('gex_frna_probeset')}")
        for key in (
            "filter-probes",
            "r1-length",
            "r2-length",
            "chemistry",
            "expect-cells",
            "force-cells",
            "no-secondary",
            "check-library-compatibility",
            "no-target-umi-filter",
            "include-introns",
        ):
            if key in gex_options:
                gex_section.append(f"{key},{option_text(gex_options[key])}")
        # create-bam defaults to true if not specified
        create_bam = gex_options.get("create-bam")
        gex_section.append(f"create-bam,{option_text(True if create_bam is None else create_bam)}")
        if ref("gex_targetpanel"):
            gex_section.append(f"target-panel,$PWD/{ref('gex_targetpanel')}")
        # CMO-related settings that live inside the [gene-expression] section
        if has_cmo:
            cmo_options = options("cmo")
            if "min-assignment-confidence" in cmo_options:
                confidence = option_text(cmo_options["min-assignment-confidence"])
                gex_section.append(f"min-assignment-confidence,{confidence}")
            if ref("cmo_reference"):
                gex_section.append(f"cmo-set,$PWD/{ref('cmo_reference')}")
            if ref("cmo_barcode_assignment"):
                gex_section.append(f"barcode-sample-assignment,$PWD/{ref('cmo_barcode_assignment')}")

    # [feature]
    fb_section = []
    if has_ab or has_crispr or has_beam:
        fb_section.append("[feature]")
        if has_ab or has_crispr:
            fb_section.append("reference,$PWD/fb_reference_copy.csv")
        if has_beam and ref("beam_antigen_panel"):
            fb_section.append(f"reference,$PWD/{ref('beam_antigen_panel')}")
        # r1/r2-length from the antibody options takes priority over the CRISPR ones
        fb_options = options("ab") if has_ab else (options("crispr") if has_crispr else {})
        for key in ("r1-length", "r2-length"):
            if key in fb_options:
                fb_section.append(f"{key},{option_text(fb_options[key])}")

    # [vdj]: all V(D)J libraries share one section; the first present library's options are used
    vdj_section = []
    if has_any_vdj:
        vdj_section.append("[vdj]")
        vdj_section.append(f"reference,$PWD/{ref('vdj')}")
        if ref("vdj_primer_index"):
            vdj_section.append(f"inner-enrichment-primers,$PWD/{ref('vdj_primer_index')}")
        vdj_options = (
            options("vdj")
            if has_vdj
            else options("vdj_t")
            if has_vdj_t
            else options("vdj_t_gd")
            if has_vdj_t_gd
            else options("vdj_b")
        )
        for key in ("r1-length", "r2-length"):
            if key in vdj_options:
                vdj_section.append(f"{key},{option_text(vdj_options[key])}")

    # [libraries]
    lib_section = ["[libraries]", "fastq_id,fastqs,lanes,feature_types"]
    for present, modality, feature_type in (
        (has_gex, "gex", "Gene Expression"),
        (has_vdj, "vdj", "VDJ"),
        (has_vdj_t, "vdj_t", "VDJ-T"),
        (has_vdj_t_gd, "vdj_t_gd", "VDJ-T-GD"),
        (has_vdj_b, "vdj_b", "VDJ-B"),
        (has_ab, "ab", "Antibody Capture"),
        (has_beam, "beam", "Antigen Capture"),
        (has_crispr, "crispr", "CRISPR Guide Capture"),
        (has_cmo, "cmo", "Multiplexing Capture"),
    ):
        if present:
            lib_section.append(f"{fastq_id(modality)},$PWD/fastq_all/{modality},,{feature_type}")

    lines = gex_section + fb_section + vdj_section + lib_section
    # Sample sections
    if has_cmo:
        lines += ["[samples]", f"$(cat {ref('cmo_barcodes')})"]
    if has_frna:
        lines += ["[samples]", f"$(cat {ref('frna_sampleinfo')})"]
    if has_ocm:
        lines += ["[samples]", f"$(cat {ref('ocm_barcodes')})"]
    if has_beam:
        lines += ["[antigen-specificity]", f"$(cat {ref('beam_control_panel')})"]
    return "\n".join(line for line in lines if line)


@dataclass
class MultiResult:
    config: File  # cellranger_multi_config.csv
    archive: File  # <prefix>.outs.tar.gz, an archive of <prefix>/outs
    results: Dir  # both


async def cellranger_multi(
    prefix: str,
    gex: MultiLibrary | None = None,
    vdj: MultiLibrary | None = None,
    ab: MultiLibrary | None = None,
    beam: MultiLibrary | None = None,
    cmo: MultiLibrary | None = None,
    crispr: MultiLibrary | None = None,
    vdj_t: MultiLibrary | None = None,
    vdj_t_gd: MultiLibrary | None = None,
    vdj_b: MultiLibrary | None = None,
    gex_reference: File | None = None,
    vdj_reference: File | None = None,
    gex_frna_probeset: File | None = None,
    gex_targetpanel: File | None = None,
    vdj_primer_index: File | None = None,
    fb_reference: File | None = None,
    beam_antigen_panel: File | None = None,
    beam_control_panel: File | None = None,
    cmo_reference: File | None = None,
    cmo_barcodes: File | None = None,
    cmo_barcode_assignment: File | None = None,
    frna_sampleinfo: File | None = None,
    ocm_barcodes: File | None = None,
    skip_renaming: bool = False,
    args: str = "",
) -> MultiResult:
    """``cellranger multi`` on one sample's libraries.

    ``gex_reference`` (from :func:`cellranger_mkref`) and ``vdj_reference``
    (from :func:`cellranger_mkvdjref`) are ``.tar.gz`` archives; the other
    references are files. A library's section is written only when its
    reference is given (as upstream).
    """
    given = {
        "gex": gex,
        "vdj": vdj,
        "vdj_t": vdj_t,
        "vdj_t_gd": vdj_t_gd,
        "vdj_b": vdj_b,
        "ab": ab,
        "beam": beam,
        "cmo": cmo,
        "crispr": crispr,
    }
    libraries = {modality: library for modality, library in given.items() if library is not None}
    for modality, library in libraries.items():
        if not library.reads or len(library.reads) % 2:
            raise ValueError(f"cellranger_multi's {modality} library needs R1/R2 pairs, got {len(library.reads)} files")
    files = {
        "gex_frna_probeset": gex_frna_probeset,
        "gex_targetpanel": gex_targetpanel,
        "vdj_primer_index": vdj_primer_index,
        "fb_reference": fb_reference,
        "beam_antigen_panel": beam_antigen_panel,
        "beam_control_panel": beam_control_panel,
        "cmo_reference": cmo_reference,
        "cmo_barcodes": cmo_barcodes,
        "cmo_barcode_assignment": cmo_barcode_assignment,
        "frna_sampleinfo": frna_sampleinfo,
        "ocm_barcodes": ocm_barcodes,
    }
    references = multi_reference_paths(
        gex_reference.name if gex_reference else None,
        vdj_reference.name if vdj_reference else None,
        {kind: file.name for kind, file in files.items() if file is not None},
    )
    config = cellranger_multi_config(prefix, libraries, references)
    results = await cellranger_multi_cmd(
        **{f"{modality}_reads": library.reads if library else None for modality, library in given.items()},
        gex_reference=gex_reference,
        vdj_reference=vdj_reference,
        **files,
        config=config,
        prefix=prefix,
        skip_renaming="true" if skip_renaming else "false",
        args=args,
    )
    tool = "cellranger multi"
    return MultiResult(
        config=await pick_file(results, "cellranger_multi_config.csv", tool),
        archive=await pick_file(results, f"{prefix}.outs.tar.gz", tool),
        results=results,
    )


@dataclass
class MkfastqResult:
    fastq: File  # fastq.tar.gz: every FASTQ under outs/fastq_path, with its path there
    reports: File | None  # Reports.tar.gz
    stats: File | None  # Stats.tar.gz
    interop: Dir  # outs/interop_path/*.bin
    results: Dir  # everything above


async def cellranger_mkfastq(samplesheet: File, run_dir: File, prefix: str, args: str = "") -> MkfastqResult:
    """Demultiplex the run folder archive ``run_dir`` (``.tar.gz``) with ``samplesheet`` (``cellranger mkfastq``)."""
    results = await cellranger_mkfastq_cmd(samplesheet=samplesheet, run_dir=run_dir, prefix=prefix, args=args)
    return MkfastqResult(
        fastq=await pick_file(results, "fastq.tar.gz", "cellranger mkfastq"),
        reports=await results.get_file("Reports.tar.gz"),
        stats=await results.get_file("Stats.tar.gz"),
        interop=Dir.from_existing_remote(results.path.rstrip("/") + "/interop"),
        results=results,
    )


@dataclass
class VdjReference:
    reference: File  # <reference_name>.tar.gz, an archive of the <reference_name> folder
    regions: File  # <reference_name>/fasta/regions.fa
    reference_json: File  # <reference_name>/reference.json


async def cellranger_mkvdjref(
    reference_name: str,
    fasta: File | None = None,
    gtf: File | None = None,
    seqs: File | None = None,
    args: str = "",
) -> VdjReference:
    """A Cell Ranger V(D)J reference from ``fasta`` and ``gtf``, or from ``seqs``: ``<reference_name>.tar.gz``."""
    if seqs is None and (fasta is None or gtf is None):
        raise ValueError("cellranger_mkvdjref needs a fasta and gtf, or seqs")
    results = await cellranger_mkvdjref_cmd(fasta=fasta, gtf=gtf, seqs=seqs, reference_name=reference_name, args=args)
    files = {}
    for name in (f"{reference_name}.tar.gz", "regions.fa", "reference.json"):
        files[name] = await results.get_file(name)
        if files[name] is None:
            raise FileNotFoundError(f"cellranger mkvdjref wrote no {name}")
    return VdjReference(
        reference=files[f"{reference_name}.tar.gz"],
        regions=files["regions.fa"],
        reference_json=files["reference.json"],
    )


@dataclass
class VdjResult:
    outs: Dir  # <prefix>/outs
    archive: File  # <prefix>.outs.tar.gz, an archive of <prefix>/outs
    airr_rearrangement: File | None  # outs/airr_rearrangement.tsv
    filtered_contig_annotations: File | None  # outs/filtered_contig_annotations.csv
    metrics_summary: File | None  # outs/metrics_summary.csv
    web_summary: File | None  # outs/web_summary.html


async def cellranger_vdj(reads: list[File], reference: File, prefix: str, args: str = "") -> VdjResult:
    """Cell Ranger V(D)J assembly and clonotypes of one sample's ``reads`` against ``reference`` (a ``.tar.gz``)."""
    if not reads:
        raise ValueError("cellranger_vdj needs the sample's FASTQs")
    results = await cellranger_vdj_cmd(reads=reads, reference=reference, prefix=prefix, args=args)
    archive = await results.get_file(f"{prefix}.outs.tar.gz")
    if archive is None:
        raise FileNotFoundError(f"cellranger vdj wrote no {prefix}.outs.tar.gz")
    outs = Dir.from_existing_remote(results.path.rstrip("/") + "/outs")
    return VdjResult(
        outs=outs,
        archive=archive,
        airr_rearrangement=await outs.get_file("airr_rearrangement.tsv"),
        filtered_contig_annotations=await outs.get_file("filtered_contig_annotations.csv"),
        metrics_summary=await outs.get_file("metrics_summary.csv"),
        web_summary=await outs.get_file("web_summary.html"),
    )
