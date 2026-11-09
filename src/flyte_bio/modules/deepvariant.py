"""deepvariant — germline small-variant calling with DeepVariant.

Exposes :func:`deepvariant_rundeepvariant`, which runs DeepVariant's whole
pipeline (``run_deepvariant``: make_examples, call_variants,
postprocess_variants) on one sample's BAM or CRAM with its index, and writes,
under upstream's names:

- ``<prefix>.vcf.gz`` with its index ``<prefix>.vcf.gz.tbi`` (or ``.csi``);
- ``<prefix>.g.vcf.gz`` with its index ``<prefix>.g.vcf.gz.tbi`` (or ``.csi``);
- ``<prefix>.visual_report.html``, when DeepVariant writes one.

``fasta`` is the reference, plain or bgzipped; its ``fai`` (and, bgzipped,
its ``gzi``) are linked next to it. ``intervals`` (a BED) restricts calling
(``--regions``) and ``par_bed`` gives the pseudo-autosomal regions of the
haploid contigs (``--par_regions_bed``). ``args`` should set
``--model_type`` (``WGS``, ``WES``, ``PACBIO``, ...), which DeepVariant
requires. The number of shards is the task's CPUs.

Its three steps also run as separate tasks, as upstream's split modules do:

- :func:`deepvariant_makeexamples` — make the examples for calling (one
  ``make_examples`` per CPU, under GNU parallel): ``<prefix>.examples.tfrecord-
  NNNNN-of-MMMMM.gz`` (each with its ``.example_info.json``) and
  ``<prefix>.gvcf.tfrecord-NNNNN-of-MMMMM.gz``, and, with the small model,
  ``<prefix>_call_variant_outputs.examples.tfrecord-NNNNN-of-MMMMM.gz`` (a
  prefix ``a.b`` gives ``a_call_variant_outputs.b.examples...``). With
  ``with_phasing`` it phases reads and writes
  ``<prefix>-read-phasing_debug-NNNNN-of-MMMMM.tsv``. The shard count is the
  task's CPUs, so it shows in every name.
- :func:`deepvariant_callvariants` — run the model over the examples:
  ``<prefix>.call-NNNNN-of-MMMMM.tfrecord.gz``. On CPU, as here; upstream
  can also run it on a GPU.
- :func:`deepvariant_postprocessvariants` — turn the calls (and the gVCF
  records, small-model calls and read phasing, when given) into
  ``<prefix>.vcf.gz`` and ``<prefix>.g.vcf.gz`` with their indexes.

The sharded inputs are passed as lists of files; the tools are given their
sharded names (``name@N.gz``), worked out from the first file's name as
upstream does. :func:`deepvariant_vcfstatsreport` draws DeepVariant's visual
report of any VCF or BCF: ``<prefix>.visual_report.html``.

For each step ``args`` goes to the tool: the model's ``--checkpoint`` (e.g.
``/opt/models/wgs``) for make_examples and call_variants, and for the small
model ``--call_small_model_examples --trained_small_model_path ...``. Values
in ``args`` aren't shell-quoted.
"""

import re
from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

# Upstream's image is docker.io/google/deepvariant:1.10.0; this is Google's own
# copy of it on gcr.io (the same image: its config digest is identical), which
# the cluster can pull.
DEEPVARIANT_IMAGE = "gcr.io/deepvariant-docker/deepvariant:1.10.0"
# Upstream's label is process_high (12 CPUs, 72 GB). This is enough for the test
# data; a whole genome needs far more CPUs (and memory to match) to finish in
# reasonable time: override with deepvariant_rundeepvariant_cmd.as_task().override(resources=...).
DEFAULT_RESOURCES = flyte.Resources(cpu=4, memory="16Gi", disk="20Gi")

deepvariant_rundeepvariant_cmd = shell.create(
    name="deepvariant_rundeepvariant",
    image=DEEPVARIANT_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "reads": File,
        "index": File,
        "intervals": File | None,
        "fasta": File,
        "fai": File,
        "gzi": File | None,
        "par_bed": File | None,
        "prefix": str,
        "args": str,
    },
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        READS=({{inputs.reads}}); INDEX=({{inputs.index}})
        FA=({{inputs.fasta}}); FAI=({{inputs.fai}})
        shopt -s nullglob
        INTERVALS=({{inputs.intervals}}); GZI=({{inputs.gzi}}); PAR=({{inputs.par_bed}})
        shopt -u nullglob
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"

        # Reads and reference beside their indexes, under their own names.
        READS_NAME=$(basename "${{READS[0]}}"); ln -s "${{READS[0]}}" "$READS_NAME"
        INDEX_NAME=$(basename "${{INDEX[0]}}"); ln -s "${{INDEX[0]}}" "$INDEX_NAME"
        FA_NAME=$(basename "${{FA[0]}}"); ln -s "${{FA[0]}}" "$FA_NAME"
        ln -s "${{FAI[0]}}" "$FA_NAME.fai"
        if [ ${{#GZI[@]}} -gt 0 ]; then ln -s "${{GZI[0]}}" "$FA_NAME.gzi"; fi
        EXTRA=()
        if [ ${{#INTERVALS[@]}} -gt 0 ]; then EXTRA+=(--regions="${{INTERVALS[0]}}"); fi
        if [ ${{#PAR[@]}} -gt 0 ]; then EXTRA+=(--par_regions_bed="${{PAR[0]}}"); fi

        export MPLCONFIGDIR=$W/.matplotlib
        mkdir -p "$MPLCONFIGDIR"

        /opt/deepvariant/bin/run_deepvariant \
            --ref="$FA_NAME" \
            --reads="$READS_NAME" \
            --output_vcf="$P.vcf.gz" \
            --output_gvcf="$P.g.vcf.gz" \
            $ARGS \
            "${{EXTRA[@]}}" \
            --intermediate_results_dir=tmp \
            --num_shards=$TASK_CPUS

        shopt -s nullglob
        for f in "$P".vcf.gz* "$P".g.vcf.gz* "$P".visual_report.htm[l]; do mv "$f" "$OUT"/; done
    """,
)

# Upstream labels make_examples and call_variants process_high, postprocess_variants process_medium and the
# report process_single; its tests run them with 4 CPUs.
REPORT_RESOURCES = flyte.Resources(cpu=1, memory="4Gi")

# Links the reads, the reference and their indexes into the working directory under their own names, and sets
# READS_NAME and FA_NAME. Expects READS, INDEX, FA, FAI and GZI arrays.
LINK_READS_AND_REFERENCE = r"""
        READS_NAME=$(basename "${READS[0]}"); ln -s "${READS[0]}" "$READS_NAME"
        ln -s "${INDEX[0]}" "$(basename "${INDEX[0]}")"
        FA_NAME=$(basename "${FA[0]}"); ln -s "${FA[0]}" "$FA_NAME"
        ln -s "${FAI[0]}" "$FA_NAME.fai"
        if [ ${#GZI[@]} -gt 0 ]; then ln -s "${GZI[0]}" "$FA_NAME.gzi"; fi
"""

# Links a list of sharded files into the working directory and prints their sharded name: the first file's
# name with "-NNNNN-of-MMMMM" before its extension replaced by "@MMMMM", as upstream works it out.
SHARDED_NAME = r"""
        sharded_name() {
            local first="" f name
            for f in "$@"; do
                name=$(basename "$f")
                ln -sf "$f" "$name"
                case "$name" in *.example_info.json) continue ;; esac
                if [ -z "$first" ]; then first=$name; fi
            done
            local base=${first%.*} extension=${first##*.}
            if ! [[ $base =~ ^(.+)-[0-9]{5}-of-([0-9]{5})$ ]]; then
                echo "sharded file name '$first' doesn't match <name>-NNNNN-of-MMMMM" >&2
                exit 1
            fi
            echo "${BASH_REMATCH[1]}@${BASH_REMATCH[2]}.$extension"
        }
"""

deepvariant_makeexamples_cmd = shell.create(
    name="deepvariant_makeexamples",
    image=DEEPVARIANT_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "reads": File,
        "index": File,
        "intervals": File | None,
        "fasta": File,
        "fai": File,
        "gzi": File | None,
        "par_bed": File | None,
        "with_phasing": bool,
        "prefix": str,
        "args": str,
    },
    defaults={"with_phasing": False, "args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        READS=({{inputs.reads}}); INDEX=({{inputs.index}})
        FA=({{inputs.fasta}}); FAI=({{inputs.fai}})
        shopt -s nullglob
        INTERVALS=({{inputs.intervals}}); GZI=({{inputs.gzi}}); PAR=({{inputs.par_bed}})
        shopt -u nullglob
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        {LINK_READS_AND_REFERENCE}
        EXTRA=()
        if [ ${{#INTERVALS[@]}} -gt 0 ]; then EXTRA+=(--regions "${{INTERVALS[0]}}"); fi
        if [ ${{#PAR[@]}} -gt 0 ]; then EXTRA+=(--par_regions_bed="${{PAR[0]}}"); fi
        PHASING=()
        if [ {{inputs.with_phasing}} = true ]; then
            PHASING=(--phase_reads --track_ref_reads --output_phase_info)
            PHASING+=(--output_local_read_phasing "./$P-read-phasing_debug@$TASK_CPUS.tsv")
        fi

        export MPLCONFIGDIR=$W/.matplotlib
        mkdir -p "$MPLCONFIGDIR"

        seq 0 $(( TASK_CPUS - 1 )) | parallel -q --halt 2 --line-buffer /opt/deepvariant/bin/make_examples \
            --mode calling \
            --ref "$FA_NAME" \
            --reads "$READS_NAME" \
            --examples "./$P.examples.tfrecord@$TASK_CPUS.gz" \
            --gvcf "./$P.gvcf.tfrecord@$TASK_CPUS.gz" \
            "${{EXTRA[@]}}" \
            $ARGS \
            "${{PHASING[@]}}" \
            --task {{}}

        shopt -s nullglob
        mv ./*.tfrecord-*-of-*.gz* ./*-read-phasing_debug-*-of-*.tsv "$OUT"/
    """,
)

deepvariant_callvariants_cmd = shell.create(
    name="deepvariant_callvariants",
    image=DEEPVARIANT_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"examples": list[File], "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        EXAMPLES=({{inputs.examples}})
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        W=$(mktemp -d); cd "$W"
        {SHARDED_NAME}
        EXAMPLES_NAME=$(sharded_name "${{EXAMPLES[@]}}") || exit 1

        export MPLCONFIGDIR=$W/.matplotlib
        mkdir -p "$MPLCONFIGDIR"

        /opt/deepvariant/bin/call_variants \
            $ARGS \
            --outfile "$P.call.tfrecord.gz" \
            --examples "$EXAMPLES_NAME"

        mv ./"$P".call-*-of-*.tfrecord.gz {{outputs.results}}/
    """,
)

deepvariant_postprocessvariants_cmd = shell.create(
    name="deepvariant_postprocessvariants",
    image=DEEPVARIANT_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "calls": list[File],
        "gvcf_records": list[File],
        "small_model_calls": list[File],
        "phase_inputs": list[File],
        "intervals": File | None,
        "fasta": File,
        "fai": File,
        "gzi": File | None,
        "prefix": str,
        "args": str,
    },
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        CALLS=({{inputs.calls}}); GVCF=({{inputs.gvcf_records}})
        FA=({{inputs.fasta}}); FAI=({{inputs.fai}})
        shopt -s nullglob
        SMALL=({{inputs.small_model_calls}}); PHASE=({{inputs.phase_inputs}})
        INTERVALS=({{inputs.intervals}}); GZI=({{inputs.gzi}})
        shopt -u nullglob
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        FA_NAME=$(basename "${{FA[0]}}"); ln -s "${{FA[0]}}" "$FA_NAME"
        ln -s "${{FAI[0]}}" "$FA_NAME.fai"
        if [ ${{#GZI[@]}} -gt 0 ]; then ln -s "${{GZI[0]}}" "$FA_NAME.gzi"; fi
        {SHARDED_NAME}
        # As upstream, the calls are named by the first file's name without its shard numbers.
        for f in "${{CALLS[@]}}"; do ln -s "$f" "$(basename "$f")"; done
        CALLS_NAME=$(basename "${{CALLS[0]}}" | sed -E 's/-[0-9]{{5}}-of-[0-9]{{5}}//')
        GVCF_NAME=$(sharded_name "${{GVCF[@]}}") || exit 1
        EXTRA=()
        if [ ${{#INTERVALS[@]}} -gt 0 ]; then EXTRA+=(--regions "${{INTERVALS[0]}}"); fi
        if [ ${{#SMALL[@]}} -gt 0 ]; then
            SMALL_NAME=$(sharded_name "${{SMALL[@]}}") || exit 1
            EXTRA+=(--small_model_cvo_records "$SMALL_NAME")
        fi
        if [ ${{#PHASE[@]}} -gt 0 ]; then
            PHASE_NAME=$(sharded_name "${{PHASE[@]}}") || exit 1
            EXTRA+=(--phased_reads_input_path "$PHASE_NAME")
        fi

        export MPLCONFIGDIR=$W/.matplotlib
        mkdir -p "$MPLCONFIGDIR"

        /opt/deepvariant/bin/postprocess_variants \
            $ARGS \
            --ref "$FA_NAME" \
            --infile "$CALLS_NAME" \
            --outfile "$P.vcf.gz" \
            --nonvariant_site_tfrecord_path "$GVCF_NAME" \
            --gvcf_outfile "$P.g.vcf.gz" \
            "${{EXTRA[@]}}" \
            --cpus $TASK_CPUS

        shopt -s nullglob
        mv "$P".vcf.gz* "$P".g.vcf.gz* "$OUT"/
    """,
)

deepvariant_vcfstatsreport_cmd = shell.create(
    name="deepvariant_vcfstatsreport",
    image=DEEPVARIANT_IMAGE,
    resources=REPORT_RESOURCES,
    inputs={"vcf": File, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        VCF=({inputs.vcf})
        P={inputs.prefix}
        ARGS={inputs.args}
        W=$(mktemp -d); cd "$W"
        NAME=$(basename "${VCF[0]}"); ln -s "${VCF[0]}" "$NAME"
        export MPLCONFIGDIR=$W/.matplotlib
        mkdir -p "$MPLCONFIGDIR"
        /opt/deepvariant/bin/vcf_stats_report \
            --input_vcf="$NAME" \
            --outfile_base "$P" \
            $ARGS
        mv "$P".visual_report.html {outputs.results}/
    """,
)

env = flyte.TaskEnvironment.from_task(
    "deepvariant",
    deepvariant_rundeepvariant_cmd.as_task(),
    deepvariant_makeexamples_cmd.as_task(),
    deepvariant_callvariants_cmd.as_task(),
    deepvariant_postprocessvariants_cmd.as_task(),
    deepvariant_vcfstatsreport_cmd.as_task(),
)


@dataclass
class DeepVariantResult:
    vcf: File  # <prefix>.vcf.gz
    vcf_index: File  # its .tbi or .csi
    gvcf: File  # <prefix>.g.vcf.gz
    gvcf_index: File  # its .tbi or .csi
    report: File | None  # <prefix>.visual_report.html, when written
    results: Dir  # everything above


async def deepvariant_rundeepvariant(
    reads: File,
    index: File,
    fasta: File,
    fai: File,
    intervals: File | None = None,
    gzi: File | None = None,
    par_bed: File | None = None,
    prefix: str = "sample",
    args: str = "",
) -> DeepVariantResult:
    """Call variants in ``reads`` (a BAM or CRAM with its ``index``) against ``fasta``, as VCF and gVCF."""
    results = await deepvariant_rundeepvariant_cmd(
        reads=reads,
        index=index,
        intervals=intervals,
        fasta=fasta,
        fai=fai,
        gzi=gzi,
        par_bed=par_bed,
        prefix=prefix,
        args=args,
    )

    async def required(*names: str) -> File:
        for name in names:
            file = await results.get_file(name)
            if file is not None:
                return file
        raise FileNotFoundError(f"DeepVariant wrote none of {', '.join(names)}")

    return DeepVariantResult(
        vcf=await required(f"{prefix}.vcf.gz"),
        vcf_index=await required(f"{prefix}.vcf.gz.tbi", f"{prefix}.vcf.gz.csi"),
        gvcf=await required(f"{prefix}.g.vcf.gz"),
        gvcf_index=await required(f"{prefix}.g.vcf.gz.tbi", f"{prefix}.g.vcf.gz.csi"),
        report=await results.get_file(f"{prefix}.visual_report.html"),
        results=results,
    )


SHARD = re.compile(r"-\d{5}-of-\d{5}")


@dataclass
class DeepVariantExamples:
    examples: list[File]  # <prefix>.examples.tfrecord-NNNNN-of-MMMMM.gz and their .example_info.json
    gvcf: list[File]  # <prefix>.gvcf.tfrecord-NNNNN-of-MMMMM.gz
    small_model_calls: list[File]  # <prefix>_call_variant_outputs.examples.tfrecord-*.gz, with the small model
    read_phase_inputs: list[File]  # <prefix>-read-phasing_debug-NNNNN-of-MMMMM.tsv, with phasing
    results: Dir  # everything above


async def deepvariant_examples(results: Dir, prefix: str) -> DeepVariantExamples:
    """Sort a make_examples run's files (from ``deepvariant_makeexamples_cmd``) into :class:`DeepVariantExamples`."""
    first, _, rest = prefix.partition(".")
    small_model = f"{first}_call_variant_outputs.{rest + '.' if rest else ''}examples.tfrecord-"
    files = sorted([f async for f in results.walk()], key=lambda f: f.name)

    def named(start: str, end: str = "") -> list[File]:
        return [f for f in files if f.name.startswith(start) and f.name.endswith(end) and SHARD.search(f.name)]

    examples = named(f"{prefix}.examples.tfrecord-", ".gz") + named(f"{prefix}.examples.tfrecord-", ".json")
    if not examples:
        raise FileNotFoundError(f"DeepVariant make_examples wrote no {prefix}.examples.tfrecord-*.gz")
    return DeepVariantExamples(
        examples=sorted(examples, key=lambda f: f.name),
        gvcf=named(f"{prefix}.gvcf.tfrecord-", ".gz"),
        small_model_calls=named(small_model, ".gz"),
        read_phase_inputs=named(f"{prefix}-read-phasing_debug-", ".tsv"),
        results=results,
    )


async def deepvariant_makeexamples(
    reads: File,
    index: File,
    fasta: File,
    fai: File,
    intervals: File | None = None,
    gzi: File | None = None,
    par_bed: File | None = None,
    with_phasing: bool = False,
    prefix: str = "sample",
    args: str = "",
) -> DeepVariantExamples:
    """Make DeepVariant's examples from ``reads`` (a BAM or CRAM with its ``index``), one shard per CPU."""
    results = await deepvariant_makeexamples_cmd(
        reads=reads,
        index=index,
        intervals=intervals,
        fasta=fasta,
        fai=fai,
        gzi=gzi,
        par_bed=par_bed,
        with_phasing=with_phasing,
        prefix=prefix,
        args=args,
    )
    return await deepvariant_examples(results, prefix)


async def deepvariant_callvariants(examples: list[File], prefix: str = "sample", args: str = "") -> list[File]:
    """Call variants in make_examples' ``examples``: ``<prefix>.call-NNNNN-of-MMMMM.tfrecord.gz``."""
    results = await deepvariant_callvariants_cmd(examples=examples, prefix=prefix, args=args)
    calls = sorted([f async for f in results.walk()], key=lambda f: f.name)
    if not calls:
        raise FileNotFoundError(f"DeepVariant call_variants wrote no {prefix}.call-*.tfrecord.gz")
    return calls


async def deepvariant_postprocessvariants(
    calls: list[File],
    gvcf_records: list[File],
    fasta: File,
    fai: File,
    small_model_calls: list[File] | None = None,
    phase_inputs: list[File] | None = None,
    intervals: File | None = None,
    gzi: File | None = None,
    prefix: str = "sample",
    args: str = "",
) -> DeepVariantResult:
    """Write the VCF and gVCF from call_variants' ``calls`` and make_examples' ``gvcf_records``."""
    results = await deepvariant_postprocessvariants_cmd(
        calls=calls,
        gvcf_records=gvcf_records,
        small_model_calls=small_model_calls or [],
        phase_inputs=phase_inputs or [],
        intervals=intervals,
        fasta=fasta,
        fai=fai,
        gzi=gzi,
        prefix=prefix,
        args=args,
    )

    async def required(*names: str) -> File:
        for name in names:
            file = await results.get_file(name)
            if file is not None:
                return file
        raise FileNotFoundError(f"DeepVariant postprocess_variants wrote none of {', '.join(names)}")

    return DeepVariantResult(
        vcf=await required(f"{prefix}.vcf.gz"),
        vcf_index=await required(f"{prefix}.vcf.gz.tbi", f"{prefix}.vcf.gz.csi"),
        gvcf=await required(f"{prefix}.g.vcf.gz"),
        gvcf_index=await required(f"{prefix}.g.vcf.gz.tbi", f"{prefix}.g.vcf.gz.csi"),
        report=None,
        results=results,
    )


async def deepvariant_vcfstatsreport(vcf: File, prefix: str = "sample", args: str = "") -> File:
    """DeepVariant's visual report of ``vcf`` (VCF or BCF, plain or compressed): ``<prefix>.visual_report.html``."""
    results = await deepvariant_vcfstatsreport_cmd(vcf=vcf, prefix=prefix, args=args)
    report = await results.get_file(f"{prefix}.visual_report.html")
    if report is None:
        raise FileNotFoundError(f"DeepVariant vcf_stats_report wrote no {prefix}.visual_report.html")
    return report
