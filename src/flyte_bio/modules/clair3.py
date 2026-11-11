"""clair3 — Clair3: germline small-variant calling from long (and short) reads with a deep-learning model.

Exposes :func:`clair3`, which runs Clair3's whole pipeline
(``run_clair3.sh``: pileup calling, full-alignment calling, merging) on one
sample's BAM with its index against an indexed reference, and writes, under
upstream's names (``<prefix>`` is prepended with no separator, as upstream):

- ``<prefix>merge_output.vcf.gz`` with its ``.tbi``;
- with ``--gvcf`` in ``args``: ``<prefix>merge_output.gvcf.gz`` with its ``.tbi``;
- with phasing (``--enable_phasing``): ``<prefix>phased_merge_output.vcf.gz``
  with its ``.tbi``.

The model is either one packaged in the image (``packaged_model``, a
directory name under ``/opt/models``, e.g. ``hifi_revio`` or
``r1041_e82_400bps_sup_v500``) or ``user_model``, a directory holding the
model's files; give exactly one. ``platform`` is ``ont``, ``hifi`` or
``ilmn``. Clair3 gets the task's CPUs as threads. ``args`` go to
``run_clair3.sh``.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

# Upstream's image, from Docker Hub (Clair3's own build, with its models under
# /opt/models; bioconda's build packages none).
CLAIR3_IMAGE = "docker.io/hkubal/clair3:v2.0.0"
# Upstream labels it process_high (12 CPUs, 72 GB); its tests ran on 4 CPUs.
# A whole genome needs more: override with clair3_cmd.override(resources=...).
DEFAULT_RESOURCES = flyte.Resources(cpu=4, memory="16Gi", disk="20Gi")

OUTPUT_NAMES = (
    "merge_output.vcf.gz",
    "merge_output.vcf.gz.tbi",
    "phased_merge_output.vcf.gz",
    "phased_merge_output.vcf.gz.tbi",
    "merge_output.gvcf.gz",
    "merge_output.gvcf.gz.tbi",
)

clair3_cmd = shell.create(
    name="clair3",
    image=CLAIR3_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "bam": File,
        "bai": File,
        "packaged_model": str,
        "user_model": Dir | None,
        "platform": str,
        "reference": File,
        "fai": File,
        "prefix": str,
        "args": str,
    },
    defaults={"packaged_model": "", "args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        BAM=({{inputs.bam}}); BAI=({{inputs.bai}})
        FA=({{inputs.reference}}); FAI=({{inputs.fai}})
        PACKAGED_MODEL={{inputs.packaged_model}}
        USER_MODEL={{inputs.user_model}}
        PLATFORM={{inputs.platform}}
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"

        # Reads and reference beside their indexes, under their own names.
        BAM_NAME=$(basename "${{BAM[0]}}"); ln -s "${{BAM[0]}}" "$BAM_NAME"
        case "${{BAI[0]}}" in
            *.csi) ln -s "${{BAI[0]}}" "$BAM_NAME.csi" ;;
            *) ln -s "${{BAI[0]}}" "$BAM_NAME.bai" ;;
        esac
        FA_NAME=$(basename "${{FA[0]}}"); ln -s "${{FA[0]}}" "$FA_NAME"
        ln -s "${{FAI[0]}}" "$FA_NAME.fai"
        if [ -n "$USER_MODEL" ] && [ -d "$USER_MODEL" ]; then
            MODEL=$(cd "$USER_MODEL" && pwd)
        else
            MODEL=/opt/models/$PACKAGED_MODEL
        fi

        run_clair3.sh \
            --bam_fn="$W/$BAM_NAME" \
            --ref_fn="$W/$FA_NAME" \
            --threads=$TASK_CPUS \
            --output="$W" \
            --platform="$PLATFORM" \
            --model_path="$MODEL" \
            $ARGS || exit 1

        for f in {" ".join(OUTPUT_NAMES)}; do
            if [ -e "$f" ]; then mv "$f" "$OUT/$P$f"; fi
        done
    """,
)

env = flyte.TaskEnvironment.from_task("clair3", clair3_cmd.as_task())


@dataclass
class Clair3Result:
    vcf: File  # <prefix>merge_output.vcf.gz
    tbi: File  # its .tbi
    phased_vcf: File | None  # <prefix>phased_merge_output.vcf.gz, with phasing
    phased_tbi: File | None  # its .tbi
    gvcf: File | None  # <prefix>merge_output.gvcf.gz, with --gvcf
    gvcf_tbi: File | None  # its .tbi
    results: Dir  # everything above


async def clair3(
    bam: File,
    bai: File,
    reference: File,
    fai: File,
    platform: str,
    packaged_model: str = "",
    user_model: Dir | None = None,
    prefix: str = "sample",
    args: str = "",
) -> Clair3Result:
    """Call variants in ``bam`` (with its index ``bai``) against ``reference`` (with its ``fai``).

    Give exactly one model: ``packaged_model`` (a name under the image's
    ``/opt/models``) or ``user_model`` (a directory of model files).
    """
    if packaged_model and user_model is not None:
        raise ValueError(f"Two models specified {user_model.path} and {packaged_model}, specify one of them.")
    if not packaged_model and user_model is None:
        raise ValueError("clair3 needs a model: packaged_model or user_model")
    results = await clair3_cmd(
        bam=bam,
        bai=bai,
        packaged_model=packaged_model,
        user_model=user_model,
        platform=platform,
        reference=reference,
        fai=fai,
        prefix=prefix,
        args=args,
    )

    async def optional(name: str) -> File | None:
        return await results.get_file(prefix + name)

    vcf, tbi = await optional("merge_output.vcf.gz"), await optional("merge_output.vcf.gz.tbi")
    if vcf is None or tbi is None:
        raise FileNotFoundError(f"Clair3 wrote no {prefix}merge_output.vcf.gz or its index")
    return Clair3Result(
        vcf=vcf,
        tbi=tbi,
        phased_vcf=await optional("phased_merge_output.vcf.gz"),
        phased_tbi=await optional("phased_merge_output.vcf.gz.tbi"),
        gvcf=await optional("merge_output.gvcf.gz"),
        gvcf_tbi=await optional("merge_output.gvcf.gz.tbi"),
        results=results,
    )
