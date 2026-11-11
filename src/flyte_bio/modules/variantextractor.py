"""variantextractor — normalise a VCF's variants with variant-extractor.

Exposes :func:`variantextractor`, which reads every variant of a VCF (small
variants, SVs and breakends) with variant-extractor and writes them back,
under the input's header, as ``<prefix>.extracted.vcf``, with a
``versions.yml``. As upstream, ``args`` may hold ``--pass-only`` (keep only
PASS variants) and ``--no-ensure-pairs`` (don't require both breakends of a
pair).

The Python code is upstream's template, written out by an unquoted bash
heredoc: Groovy's ``\\\\`` escape means the same to bash, so only the dotted
placeholder is renamed (``${task.process}`` → ``${task_process}``), and the
f-string field ``{variant_extractor.__version__}`` is written with an empty
``${NOTHING}`` before its closing brace, since the shell task would read it as
a placeholder.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

VARIANTEXTRACTOR_IMAGE = "quay.io/biocontainers/variant-extractor:5.1.0--pyh106432d_0"
# Upstream's process_single.
DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory="4Gi")

# Upstream's templates/variantextractor.py, for an unquoted heredoc (see the module docstring).
VARIANTEXTRACTOR_TEMPLATE = r"""#!/usr/bin/env python3

import pysam
import variant_extractor
from variant_extractor import VariantExtractor

vcf_in = pysam.VariantFile("${vcf}")
header = str(vcf_in.header).rstrip("\\n")
vcf_in.close()

pass_only_str = "${pass_only}"
ensure_pairs_str = "${ensure_pairs}"

extractor = VariantExtractor(
    "${vcf}",
    pass_only=pass_only_str == "True",
    ensure_pairs=ensure_pairs_str == "True",
)

with open("${prefix}.extracted.vcf", "w") as out:
    out.write(header + "\\n")
    for variant_record in extractor:
        out.write(str(variant_record) + "\\n")

extractor.close()

with open("versions.yml", "w") as f:
    f.write('"${task_process}":\\n')
    f.write(f"    variant-extractor: {variant_extractor.__version__${NOTHING}}\\n")
"""

variantextractor_cmd = shell.create(
    name="variantextractor",
    image=VARIANTEXTRACTOR_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"vcf": File, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        VCF=({inputs.vcf})
        prefix={inputs.prefix}
        ARGS={inputs.args}
        task_process=VARIANTEXTRACTOR
        NOTHING=
        pass_only=False; ensure_pairs=True
        if [[ " $ARGS " == *" --pass-only "* ]]; then pass_only=True; fi
        if [[ " $ARGS " == *" --no-ensure-pairs "* ]]; then ensure_pairs=False; fi
        W=$(mktemp -d); cd "$W"
        ln -s "${VCF[0]}" .
        vcf=$(basename "${VCF[0]}")
        cat > variantextractor.py <<END_OF_TEMPLATE
"""
    + VARIANTEXTRACTOR_TEMPLATE
    + r"""END_OF_TEMPLATE
        python3 variantextractor.py || exit 1
        mv "$prefix.extracted.vcf" versions.yml {outputs.results}/
    """,
)

env = flyte.TaskEnvironment.from_task("variantextractor", variantextractor_cmd.as_task())


@dataclass
class VariantExtractorResult:
    vcf: File  # <prefix>.extracted.vcf
    versions: File  # versions.yml


async def variantextractor(vcf: File, prefix: str, args: str = "") -> VariantExtractorResult:
    """``vcf``'s variants as variant-extractor reads them: ``<prefix>.extracted.vcf``."""
    results = await variantextractor_cmd(vcf=vcf, prefix=prefix, args=args)
    extracted = await results.get_file(f"{prefix}.extracted.vcf")
    versions = await results.get_file("versions.yml")
    if extracted is None or versions is None:
        raise FileNotFoundError(f"variant-extractor wrote no {prefix}.extracted.vcf or versions.yml")
    return VariantExtractorResult(vcf=extracted, versions=versions)
