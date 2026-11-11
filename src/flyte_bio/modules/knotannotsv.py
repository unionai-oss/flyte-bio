"""knotannotsv — knotAnnotSV reports of AnnotSV structural-variant annotations.

Exposes :func:`knotannotsv`, which turns an AnnotSV annotation table (the
``.tsv`` from :func:`flyte_bio.modules.annotsv.annotsv_annotsv`) into an
interactive HTML report (``knotAnnotSV.pl``), or with ``excel`` an Excel
workbook with macros (``knotAnnotSV2XL.pl``): ``<prefix>.html`` or
``<prefix>.xlsm``.

As upstream, the tool names its output after its input, so the table is
copied to ``<prefix>.tsv`` first; the configuration is the one bundled with
the tool (``share/knotAnnotSV/config_AnnotSV.yaml``). Options such as
``--genomeBuild hg38`` go in ``args``.
"""

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

KNOTANNOTSV_IMAGE = "quay.io/biocontainers/knotannotsv:1.1.5--hdfd78af_0"

# Upstream labels it process_medium; the tool is single-threaded.
DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory="6Gi")

knotannotsv_cmd = shell.create(
    name="knotannotsv",
    image=KNOTANNOTSV_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"annotsv_tsv": File, "excel": bool, "prefix": str, "args": str},
    defaults={"excel": False, "args": ""},
    outputs={"results": Dir},
    script=r"""
        TSV=({inputs.annotsv_tsv})
        P={inputs.prefix}
        ARGS={inputs.args}
        W=$(mktemp -d); cd "$W"
        if [ {inputs.excel} = true ]; then SCRIPT=knotAnnotSV2XL.pl; EXT=xlsm; else SCRIPT=knotAnnotSV.pl; EXT=html; fi
        CONFIG="${CONDA_PREFIX:-/usr/local}/share/knotAnnotSV/config_AnnotSV.yaml"
        # The tool names its output after its input (sample.tsv -> sample.html), so the input takes the prefix.
        cp "${TSV[0]}" "$P.tsv"
        $SCRIPT \
            $ARGS \
            --configFile "$CONFIG" \
            --annotSVfile "$P.tsv" || exit 1
        mv "$P.$EXT" {outputs.results}/
    """,
)

env = flyte.TaskEnvironment.from_task("knotannotsv", knotannotsv_cmd.as_task())


async def knotannotsv(annotsv_tsv: File, prefix: str, excel: bool = False, args: str = "") -> File:
    """Report ``annotsv_tsv`` as ``<prefix>.html``, or ``<prefix>.xlsm`` with ``excel``."""
    results = await knotannotsv_cmd(annotsv_tsv=annotsv_tsv, excel=excel, prefix=prefix, args=args)
    name = f"{prefix}.{'xlsm' if excel else 'html'}"
    report = await results.get_file(name)
    if report is None:
        raise FileNotFoundError(f"knotAnnotSV wrote no {name}")
    return report
