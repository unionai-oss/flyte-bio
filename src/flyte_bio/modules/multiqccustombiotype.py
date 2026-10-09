"""multiqccustombiotype — featureCounts biotype counts as MultiQC custom content.

Exposes :func:`multiqc_custom_biotype`, which turns a featureCounts table
grouped by biotype into ``<prefix>.biotype_counts_mqc.tsv`` (the bar plot)
and ``<prefix>.biotype_counts_rrna_mqc.tsv`` (% rRNA for the general stats
table). Pure Python, run as a shell task on a stock python biocontainer with
the vendored script — and, by default, upstream's MultiQC header — staged in
as File inputs.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.scripts import path

MULTIQCCUSTOMBIOTYPE_IMAGE = "quay.io/biocontainers/python:3.12.12"

DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory="2Gi")

multiqc_custom_biotype_cmd = shell.create(
    name="multiqc_custom_biotype",
    image=MULTIQCCUSTOMBIOTYPE_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"script": File, "count": File, "header": File, "prefix": str, "sample": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        SCRIPT=({inputs.script}); COUNT=({inputs.count}); HEADER=({inputs.header})
        ARGS={inputs.args}
        cd {outputs.results}
        python "${SCRIPT[0]}" --count "${COUNT[0]}" --header "${HEADER[0]}" \
            --prefix {inputs.prefix} --sample {inputs.sample} $ARGS
    """,
)


env = flyte.TaskEnvironment.from_task(
    "multiqccustombiotype",
    multiqc_custom_biotype_cmd.as_task(),
)


@dataclass
class BiotypeQC:
    counts_mqc: File  # <prefix>.biotype_counts_mqc.tsv
    rrna_mqc: File | None  # <prefix>.biotype_counts_rrna_mqc.tsv (absent when no counts)


async def multiqc_custom_biotype(
    count: File,
    prefix: str,
    sample: str | None = None,
    header: File | None = None,
    args: str = "",
) -> BiotypeQC:
    """Format a featureCounts biotype table for MultiQC (``header`` defaults to upstream's)."""
    script = await File.from_local(str(path("multiqccustombiotype.py")))
    if header is None:
        header = await File.from_local(str(path("biotypes_header.txt")))
    results = await multiqc_custom_biotype_cmd(
        script=script, count=count, header=header, prefix=prefix, sample=sample or prefix, args=args
    )
    counts_mqc = await results.get_file(f"{prefix}.biotype_counts_mqc.tsv")
    if counts_mqc is None:
        raise FileNotFoundError(f"no {prefix}.biotype_counts_mqc.tsv")
    return BiotypeQC(counts_mqc=counts_mqc, rrna_mqc=await results.get_file(f"{prefix}.biotype_counts_rrna_mqc.tsv"))
