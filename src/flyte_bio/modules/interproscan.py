"""interproscan — protein domain and family annotation with InterProScan.

Exposes :func:`interproscan`, which scans a protein FASTA (plain or gzipped)
with InterProScan 5 and returns its ``<prefix>.tsv``, ``.xml``, ``.json`` and
``.gff3`` (the formats ``args`` asks for with ``--formats``; all four by
default). As upstream, the sequences are first cleaned of stop codons, which
InterProScan rejects: a trailing ``*`` is dropped and internal ones become
``X``, so coordinates still match the input.

``data`` is InterProScan's data directory (the member databases' HMMs and
models) as a ``.tar.gz``: it's tens of GiB for a real release, too many files
to stage as a Dir. It's extracted into the working directory (an archive
holding a single top-level folder is entered, as upstream's untar does), and a
copy of the container's ``interproscan.properties`` points InterProScan at it.
Without ``data``, InterProScan runs with the small sample data in the
container, which is for testing only. Pick the member databases with
``--applications`` in ``args``. A full release also needs its data's disk
(``disk=``) and more memory than the defaults: raise them with
``interproscan_cmd.override(resources=...)``.

As upstream, InterProScan runs with ``--cpu`` set to the task's CPUs and a
single CPU by default: it fails with more.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

INTERPROSCAN_IMAGE = "quay.io/biocontainers/interproscan:5.59_91.0--hec16e2b_1"
# Upstream: process_single; its tests give it 12 GB.
DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory="12Gi", disk="20Gi")

interproscan_cmd = shell.create(
    name="interproscan",
    image=INTERPROSCAN_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"fasta": File, "data": File | None, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        FA=({{inputs.fasta}})
        shopt -s nullglob; DATA_ARCHIVE=({{inputs.data}}); shopt -u nullglob
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"

        # InterProScan rejects '*': drop the trailing stop, mask internal ones as X.
        NAME=$(basename "${{FA[0]}}")
        if [ "${{NAME##*.}}" = gz ]; then
            NAME=${{NAME%.gz}}
            gzip -c -d "${{FA[0]}}" | sed '/^>/!{{s/\*$//; s/\*/X/g}}' > "$NAME"
        else
            NAME="stripped_$NAME"
            sed '/^>/!{{s/\*$//; s/\*/X/g}}' "${{FA[0]}}" > "$NAME"
        fi

        if [ ${{#DATA_ARCHIVE[@]}} -gt 0 ]; then
            mkdir unpacked
            tar -xzf "${{DATA_ARCHIVE[0]}}" -C unpacked
            DATA="$W/unpacked"
            ENTRIES=(unpacked/*)
            if [ ${{#ENTRIES[@]}} -eq 1 ] && [ -d "${{ENTRIES[0]}}" ]; then DATA="$W/${{ENTRIES[0]}}"; fi
            # Point a copy of interproscan.properties at the data.
            INTERPROSCAN_DIR="$(dirname "$(dirname "$(which interproscan.sh)")")"
            INTERPROSCAN_PROPERTIES="$(find "$INTERPROSCAN_DIR/share" -name interproscan.properties | head -n 1)"
            cp "$INTERPROSCAN_PROPERTIES" .
            sed -i \
                -e "s|^bin\.directory=.*|bin.directory=$(dirname "$INTERPROSCAN_PROPERTIES")/bin|" \
                -e "s|^data\.directory=.*|data.directory=$DATA|" \
                interproscan.properties
            # Must be absolute: interproscan.sh changes into its install directory before reading it.
            export INTERPROSCAN_CONF="$W/interproscan.properties"
        fi

        interproscan.sh \
            --cpu $TASK_CPUS \
            --input "$NAME" \
            $ARGS \
            --output-file-base "$OUT/$P" || exit 1
    """,
)

env = flyte.TaskEnvironment.from_task("interproscan", interproscan_cmd.as_task())


@dataclass
class InterproscanResult:
    tsv: File | None
    xml: File | None
    json: File | None
    gff3: File | None
    results: Dir  # everything InterProScan wrote


async def interproscan(
    fasta: File, data: File | None = None, prefix: str = "sample", args: str = ""
) -> InterproscanResult:
    """Annotate the proteins in ``fasta``, against ``data`` (a ``.tar.gz`` data directory) if given.

    Each output is None if ``--formats`` in ``args`` leaves it out.
    """
    results = await interproscan_cmd(fasta=fasta, data=data, prefix=prefix, args=args)
    return InterproscanResult(
        tsv=await results.get_file(f"{prefix}.tsv"),
        xml=await results.get_file(f"{prefix}.xml"),
        json=await results.get_file(f"{prefix}.json"),
        gff3=await results.get_file(f"{prefix}.gff3"),
        results=results,
    )
