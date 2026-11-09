"""kofamscan — KEGG Orthology assignment with KofamScan.

Exposes :func:`kofamscan`, which assigns KEGG Orthology (KO) numbers to the
proteins of a FASTA with ``exec_annotation``: HMM searches against the KOfam
profiles, scored by the per-KO thresholds in ``ko_list``. It writes
``<prefix>.tsv`` when ``args`` has ``--format detail-tsv``, else
``<prefix>.txt`` (KofamScan's ``detail`` table, or whichever ``--format``
``args`` names).

``profiles`` is a ``.tar.gz`` of the profile directory, as KEGG distributes it
(``profiles.tar.gz``, which unpacks to ``profiles/`` with one ``.hmm`` per KO
and the ``prokaryote.hal`` and ``eukaryote.hal`` subsets). It's an archive
rather than a Dir because the full set has tens of thousands of files, more
than a task's input staging copes with. ``hal`` picks a subset in it (e.g.
``prokaryote.hal``); by default every profile is searched. ``ko_list`` may be
gzipped, as KEGG distributes it. KofamScan runs with ``--cpu`` of the task's
CPUs.
"""

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

KOFAMSCAN_IMAGE = "quay.io/biocontainers/kofamscan:1.3.0--hdfd78af_2"
# Upstream's process_medium. The full profile set is a few GiB unpacked.
DEFAULT_RESOURCES = flyte.Resources(cpu=4, memory="8Gi", disk="20Gi")

kofamscan_cmd = shell.create(
    name="kofamscan",
    image=KOFAMSCAN_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"fasta": File, "profiles": File, "ko_list": File, "hal": str, "out_name": str, "args": str},
    defaults={"hal": "", "args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        FA=({{inputs.fasta}})
        PROFILES_ARCHIVE=({{inputs.profiles}})
        KO_LIST=({{inputs.ko_list}})
        HAL={{inputs.hal}}
        NAME={{inputs.out_name}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        mkdir profiles_db
        tar -xzf "${{PROFILES_ARCHIVE[0]}}" -C profiles_db || exit 1
        # The profile directory: the archive's only top-level directory, else its root.
        ENTRIES=(profiles_db/*)
        PROFILES=profiles_db
        if [ ${{#ENTRIES[@]}} -eq 1 ] && [ -d "${{ENTRIES[0]}}" ]; then PROFILES="${{ENTRIES[0]}}"; fi
        if [ -n "$HAL" ]; then PROFILES="$PROFILES/$HAL"; fi
        if [ "${{KO_LIST[0]##*.}}" = gz ]; then
            gzip -c -d "${{KO_LIST[0]}}" > ko_list
            KO=ko_list
        else
            KO="${{KO_LIST[0]}}"
        fi
        exec_annotation \
            -p "$PROFILES" \
            -k "$KO" \
            $ARGS \
            --cpu $TASK_CPUS \
            -o "$OUT/$NAME" \
            "${{FA[0]}}" || exit 1
    """,
)

env = flyte.TaskEnvironment.from_task("kofamscan", kofamscan_cmd.as_task())


async def kofamscan(
    fasta: File, profiles: File, ko_list: File, prefix: str = "sample", hal: str = "", args: str = ""
) -> File:
    """KO assignments for the proteins in ``fasta``: ``<prefix>.tsv`` with ``--format detail-tsv``, else ``.txt``."""
    name = f"{prefix}.tsv" if "--format detail-tsv" in args else f"{prefix}.txt"
    results = await kofamscan_cmd(fasta=fasta, profiles=profiles, ko_list=ko_list, hal=hal, out_name=name, args=args)
    out = await results.get_file(name)
    if out is None:
        raise FileNotFoundError(f"kofamscan wrote no {name}")
    return out
