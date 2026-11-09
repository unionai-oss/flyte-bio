"""kmcp — k-mer-based metagenomic profiling databases.

Exposes the two steps of building a KMCP database, as upstream:

- :func:`kmcp_compute` computes k-mers of reference sequences (``kmcp
  compute``): FASTA/FASTQ files, linked into ``genomes/``, or a directory of
  them used as ``genomes/``. It writes one ``.unik`` file per reference and a
  summary, ``_info.txt``, into ``<prefix>/``; as a database holds many
  references, that directory comes back as ``<prefix>.tar.gz``, beside
  ``_info.txt``.
- :func:`kmcp_index` indexes that archive (``kmcp index``) into
  ``<prefix>/`` (``R001/__db.yml``, ``__name_mapping.tsv``,
  ``_block*.uniki``, ...), returned as ``<prefix>.tar.gz`` with the log
  ``<prefix>.log``. Given taxonomy dump files and/or a sequence-to-taxon map,
  it copies them into ``<prefix>/taxdmp/`` for ``kmcp profile``. As
  ``_info.txt`` names the ``.unik`` files by the compute directory's name,
  the archive is extracted under that name, so the two prefixes must differ.

``--threads`` is the task's CPUs; options go in ``args``.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

KMCP_IMAGE = "quay.io/biocontainers/kmcp:0.9.4--h9ee0642_0"
# Upstream labels both process_medium.
DEFAULT_RESOURCES = flyte.Resources(cpu=4, memory="8Gi")

kmcp_compute_cmd = shell.create(
    name="kmcp_compute",
    image=KMCP_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"sequences": list[File], "sequences_dir": Dir | None, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        shopt -s nullglob; SEQUENCES=({{inputs.sequences}}); shopt -u nullglob
        SEQUENCES_DIR={{inputs.sequences_dir}}
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        # An unset Dir input renders as a path that doesn't exist (or is empty).
        if [ -d "$SEQUENCES_DIR" ] && [ -n "$(ls -A "$SEQUENCES_DIR")" ]; then
            ln -s "$SEQUENCES_DIR" genomes
        else
            mkdir genomes
            for F in "${{SEQUENCES[@]}}"; do ln -s "$F" genomes/"$(basename "$F")"; done
        fi
        kmcp \
            compute \
            $ARGS \
            --threads $TASK_CPUS \
            --out-dir "$P"/ \
            --in-dir genomes/ \
            || exit 1
        cp "$P/_info.txt" "$OUT"/
        tar -czf "$OUT/$P.tar.gz" "$P"
    """,
)

kmcp_index_cmd = shell.create(
    name="kmcp_index",
    image=KMCP_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"compute": File, "taxdmp": list[File], "seq2taxidmap": list[File], "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        COMPUTE=({{inputs.compute}})
        shopt -s nullglob; TAXDMP=({{inputs.taxdmp}} {{inputs.seq2taxidmap}}); shopt -u nullglob
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        # _info.txt names the .unik files by the compute directory's own name, so keep it.
        tar -xzf "${{COMPUTE[0]}}"
        COMPUTE_DIR=$(tar -tzf "${{COMPUTE[0]}}" | awk -F/ 'NR == 1 {{ print $1 }}')
        kmcp \
            index \
            --in-dir "$COMPUTE_DIR" \
            $ARGS \
            --threads $TASK_CPUS \
            --log "$P.log" \
            --out-dir "$P" \
            || exit 1
        # Optionally copy over taxonomy files if provided for downstream module
        if [ ${{#TAXDMP[@]}} -gt 0 ]; then
            mkdir "$P/taxdmp"
            for F in "${{TAXDMP[@]}}"; do cp "$F" "$P/taxdmp/$(basename "$F")"; done
        fi
        mv "$P.log" "$OUT"/
        tar -czf "$OUT/$P.tar.gz" "$P"
    """,
)

env = flyte.TaskEnvironment.from_task("kmcp", kmcp_compute_cmd.as_task(), kmcp_index_cmd.as_task())


@dataclass
class KmcpComputeResult:
    archive: File  # <prefix>.tar.gz: the directory <prefix>/ of .unik files and _info.txt
    info: File  # _info.txt


async def kmcp_compute(sequences: File | list[File] | Dir, prefix: str = "sample", args: str = "") -> KmcpComputeResult:
    """Compute k-mers of ``sequences`` (files, or a directory of them)."""
    if isinstance(sequences, Dir):
        files, directory = [], sequences
    else:
        files, directory = (sequences if isinstance(sequences, list) else [sequences]), None
    results = await kmcp_compute_cmd(sequences=files, sequences_dir=directory, prefix=prefix, args=args)
    archive = await results.get_file(f"{prefix}.tar.gz")
    info = await results.get_file("_info.txt")
    if archive is None or info is None:
        raise FileNotFoundError(f"kmcp compute wrote no {prefix}.tar.gz or _info.txt")
    return KmcpComputeResult(archive=archive, info=info)


@dataclass
class KmcpIndexResult:
    archive: File  # <prefix>.tar.gz: the database directory <prefix>/
    log: File  # <prefix>.log


async def kmcp_index(
    compute: File,
    taxdmp: list[File] | None = None,
    seq2taxidmap: File | None = None,
    prefix: str = "sample",
    args: str = "",
) -> KmcpIndexResult:
    """Index :func:`kmcp_compute`'s archive; taxonomy files are copied into ``<prefix>/taxdmp/``."""
    if compute.name == f"{prefix}.tar.gz":
        raise ValueError(f"kmcp index's prefix {prefix!r} must differ from kmcp compute's")
    results = await kmcp_index_cmd(
        compute=compute,
        taxdmp=taxdmp or [],
        seq2taxidmap=[seq2taxidmap] if seq2taxidmap is not None else [],
        prefix=prefix,
        args=args,
    )
    archive = await results.get_file(f"{prefix}.tar.gz")
    log = await results.get_file(f"{prefix}.log")
    if archive is None or log is None:
        raise FileNotFoundError(f"kmcp index wrote no {prefix}.tar.gz or {prefix}.log")
    return KmcpIndexResult(archive=archive, log=log)
