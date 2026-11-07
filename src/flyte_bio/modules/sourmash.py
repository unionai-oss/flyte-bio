"""sourmash — MinHash sketches of sequence files.

Exposes :func:`sourmash_sketch`, which sketches one or more sequence files
(FASTA/FASTQ, optionally gzipped) into ``<prefix>.sig``. Each file is staged as
``library/<name>``, as upstream, so signatures record that file name. With
``merge_sigs`` all files go into one signature named ``<prefix>``; otherwise
there's one signature per file, all in ``<prefix>.sig``.

``args`` sets what to sketch; the default is upstream's ``dna --param-string
scaled=1000,k=21,k=31,k=51,abund`` (its quotes dropped, as ``args`` isn't
shell-parsed).
"""

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

SOURMASH_IMAGE = "quay.io/biocontainers/sourmash:4.9.4--hdfd78af_0"
DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory="4Gi")
DEFAULT_ARGS = "dna --param-string scaled=1000,k=21,k=31,k=51,abund"

sourmash_sketch_cmd = shell.create(
    name="sourmash_sketch",
    image=SOURMASH_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"library": list[File], "merge_sigs": bool, "prefix": str, "args": str},
    defaults={"merge_sigs": False, "args": DEFAULT_ARGS},
    outputs={"results": Dir},
    script=r"""
        LIBRARY=({inputs.library})
        P={inputs.prefix}
        ARGS={inputs.args}
        OUT={outputs.results}
        W=$(mktemp -d); cd "$W"
        mkdir library
        for F in "${LIBRARY[@]}"; do ln -s "$F" library/; done
        find -L library/ -type f > library.txt
        MERGE=()
        if [ {inputs.merge_sigs} = true ]; then MERGE=(--merge "$P"); fi
        sourmash sketch \
            $ARGS \
            "${MERGE[@]}" \
            --output "$P.sig" \
            --from-file library.txt
        mv "$P.sig" "$OUT"/
    """,
)

env = flyte.TaskEnvironment.from_task("sourmash", sourmash_sketch_cmd.as_task())


async def sourmash_sketch(
    library: list[File], prefix: str = "sample", merge_sigs: bool = False, args: str = DEFAULT_ARGS
) -> File:
    """The sourmash signatures ``<prefix>.sig`` of ``library``'s files, merged into one with ``merge_sigs``."""
    results = await sourmash_sketch_cmd(library=library, merge_sigs=merge_sigs, prefix=prefix, args=args)
    signatures = await results.get_file(f"{prefix}.sig")
    if signatures is None:
        raise FileNotFoundError(f"sourmash wrote no {prefix}.sig")
    return signatures
