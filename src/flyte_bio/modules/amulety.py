"""amulety — protein language model embeddings of adaptive immune receptor (BCR/TCR) sequences.

Exposes:

- :func:`amulety_translate`, which translates an AIRR rearrangement table's
  sequences to amino acids with IgBLAST (``amulety translate-igblast``),
  writing ``<prefix>_translated.tsv``. ``reference_igblast`` is an IgBLAST
  reference folder (``internal_data/``, ``optional_file/``, the germline
  databases) as a tar archive (``.tar`` or ``.tar.gz``): it holds hundreds of
  files, more than a Dir input stages reliably. An archive of a single
  top-level folder is read from that folder.
- :func:`amulety_embed`, which embeds the (translated) table's ``chain``
  (``H``, ``L``, ``HL``, ...) with a protein language model (``antiberty``,
  ``ablang``, ``esm2``, ...): ``<prefix>.tsv`` holds the embeddings and
  ``<prefix>_metadata.tsv`` the sequences' metadata. As upstream, ``prefix``
  must not name the input's own file.

Models that aren't bundled with the image are downloaded from Hugging Face at
run time, into a cache in the task's working directory (``HOME``, ``HF_HOME``
and ``TRANSFORMERS_CACHE``). Upstream runs the embedding on a GPU where there
is one; here it runs on the CPU. PyTorch's threads follow the task's CPUs.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

AMULETY_IMAGE = "community.wave.seqera.io/library/igblast_curl_python_transformers_pruned:05685e2c81024d42"
# Upstream labels translate process_single and embed process_medium (plus a GPU).
TRANSLATE_RESOURCES = flyte.Resources(cpu=1, memory="4Gi")
EMBED_RESOURCES = flyte.Resources(cpu=4, memory="16Gi", disk="20Gi")

# Writable caches for the tools and the models they download.
WRITABLE_HOME = r"""
        W=$(mktemp -d); cd "$W"
        export HOME="$W" HF_HOME="$W/cache" TRANSFORMERS_CACHE="$W/cache" MPLCONFIGDIR="$W/.config"
""".strip()

amulety_translate_cmd = shell.create(
    name="amulety_translate",
    image=AMULETY_IMAGE,
    resources=TRANSLATE_RESOURCES,
    inputs={"tsv": File, "reference_igblast": File, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        TSV=({{inputs.tsv}})
        REF_ARCHIVE=({{inputs.reference_igblast}})
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        {WRITABLE_HOME}

        mkdir reference translated
        tar -xf "${{REF_ARCHIVE[0]}}" -C reference
        set -- reference/*
        if [ $# -eq 1 ] && [ -d "$1" ]; then REF="$W/$1"; else REF="$W/reference"; fi

        export IGDATA="$REF"
        amulety \
            translate-igblast \
            --nproc $TASK_CPUS \
            $ARGS \
            --input-file "${{TSV[0]}}" \
            --output-dir translated \
            --reference-dir "$REF" || exit 1

        mv translated/*_translated.tsv "$OUT/${{P}}_translated.tsv"
    """,
)

amulety_embed_cmd = shell.create(
    name="amulety_embed",
    image=AMULETY_IMAGE,
    resources=EMBED_RESOURCES,
    inputs={"tsv": File, "chain": str, "model": str, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        TSV=({{inputs.tsv}})
        CHAIN={{inputs.chain}}
        MODEL={{inputs.model}}
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        {WRITABLE_HOME}
        export OMP_NUM_THREADS=$TASK_CPUS MKL_NUM_THREADS=$TASK_CPUS

        amulety \
            embed \
            $ARGS \
            --input-airr "${{TSV[0]}}" \
            --chain "$CHAIN" \
            --model "$MODEL" \
            --output-file-path "$P.tsv" || exit 1

        mv "$P.tsv" ./*metadata.tsv "$OUT"/
    """,
)

env = flyte.TaskEnvironment.from_task("amulety", amulety_translate_cmd.as_task(), amulety_embed_cmd.as_task())


def file_name(f: File) -> str:
    return f.path.rstrip("/").rsplit("/", 1)[-1]


async def amulety_translate(tsv: File, reference_igblast: File, prefix: str, args: str = "") -> File:
    """``tsv``'s sequences translated with IgBLAST against ``reference_igblast`` (a tar archive)."""
    results = await amulety_translate_cmd(tsv=tsv, reference_igblast=reference_igblast, prefix=prefix, args=args)
    name = f"{prefix}_translated.tsv"
    out = await results.get_file(name)
    if out is None:
        raise FileNotFoundError(f"amulety translate wrote no {name}")
    return out


@dataclass
class AmuletyEmbedding:
    embedding: File  # <prefix>.tsv
    metadata: File  # <prefix>_metadata.tsv
    results: Dir  # both


async def amulety_embed(tsv: File, chain: str, model: str, prefix: str, args: str = "") -> AmuletyEmbedding:
    """Embeddings of ``tsv``'s ``chain`` sequences with ``model``: ``<prefix>.tsv`` and its metadata table."""
    if file_name(tsv) == f"{prefix}.tsv":
        raise ValueError('Input and output names are the same, use "prefix" to disambiguate!')
    results = await amulety_embed_cmd(tsv=tsv, chain=chain, model=model, prefix=prefix, args=args)
    embedding = await results.get_file(f"{prefix}.tsv")
    if embedding is None:
        raise FileNotFoundError(f"amulety embed wrote no {prefix}.tsv")
    metadata = [f for f in await results.list_files() if file_name(f).endswith("metadata.tsv")]
    if not metadata:
        raise FileNotFoundError("amulety embed wrote no *metadata.tsv")
    return AmuletyEmbedding(embedding=embedding, metadata=metadata[0], results=results)
