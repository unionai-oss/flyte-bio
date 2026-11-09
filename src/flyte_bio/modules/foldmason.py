"""foldmason — multiple protein structure alignment with FoldMason.

Exposes:

- :func:`foldmason_createdb` — a FoldMason (Foldseek-format) database from
  structure files: a Dir of ``<prefix>``, ``<prefix>.{dbtype,index,lookup,source}``
  and ``<prefix>_{ca,h,ss}{,.dbtype,.index}``.
- :func:`foldmason_easymsa` — align structures: ``<prefix>_3di.fa`` (the 3Di
  alignment) and ``<prefix>_aa.fa`` (amino acids), each gzipped with pigz
  with ``compress``. With ``tree`` (Newick), FoldMason uses it as the guide
  tree (``--guide-tree``).
- :func:`foldmason_msa2lddtreport` — an HTML report of an alignment's LDDT
  scores against the structures in a database: ``<prefix>.html``.

All run with ``--threads`` set to the task's CPUs; FoldMason's options go in
``args``.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

# Upstream's Docker images: FoldMason, and FoldMason with pigz for easy-msa.
FOLDMASON_IMAGE = "community.wave.seqera.io/library/foldmason:4.dd3c235--375019dcc3c5fe6f"
FOLDMASON_PIGZ_IMAGE = "community.wave.seqera.io/library/foldmason_pigz:289a5a7a978299df"
# Upstream: createdb process_low, easymsa and msa2lddtreport process_medium.
# The CPUs match what upstream's tests run with.
CREATEDB_RESOURCES = flyte.Resources(cpu=2, memory="4Gi")
DEFAULT_RESOURCES = flyte.Resources(cpu=4, memory="8Gi")

foldmason_createdb_cmd = shell.create(
    name="foldmason_createdb",
    image=FOLDMASON_IMAGE,
    resources=CREATEDB_RESOURCES,
    inputs={"structures": list[File], "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"db": Dir},
    script=rf"""
        {TASK_LIMITS}
        STRUCTURES=({{inputs.structures}})
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.db}}
        W=$(mktemp -d); cd "$W"

        for F in "${{STRUCTURES[@]}}"; do ln -s "$F" "$(basename "$F")"; done
        foldmason createdb \
            "${{STRUCTURES[@]##*/}}" \
            "$P" \
            $ARGS \
            --threads $TASK_CPUS || exit 1
        mv "$P"* "$OUT"/
    """,
)

foldmason_easymsa_cmd = shell.create(
    name="foldmason_easymsa",
    image=FOLDMASON_PIGZ_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"structures": list[File], "tree": File | None, "prefix": str, "compress": bool, "args": str},
    defaults={"compress": False, "args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        STRUCTURES=({{inputs.structures}})
        shopt -s nullglob; TREE=({{inputs.tree}}); shopt -u nullglob
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"

        for F in "${{STRUCTURES[@]}}"; do ln -s "$F" "$(basename "$F")"; done
        OPTIONS_TREE=()
        if [ ${{#TREE[@]}} -gt 0 ]; then OPTIONS_TREE=(--guide-tree "${{TREE[0]}}"); fi
        foldmason easy-msa \
            "${{STRUCTURES[@]##*/}}" \
            "$P" \
            tmp \
            "${{OPTIONS_TREE[@]}}" \
            $ARGS \
            --threads $TASK_CPUS || exit 1

        if [ {{inputs.compress}} = true ]; then
            pigz -p $TASK_CPUS "$P"_3di.fa
            pigz -p $TASK_CPUS "$P"_aa.fa
        fi
        mv "$P"_3di.fa* "$P"_aa.fa* "$OUT"/
    """,
)

foldmason_msa2lddtreport_cmd = shell.create(
    name="foldmason_msa2lddtreport",
    image=FOLDMASON_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "msa": File,
        "db": Dir,
        "db_name": str,
        "structures": list[File],
        "tree": File | None,
        "prefix": str,
        "args": str,
    },
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        MSA=({{inputs.msa}})
        DB={{inputs.db}}
        STRUCTURES=({{inputs.structures}})
        shopt -s nullglob; TREE=({{inputs.tree}}); shopt -u nullglob
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"

        # As upstream, the database's files and the structures sit side by side
        # in the working directory, and the database is named by its base name.
        for F in "$DB"/* "${{STRUCTURES[@]}}"; do ln -s "$F" "$(basename "$F")"; done
        OPTIONS_TREE=()
        if [ ${{#TREE[@]}} -gt 0 ]; then OPTIONS_TREE=(--guide-tree "${{TREE[0]}}"); fi
        foldmason msa2lddtreport \
            {{inputs.db_name}} \
            "${{MSA[0]}}" \
            "$OUT/$P.html" \
            $ARGS \
            "${{OPTIONS_TREE[@]}}" \
            --threads $TASK_CPUS || exit 1
    """,
)

# One image per environment: easy-msa's adds pigz.
foldmason_env = flyte.TaskEnvironment.from_task(
    "foldmason_tools", foldmason_createdb_cmd.as_task(), foldmason_msa2lddtreport_cmd.as_task()
)
foldmason_pigz_env = flyte.TaskEnvironment.from_task("foldmason_pigz", foldmason_easymsa_cmd.as_task())
env = flyte.TaskEnvironment(name="foldmason", depends_on=[foldmason_env, foldmason_pigz_env])


@dataclass
class FoldmasonMsaResult:
    msa_3di: File  # <prefix>_3di.fa(.gz): the alignment in the 3Di alphabet
    msa_aa: File  # <prefix>_aa.fa(.gz): the alignment in amino acids


async def foldmason_createdb(structures: File | list[File], prefix: str = "sample", args: str = "") -> Dir:
    """A FoldMason database of ``structures`` (PDB or mmCIF): a Dir of ``<prefix>*``."""
    if isinstance(structures, File):
        structures = [structures]
    return await foldmason_createdb_cmd(structures=structures, prefix=prefix, args=args)


async def foldmason_easymsa(
    structures: list[File],
    tree: File | None = None,
    prefix: str = "sample",
    compress: bool = False,
    args: str = "",
) -> FoldmasonMsaResult:
    """Align ``structures`` (with ``tree`` as the guide tree, if given): the 3Di and amino-acid alignments."""
    results = await foldmason_easymsa_cmd(structures=structures, tree=tree, prefix=prefix, compress=compress, args=args)
    suffix = ".gz" if compress else ""

    async def pick(name: str) -> File:
        found = await results.get_file(name)
        if found is None:
            raise FileNotFoundError(f"foldmason wrote no {name}")
        return found

    return FoldmasonMsaResult(
        msa_3di=await pick(f"{prefix}_3di.fa{suffix}"), msa_aa=await pick(f"{prefix}_aa.fa{suffix}")
    )


async def foldmason_msa2lddtreport(
    msa: File,
    db: Dir,
    db_name: str,
    structures: list[File],
    tree: File | None = None,
    prefix: str = "sample",
    args: str = "",
) -> File:
    """LDDT report of ``msa`` against the structures in ``db`` (named ``db_name``): ``<prefix>.html``."""
    results = await foldmason_msa2lddtreport_cmd(
        msa=msa, db=db, db_name=db_name, structures=structures, tree=tree, prefix=prefix, args=args
    )
    found = await results.get_file(f"{prefix}.html")
    if found is None:
        raise FileNotFoundError(f"foldmason wrote no {prefix}.html")
    return found
