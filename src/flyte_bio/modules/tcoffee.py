"""tcoffee — multiple sequence alignment, evaluation and utilities with T-Coffee.

Exposes:

- :func:`tcoffee_align` — align the sequences in ``fasta`` (``t_coffee -seq``):
  ``<prefix>.aln``, or ``<prefix>.aln.gz`` (gzipped with pigz) with
  ``compress``, in FASTA unless ``args`` sets ``-output``. With ``tree``
  (Newick) it is the guide tree (``-usetree``); with ``template`` it is the
  template file (``-template_file``) naming the ``accessory`` files (e.g. PDB
  structures). An ``-out_lib`` in ``args`` also returns the library
  (``*.*lib``).
- :func:`tcoffee_regressive` — the same with the regressive algorithm
  (``t_coffee -reg``, ``-reg_tree``, ``-reg_thread``).
- :func:`tcoffee_consensus` — the consensus of several alignments of the same
  sequences (``t_coffee -aln``, M-Coffee style): ``<prefix>.aln(.gz)``, plus
  any ``*.score_html`` / ``*.sp_ascii`` evaluations ``args`` asks for.
- :func:`tcoffee_alncompare` — compare an alignment (may be gzipped) with a
  reference (``aln_compare``, ``-compare_mode tc`` unless ``args`` sets
  another mode): ``<prefix>.scores``, a CSV of the sample's metadata and the
  score.
- :func:`tcoffee_tcs` — the TCS (transitive consistency score) of an
  alignment (may be gzipped), optionally from a T-Coffee library:
  ``<prefix>.tcs`` (the score_ascii evaluation) and ``<prefix>.scores`` (the
  metadata and the overall score, as a CSV).
- :func:`tcoffee_irmsd` — the iRMSD of an alignment (may be gzipped) given
  the structures a template file names: ``<prefix>.irmsd``.
- :func:`tcoffee_extractfrompdb` — reformat a PDB file with
  ``extract_from_pdb``: ``<prefix>.pdb``.
- :func:`tcoffee_seqreformat` — convert or analyse sequences or alignments with
  ``seq_reformat`` (``-output=...`` in ``args``): ``<prefix>.txt``.

T-Coffee writes temporary files and a cache under its home directory, so, as
upstream, each task runs in a scratch directory with ``HOME``, ``TEMP`` and
``TMP_4_TCOFFEE`` pointing there. Inputs a tool finds by name (structures named
in a template file) are linked into it. Threaded tools run with the task's
CPUs.

``args`` is split on whitespace. Where upstream's script adds ``-output
fasta_aln`` unless "``-output`` is in ``args``", its Groovy test compares whole
strings and so always adds it; here it's added only when ``args`` has no
``-output``, as intended. Upstream's alncompare appends ``-compare_mode tc``
to ``args`` with no space between; here there is one.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

# Upstream's Docker images, all T-Coffee 13.46.0.919e8c6b: with TM-align and
# pigz (align), with pigz (most tools), and alone (extractfrompdb,
# seqreformat).
TCOFFEE_TMALIGN_IMAGE = "community.wave.seqera.io/library/t-coffee_tmalign_pigz:c536fc45adb585e9"
TCOFFEE_PIGZ_IMAGE = "community.wave.seqera.io/library/t-coffee_pigz:aa6b4f1ff4adb194"
TCOFFEE_IMAGE = "community.wave.seqera.io/library/t-coffee:13.46.0.919e8c6b--6242390d545250f1"
# Upstream: align, regressive, consensus, alncompare and irmsd process_medium;
# tcs, extractfrompdb and seqreformat process_low.
DEFAULT_RESOURCES = flyte.Resources(cpu=4, memory="8Gi")
LOW_RESOURCES = flyte.Resources(cpu=2, memory="4Gi")

# Starts every script: a scratch working directory that is also T-Coffee's home.
SCRATCH = r"""
        W=$(mktemp -d); cd "$W"
        export TEMP='./'
        export TMP_4_TCOFFEE="./"
        export HOME="$W"
""".strip()

tcoffee_align_cmd = shell.create(
    name="tcoffee_align",
    image=TCOFFEE_TMALIGN_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "fasta": File,
        "tree": File | None,
        "template": File | None,
        "accessory": list[File],
        "prefix": str,
        "compress": bool,
        "args": str,
    },
    defaults={"compress": False, "args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        set -o pipefail
        FA=({{inputs.fasta}})
        shopt -s nullglob
        TREE=({{inputs.tree}}); TEMPLATE=({{inputs.template}}); ACCESSORY=({{inputs.accessory}})
        shopt -u nullglob
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        {SCRATCH}

        for F in "${{TEMPLATE[@]}}" "${{ACCESSORY[@]}}"; do ln -s "$F" "$(basename "$F")"; done
        TREE_ARGS=(); TEMPLATE_ARGS=(); DEFAULT_OUT_FORMAT=(-output fasta_aln)
        if [ ${{#TREE[@]}} -gt 0 ]; then TREE_ARGS=(-usetree "${{TREE[0]}}"); fi
        if [ ${{#TEMPLATE[@]}} -gt 0 ]; then TEMPLATE_ARGS=(-template_file "$(basename "${{TEMPLATE[0]}}")"); fi
        if [[ " $ARGS" == *" -output"* ]]; then DEFAULT_OUT_FORMAT=(); fi

        if [ {{inputs.compress}} = true ]; then
            t_coffee -seq "${{FA[0]}}" \
                "${{TREE_ARGS[@]}}" \
                "${{TEMPLATE_ARGS[@]}}" \
                $ARGS \
                "${{DEFAULT_OUT_FORMAT[@]}}" \
                -thread $TASK_CPUS \
                -outfile stdout \
                | pigz -cp $TASK_CPUS > "$P.aln.gz" || exit 1
        else
            t_coffee -seq "${{FA[0]}}" \
                "${{TREE_ARGS[@]}}" \
                "${{TEMPLATE_ARGS[@]}}" \
                $ARGS \
                "${{DEFAULT_OUT_FORMAT[@]}}" \
                -thread $TASK_CPUS \
                -outfile "$P.aln" || exit 1
        fi
        shopt -s nullglob
        mv "$P".aln* ./*.*lib "$OUT"/
    """,
)

tcoffee_regressive_cmd = shell.create(
    name="tcoffee_regressive",
    image=TCOFFEE_PIGZ_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "fasta": File,
        "tree": File | None,
        "template": File | None,
        "accessory": list[File],
        "prefix": str,
        "compress": bool,
        "args": str,
    },
    defaults={"compress": False, "args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        FA=({{inputs.fasta}})
        shopt -s nullglob
        TREE=({{inputs.tree}}); TEMPLATE=({{inputs.template}}); ACCESSORY=({{inputs.accessory}})
        shopt -u nullglob
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        {SCRATCH}

        for F in "${{TEMPLATE[@]}}" "${{ACCESSORY[@]}}"; do ln -s "$F" "$(basename "$F")"; done
        TREE_ARGS=(); TEMPLATE_ARGS=(); DEFAULT_OUT_FORMAT=(-output fasta_aln)
        if [ ${{#TREE[@]}} -gt 0 ]; then TREE_ARGS=(-reg_tree "${{TREE[0]}}"); fi
        if [ ${{#TEMPLATE[@]}} -gt 0 ]; then TEMPLATE_ARGS=(-template_file "$(basename "${{TEMPLATE[0]}}")"); fi
        if [[ " $ARGS" == *" -output"* ]]; then DEFAULT_OUT_FORMAT=(); fi
        OUTFILE="$P.aln"
        if [ {{inputs.compress}} = true ]; then OUTFILE=stdout; fi

        t_coffee -reg \
            -seq "${{FA[0]}}" \
            "${{TREE_ARGS[@]}}" \
            "${{TEMPLATE_ARGS[@]}}" \
            $ARGS \
            "${{DEFAULT_OUT_FORMAT[@]}}" \
            -reg_thread $TASK_CPUS \
            -outfile $OUTFILE || exit 1

        if [ {{inputs.compress}} = true ]; then
            pigz -cp $TASK_CPUS < stdout > "$P.aln.gz"
            rm stdout
        fi
        mv "$P".aln* "$OUT"/
    """,
)

tcoffee_consensus_cmd = shell.create(
    name="tcoffee_consensus",
    image=TCOFFEE_PIGZ_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"alignments": list[File], "tree": File | None, "prefix": str, "compress": bool, "args": str},
    defaults={"compress": False, "args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        set -o pipefail
        ALNS=({{inputs.alignments}})
        shopt -s nullglob; TREE=({{inputs.tree}}); shopt -u nullglob
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        {SCRATCH}

        TREE_ARGS=()
        if [ ${{#TREE[@]}} -gt 0 ]; then TREE_ARGS=(-usetree "${{TREE[0]}}"); fi
        if [ {{inputs.compress}} = true ]; then
            t_coffee -aln "${{ALNS[@]}}" \
                "${{TREE_ARGS[@]}}" \
                $ARGS \
                -thread $TASK_CPUS \
                -outfile stdout \
                | pigz -cp $TASK_CPUS > "$P.aln.gz" || exit 1
        else
            t_coffee -aln "${{ALNS[@]}}" \
                "${{TREE_ARGS[@]}}" \
                $ARGS \
                -thread $TASK_CPUS \
                -outfile "$P.aln" || exit 1
        fi
        shopt -s nullglob
        mv ./*.aln ./*.aln.gz ./*.score_html ./*.sp_ascii "$OUT"/
    """,
)

tcoffee_alncompare_cmd = shell.create(
    name="tcoffee_alncompare",
    image=TCOFFEE_PIGZ_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "msa": File,
        "ref_msa": File,
        "header": str,
        "values": str,
        "metric": str,
        "prefix": str,
        "args": str,
    },
    outputs={"results": Dir},
    script=rf"""
        MSA=({{inputs.msa}})
        REF=({{inputs.ref_msa}})
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        HEADER={{inputs.header}}
        VALUES={{inputs.values}}
        METRIC={{inputs.metric}}
        OUT={{outputs.results}}
        {SCRATCH}

        # check whether it is compressed
        if [[ "${{MSA[0]}}" == *.gz ]]; then
            unpigz -c "${{MSA[0]}}" > uncompressed_msa.fa
        else
            cp "${{MSA[0]}}" uncompressed_msa.fa
        fi

        t_coffee -other_pg aln_compare \
            -al1 "${{REF[0]}}" \
            -al2 uncompressed_msa.fa \
            $ARGS \
            | grep -v "seq1" | grep -v '*' | \
            awk '{{ print $4}}' ORS="\t" \
            >> "scores.txt"

        # Add metadata info to output file
        echo "$HEADER,$METRIC" > "$OUT/$P.scores"

        # Add values
        scores=$(awk '{{sub(/[[:space:]]+$/, "")}} 1' scores.txt | tr -s '[:blank:]' ',')
        echo "$VALUES,$scores" >> "$OUT/$P.scores"
    """,
)

tcoffee_tcs_cmd = shell.create(
    name="tcoffee_tcs",
    image=TCOFFEE_PIGZ_IMAGE,
    resources=LOW_RESOURCES,
    inputs={"msa": File, "lib": File | None, "header": str, "values": str, "prefix": str},
    outputs={"results": Dir},
    script=rf"""
        MSA=({{inputs.msa}})
        shopt -s nullglob; LIB=({{inputs.lib}}); shopt -u nullglob
        P={{inputs.prefix}}
        HEADER={{inputs.header}}
        VALUES={{inputs.values}}
        OUT={{outputs.results}}
        {SCRATCH}

        filename=$(basename "${{MSA[0]}}")
        if [[ "$filename" == *.gz ]] ; then
            unpigz -c "${{MSA[0]}}" > "${{filename%.gz}}"
            filename="${{filename%.gz}}"
        else
            ln -s "${{MSA[0]}}" "$filename"
        fi

        # Bad hack to circumvent t_coffee bug
        # Issue described already in: https://github.com/cbcrg/tcoffee/issues/3
        # Add an A in front of filename if the file begins with A
        first_letter_filename=${{filename:0:1}}
        if [ "$first_letter_filename" == "A" ]; then input="A"$filename; else input=$filename;  fi

        LIB_ARG=()
        if [ ${{#LIB[@]}} -gt 0 ]; then LIB_ARG=(-lib "${{LIB[0]}}"); fi
        t_coffee -infile $input \
            -evaluate -output=score_ascii \
            "${{LIB_ARG[@]}}" \
            -outfile "$P.tcs" || exit 1

        # Add metadata info to output file
        echo "$HEADER,TCS" > "$P.scores"

        # Add values
        scores=$(grep 'SCORE=' "$P.tcs" | cut -d '=' -f 2 )
        echo "$VALUES,$scores" >> "$P.scores"
        mv "$P.tcs" "$P.scores" "$OUT"/
    """,
)

tcoffee_irmsd_cmd = shell.create(
    name="tcoffee_irmsd",
    image=TCOFFEE_PIGZ_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"msa": File, "template": File, "structures": list[File], "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        MSA=({{inputs.msa}})
        TEMPLATE=({{inputs.template}})
        STRUCTURES=({{inputs.structures}})
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        {SCRATCH}

        for F in "${{TEMPLATE[0]}}" "${{STRUCTURES[@]}}"; do ln -s "$F" "$(basename "$F")"; done
        NAME=$(basename "${{MSA[0]}}")
        if [[ "$NAME" == *.gz ]] ; then
            unpigz -c "${{MSA[0]}}" > "${{NAME%.gz}}"
        else
            ln -s "${{MSA[0]}}" "$NAME"
        fi

        t_coffee -other_pg irmsd \
            "${{NAME%.gz}}" \
            $ARGS \
            -template_file "$(basename "${{TEMPLATE[0]}}")" > "$OUT/$P.irmsd" || exit 1
    """,
)

tcoffee_extractfrompdb_cmd = shell.create(
    name="tcoffee_extractfrompdb",
    image=TCOFFEE_IMAGE,
    resources=LOW_RESOURCES,
    inputs={"pdb": File, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        PDB=({{inputs.pdb}})
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        {SCRATCH}

        NAME=$(basename "${{PDB[0]}}")
        ln -s "${{PDB[0]}}" "$NAME"
        t_coffee -other_pg extract_from_pdb \
            -infile "$NAME" \
            $ARGS \
            > "$OUT/$P.pdb" || exit 1
    """,
)

tcoffee_seqreformat_cmd = shell.create(
    name="tcoffee_seqreformat",
    image=TCOFFEE_IMAGE,
    resources=LOW_RESOURCES,
    inputs={"infile": File, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        IN=({{inputs.infile}})
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        {SCRATCH}

        NAME=$(basename "${{IN[0]}}")
        ln -s "${{IN[0]}}" "$NAME"
        t_coffee -other_pg seq_reformat \
            -in "$NAME" \
            $ARGS \
            > "$OUT/$P.txt" || exit 1
    """,
)

# One image per environment: align's adds TM-align, extractfrompdb's and
# seqreformat's lack pigz.
tcoffee_tmalign_env = flyte.TaskEnvironment.from_task("tcoffee_tmalign", tcoffee_align_cmd.as_task())
tcoffee_pigz_env = flyte.TaskEnvironment.from_task(
    "tcoffee_pigz",
    tcoffee_regressive_cmd.as_task(),
    tcoffee_consensus_cmd.as_task(),
    tcoffee_alncompare_cmd.as_task(),
    tcoffee_tcs_cmd.as_task(),
    tcoffee_irmsd_cmd.as_task(),
)
tcoffee_tools_env = flyte.TaskEnvironment.from_task(
    "tcoffee_tools", tcoffee_extractfrompdb_cmd.as_task(), tcoffee_seqreformat_cmd.as_task()
)
env = flyte.TaskEnvironment(name="tcoffee", depends_on=[tcoffee_tmalign_env, tcoffee_pigz_env, tcoffee_tools_env])


@dataclass
class TcoffeeAlignResult:
    alignment: File  # <prefix>.aln(.gz)
    lib: File | None  # the T-Coffee library, when args asks for one (-out_lib)


@dataclass
class TcoffeeConsensusResult:
    alignment: File  # <prefix>.aln(.gz)
    evaluations: list[File]  # *.score_html / *.sp_ascii, when args asks for them


@dataclass
class TcoffeeTcsResult:
    tcs: File  # <prefix>.tcs: the score_ascii evaluation
    scores: File  # <prefix>.scores: the metadata and the TCS, as a CSV


def base_name(file: File) -> str:
    """The file's name without its last extension, like Nextflow's ``baseName``."""
    name = file.name
    return name.rsplit(".", 1)[0] if "." in name else name


def meta_csv(meta: dict[str, str]) -> tuple[str, str]:
    """The metadata's keys and values, each joined with commas."""
    return ",".join(meta), ",".join(str(v) for v in meta.values())


async def pick(results: Dir, name: str) -> File:
    found = await results.get_file(name)
    if found is None:
        raise FileNotFoundError(f"t_coffee wrote no {name}")
    return found


async def tcoffee_align(
    fasta: File,
    tree: File | None = None,
    template: File | None = None,
    accessory: list[File] | None = None,
    prefix: str = "sample",
    compress: bool = False,
    args: str = "",
) -> TcoffeeAlignResult:
    """Align ``fasta`` (with a guide tree, or a template naming ``accessory`` files): the alignment and any library."""
    results = await tcoffee_align_cmd(
        fasta=fasta,
        tree=tree,
        template=template,
        accessory=accessory or [],
        prefix=prefix,
        compress=compress,
        args=args,
    )
    alignment = await pick(results, f"{prefix}.aln.gz" if compress else f"{prefix}.aln")
    libs = [f async for f in results.walk() if f.name != alignment.name]
    return TcoffeeAlignResult(alignment=alignment, lib=libs[0] if libs else None)


async def tcoffee_regressive(
    fasta: File,
    tree: File | None = None,
    template: File | None = None,
    accessory: list[File] | None = None,
    prefix: str = "sample",
    compress: bool = False,
    args: str = "",
) -> File:
    """Align ``fasta`` with the regressive algorithm: ``<prefix>.aln`` or ``<prefix>.aln.gz``."""
    results = await tcoffee_regressive_cmd(
        fasta=fasta,
        tree=tree,
        template=template,
        accessory=accessory or [],
        prefix=prefix,
        compress=compress,
        args=args,
    )
    return await pick(results, f"{prefix}.aln.gz" if compress else f"{prefix}.aln")


async def tcoffee_consensus(
    alignments: list[File],
    tree: File | None = None,
    prefix: str = "sample",
    compress: bool = False,
    args: str = "",
) -> TcoffeeConsensusResult:
    """The consensus of ``alignments`` (distinct file names): the alignment and any evaluations."""
    results = await tcoffee_consensus_cmd(alignments=alignments, tree=tree, prefix=prefix, compress=compress, args=args)
    alignment = await pick(results, f"{prefix}.aln.gz" if compress else f"{prefix}.aln")
    evaluations = sorted(
        [f async for f in results.walk() if f.name.endswith((".score_html", ".sp_ascii"))], key=lambda f: f.name
    )
    return TcoffeeConsensusResult(alignment=alignment, evaluations=evaluations)


async def tcoffee_alncompare(
    msa: File,
    ref_msa: File,
    meta: dict[str, str] | None = None,
    prefix: str | None = None,
    args: str = "",
) -> File:
    """Score ``msa`` against ``ref_msa``: ``<prefix>.scores`` (prefix defaults to the msa's base name).

    The CSV's columns are ``meta``'s keys (default ``id``: the prefix) and the
    compare mode's name; its row is their values and the score.
    """
    prefix = prefix or base_name(msa)
    args = args.strip()
    if "compare_mode" not in args:
        args = f"{args} -compare_mode tc".strip()
    metric = args.split("compare_mode ")[1].split(" ")[0]
    header, values = meta_csv(meta or {"id": prefix})
    results = await tcoffee_alncompare_cmd(
        msa=msa, ref_msa=ref_msa, header=header, values=values, metric=metric, prefix=prefix, args=args
    )
    return await pick(results, f"{prefix}.scores")


async def tcoffee_tcs(
    msa: File, lib: File | None = None, meta: dict[str, str] | None = None, prefix: str = "sample"
) -> TcoffeeTcsResult:
    """The TCS of ``msa`` (from ``lib``, if given): ``<prefix>.tcs`` and ``<prefix>.scores``.

    The CSV's columns are ``meta``'s keys (default ``id``: the prefix) and
    ``TCS``.
    """
    header, values = meta_csv(meta or {"id": prefix})
    results = await tcoffee_tcs_cmd(msa=msa, lib=lib, header=header, values=values, prefix=prefix)
    return TcoffeeTcsResult(tcs=await pick(results, f"{prefix}.tcs"), scores=await pick(results, f"{prefix}.scores"))


async def tcoffee_irmsd(
    msa: File, template: File, structures: list[File], prefix: str | None = None, args: str = ""
) -> File:
    """The iRMSD of ``msa`` over the ``structures`` ``template`` names: ``<prefix>.irmsd``.

    The prefix defaults to the msa's base name.
    """
    prefix = prefix or base_name(msa)
    results = await tcoffee_irmsd_cmd(msa=msa, template=template, structures=structures, prefix=prefix, args=args)
    return await pick(results, f"{prefix}.irmsd")


async def tcoffee_extractfrompdb(pdb: File, prefix: str = "sample", args: str = "") -> File:
    """``pdb`` reformatted by ``extract_from_pdb``: ``<prefix>.pdb``."""
    results = await tcoffee_extractfrompdb_cmd(pdb=pdb, prefix=prefix, args=args)
    return await pick(results, f"{prefix}.pdb")


async def tcoffee_seqreformat(infile: File, prefix: str = "sample", args: str = "") -> File:
    """``infile`` converted or analysed by ``seq_reformat``: ``<prefix>.txt``."""
    results = await tcoffee_seqreformat_cmd(infile=infile, prefix=prefix, args=args)
    return await pick(results, f"{prefix}.txt")
