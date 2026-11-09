"""hmmer — search profile HMMs against a sequence database with hmmsearch.

Exposes :func:`hmmsearch`, which searches the profile(s) in ``hmm`` against the
sequences in ``seqdb``. Both may be gzipped: hmmsearch reads ``.gz`` inputs
itself. Outputs keep upstream's names, all gzipped with ``gzip --no-name``:

- the human-readable report: ``<prefix>.txt.gz`` (always);
- with ``write_align``: the hits' multiple alignment, ``<prefix>.sto.gz``
  (Stockholm);
- with ``write_target``: the per-target table, ``<prefix>.tbl.gz``;
- with ``write_domain``: the per-domain table, ``<prefix>.domtbl.gz``.

The report ends with run-time statistics, and the tables with the run's
directory and date, so only the alignment is byte-reproducible.

Also exposes:

- :func:`hmmalign` — align the sequences in ``fasta`` to the profile in
  ``hmm``: ``<prefix>.sto.gz`` (Stockholm, gzipped).
- :func:`hmmbuild` — build a profile named ``<prefix>`` from an alignment
  (optionally with a substitution matrix file, ``mxfile``):
  ``<prefix>.hmm.gz`` and hmmbuild's report, ``<prefix>.hmmbuild.txt``. Both
  record the run (date, CPU time), so neither is byte-reproducible.
- :func:`hmmrank` — rank each target's hits across profiles, from Parquet
  forms of hmmsearch's per-target table (and optionally its per-domain table),
  as upstream builds them with its table formatter and
  :func:`flyte_bio.modules.duckdb.duckdb_table2parquet`. One row per
  (profile, target), ranked per target by score, then e-value, profile and
  model name; with the domain table, each hit also gets the outer bounds,
  covered length and island count of the union of its domains in hmm, ali and
  env coordinates (``NA`` for a hit with no domain rows). Written with DuckDB
  as ``<prefix>.hmmrank.tsv.gz``. The tables are read, and the result written,
  under fixed names, so no name or prefix reaches the SQL. A prefix containing
  a newline is refused: Flyte's object store can't hold such a file name.
- :func:`esl_sfetch_index` — Easel's SSI index of the sequence file
  ``seqfile`` (``esl-sfetch --index``), ``<seqfile name>.ssi``, for fast
  retrieval of sequences by name. The index records the file's name, so the
  file is indexed under its own name.
- :func:`esl_alimask` — remove columns from a multiple alignment with Easel's
  ``esl-alimask`` (``args``: ``--rf-is-mask``, ``-g``, ``-p`` and the rest, or
  a mask file, ``mask``): ``<prefix>.masked.sthlm.gz``. Each ``write_*`` flag
  also writes one of esl-alimask's mask files (``--fmask-rf``,
  ``--fmask-all``, ``--gmask-rf``, ``--gmask-all``, ``--pmask-rf``,
  ``--pmask-all``) as ``<prefix>.<mask>.gz``. Everything is gzipped as
  upstream (``gzip``).
- :func:`esl_reformat` — convert a sequence file or alignment to another
  format with Easel's ``esl-reformat``. ``args`` ends with the output format
  (``afa``, ``-u phylip``, ``--gapsym=- afa``), which names the output,
  ``<prefix>.<format>.gz`` (``sequences`` with no ``args``). Upstream's
  post-processing step, a shell pipe, is here a ``sed`` script
  (``sed_script``) the output is passed through, such as ``/^>/!s/-//g`` to
  remove gaps.
- :func:`format_tsv` — HMMER's ``--tblout`` or ``--domtblout`` tables
  (``table_format``) as one TSV, ``<prefix>.<table_format>.tsv.gz``
  (``gzip -n``): a header row, then each table's rows, with a first
  ``profile`` column holding that table's label from ``labels``. The column
  count comes from each table's dash row; the free-text description (the
  rest of the line) is the last column. A table with no column-count header,
  a different layout (such as nhmmer's) or a second run's rows after the
  first run's footer is an error.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

HMMER_IMAGE = "community.wave.seqera.io/library/hmmer:3.4--ab3efa7c16f46920"
DEFAULT_RESOURCES = flyte.Resources(cpu=4, memory="8Gi")
# Upstream's process_single (hmmalign) and process_low (hmmbuild).
SINGLE_RESOURCES = flyte.Resources(cpu=1, memory="2Gi")
BUILD_RESOURCES = flyte.Resources(cpu=2, memory="4Gi")

HMMRANK_IMAGE = "community.wave.seqera.io/library/duckdb-cli:1.5.5--dca0d52c5239b678"
RANK_RESOURCES = flyte.Resources(cpu=2, memory="4Gi")

# Upstream's image for esl-sfetch --index: HMMER 3.4 (Easel 0.49) and coreutils.
SFETCH_IMAGE = "community.wave.seqera.io/library/hmmer_coreutils:e129de528776aa4c"

hmmsearch_cmd = shell.create(
    name="hmmsearch",
    image=HMMER_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "hmm": File,
        "seqdb": File,
        "prefix": str,
        "write_align": bool,
        "write_target": bool,
        "write_domain": bool,
        "args": str,
    },
    defaults={"write_align": False, "write_target": False, "write_domain": False, "args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        HMM=({{inputs.hmm}})
        SEQDB=({{inputs.seqdb}})
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"

        EXTRA=()
        TO_GZIP=("$P.txt")
        if [ {{inputs.write_align}} = true ]; then EXTRA+=(-A "$P.sto"); TO_GZIP+=("$P.sto"); fi
        if [ {{inputs.write_target}} = true ]; then EXTRA+=(--tblout "$P.tbl"); TO_GZIP+=("$P.tbl"); fi
        if [ {{inputs.write_domain}} = true ]; then EXTRA+=(--domtblout "$P.domtbl"); TO_GZIP+=("$P.domtbl"); fi

        hmmsearch \
            $ARGS \
            --cpu $TASK_CPUS \
            -o "$P.txt" \
            "${{EXTRA[@]}}" \
            "${{HMM[0]}}" \
            "${{SEQDB[0]}}"

        gzip --no-name "${{TO_GZIP[@]}}"
        mv ./*.gz "$OUT"/
    """,
)

hmmalign_cmd = shell.create(
    name="hmmalign",
    image=HMMER_IMAGE,
    resources=SINGLE_RESOURCES,
    inputs={"fasta": File, "hmm": File, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        set -o pipefail
        FASTA=({inputs.fasta})
        HMM=({inputs.hmm})
        P={inputs.prefix}
        ARGS={inputs.args}
        OUT={outputs.results}
        hmmalign $ARGS "${HMM[0]}" "${FASTA[0]}" | gzip -c > "$OUT/$P.sto.gz"
    """,
)

hmmbuild_cmd = shell.create(
    name="hmmbuild",
    image=HMMER_IMAGE,
    resources=BUILD_RESOURCES,
    inputs={"alignment": File, "mxfile": File | None, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        ALN=({{inputs.alignment}})
        shopt -s nullglob; MX=({{inputs.mxfile}}); shopt -u nullglob
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        MXOPT=()
        if [ ${{#MX[@]}} -gt 0 ]; then MXOPT=(--mxfile "${{MX[0]}}"); fi
        hmmbuild \
            $ARGS \
            --cpu $TASK_CPUS \
            -n "$P" \
            -o "$P.hmmbuild.txt" \
            "${{MXOPT[@]}}" \
            "$P.hmm" \
            "${{ALN[0]}}"
        gzip "$P.hmm"
        mv "$P.hmm.gz" "$P.hmmbuild.txt" "$OUT"/
    """,
)


def coordinate_islands_sql(coords: str) -> str:
    """Reduce one coordinate set's (hmm/ali/env) domain rows to one row per (target, profile, query).

    The outer bounds of the hit, plus the size and count of the union of its
    domains: a "merge overlapping intervals" scan. Sorted by (from, to), a
    window MAX tracks how far right the earlier rows reach; a start past that
    (plus one, as coordinates are inclusive) opens a new island. Both windows
    order by (from, to) so a tied start resolves the same way in each.
    """
    return f"""
CREATE TEMP TABLE {coords}_prev AS
SELECT accno, profile, query, {coords}_from AS f, {coords}_to AS t,
    MAX({coords}_to) OVER (
        PARTITION BY accno, profile, query ORDER BY {coords}_from, {coords}_to
        ROWS BETWEEN UNBOUNDED PRECEDING AND 1 PRECEDING
    ) AS prev_cummax
FROM dom_raw;
CREATE TEMP TABLE {coords}_isl AS
SELECT accno, profile, query, f, t,
    SUM(CASE WHEN f > COALESCE(prev_cummax, 0) + 1 THEN 1 ELSE 0 END) OVER (
        PARTITION BY accno, profile, query ORDER BY f, t
        ROWS BETWEEN UNBOUNDED PRECEDING AND CURRENT ROW
    ) AS island
FROM {coords}_prev;
CREATE TEMP TABLE {coords}_isl_agg AS
SELECT accno, profile, query, island, MIN(f) AS s, MAX(t) AS e
FROM {coords}_isl
GROUP BY accno, profile, query, island;
CREATE TEMP TABLE {coords}_coords AS
SELECT accno, profile, query,
    MIN(s) AS {coords}_from, MAX(e) AS {coords}_to,
    SUM(e - s + 1) AS {coords}_len, COUNT(*) AS {coords}_n_islands
FROM {coords}_isl_agg
GROUP BY accno, profile, query;
"""


# Rank 1 is each target's best profile; ties fall to profile, then model name.
RANKED_SQL = """
CREATE TEMP TABLE ranked AS
SELECT
    profile, target_name AS accno, query_name AS profile_desc,
    full_evalue AS evalue, full_score AS score,
    ROW_NUMBER() OVER (
        PARTITION BY target_name
        ORDER BY full_score DESC, full_evalue ASC, profile ASC, query_name ASC
    ) AS rank
FROM read_parquet('tblout.parquet');
"""

# Each hit's domain coordinates, keyed on the model (query) name as well as the
# profile, since one profile file may hold several models.
DOMAINS_SQL = (
    """
CREATE TEMP TABLE dom_raw AS
SELECT target_name AS accno, profile, query_name AS query, target_length AS tlen, query_length AS qlen,
    hmm_from, hmm_to, ali_from, ali_to, env_from, env_to
FROM read_parquet('domtblout.parquet');
"""
    + "".join(coordinate_islands_sql(coords) for coords in ("hmm", "ali", "env"))
    + """
CREATE TEMP TABLE dom_base AS
SELECT DISTINCT accno, profile, query, tlen, qlen
FROM dom_raw;
CREATE TEMP TABLE domain_coords AS
SELECT dom_base.accno, dom_base.profile, dom_base.query, dom_base.tlen, dom_base.qlen,
    hmm_coords.hmm_from, hmm_coords.hmm_to, hmm_coords.hmm_len, hmm_coords.hmm_n_islands,
    ali_coords.ali_from, ali_coords.ali_to, ali_coords.ali_len, ali_coords.ali_n_islands,
    env_coords.env_from, env_coords.env_to, env_coords.env_len, env_coords.env_n_islands
FROM dom_base
LEFT JOIN hmm_coords USING (accno, profile, query)
LEFT JOIN ali_coords USING (accno, profile, query)
LEFT JOIN env_coords USING (accno, profile, query);
"""
)

# The final ORDER BY fixes the row order, which DuckDB doesn't otherwise guarantee.
COPY_OPTIONS = "(FORMAT CSV, DELIMITER '\\t', HEADER, COMPRESSION 'gzip', NULLSTR 'NA')"
RANKS_COPY_SQL = f"COPY (SELECT * FROM ranked ORDER BY accno, rank\n) TO 'out.tsv.gz' {COPY_OPTIONS};\n"
RANKS_WITH_DOMAINS_COPY_SQL = f"""COPY (
SELECT ranked.profile, ranked.accno, ranked.profile_desc, ranked.evalue, ranked.score, ranked.rank,
    domain_coords.tlen, domain_coords.qlen,
    domain_coords.hmm_from, domain_coords.hmm_to, domain_coords.hmm_len, domain_coords.hmm_n_islands,
    domain_coords.ali_from, domain_coords.ali_to, domain_coords.ali_len, domain_coords.ali_n_islands,
    domain_coords.env_from, domain_coords.env_to, domain_coords.env_len, domain_coords.env_n_islands
FROM ranked
LEFT JOIN domain_coords
    ON ranked.accno = domain_coords.accno AND ranked.profile = domain_coords.profile
    AND ranked.profile_desc = domain_coords.query
ORDER BY ranked.accno, ranked.rank
) TO 'out.tsv.gz' {COPY_OPTIONS};
"""

# The image puts duckdb on PATH only in its entrypoint, which a shell task
# bypasses, hence the PATH line. The SQL is fixed text in quoted heredocs; only
# the limits are expanded into it.
hmmrank_cmd = shell.create(
    name="hmmrank",
    image=HMMRANK_IMAGE,
    resources=RANK_RESOURCES,
    inputs={"tblout": File, "domtblout": File | None, "prefix": str},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        export PATH="/opt/wave/.pixi/envs/default/bin:$PATH"
        TBL=({{inputs.tblout}})
        shopt -s nullglob; DOM=({{inputs.domtblout}}); shopt -u nullglob
        P={{inputs.prefix}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        ln -s "${{TBL[0]}}" tblout.parquet
        {{
            echo "SET threads=$TASK_CPUS;"
            echo "SET memory_limit='${{TASK_MEMORY_MB}}MB';"
            echo "SET temp_directory='.';"
            cat <<'SQL'
"""
    + RANKED_SQL
    + r"""SQL
            if [ ${#DOM[@]} -gt 0 ]; then
                ln -s "${DOM[0]}" domtblout.parquet
                cat <<'SQL'
"""
    + DOMAINS_SQL
    + RANKS_WITH_DOMAINS_COPY_SQL
    + r"""SQL
            else
                cat <<'SQL'
"""
    + RANKS_COPY_SQL
    + r"""SQL
            fi
        } > hmmrank.sql
        duckdb < hmmrank.sql
        mv out.tsv.gz "$OUT/$P.hmmrank.tsv.gz"
    """,
)

esl_sfetch_index_cmd = shell.create(
    name="esl_sfetch_index",
    image=SFETCH_IMAGE,
    resources=SINGLE_RESOURCES,
    inputs={"seqfile": File},
    outputs={"results": Dir},
    script=r"""
        SEQ=({inputs.seqfile})
        OUT={outputs.results}
        W=$(mktemp -d); cd "$W"
        NAME=$(basename "${SEQ[0]}")
        ln -s "${SEQ[0]}" "$NAME"
        esl-sfetch --index "$NAME" || exit 1
        mv "$NAME.ssi" "$OUT"/
    """,
)

MASKS = ("fmask-rf", "fmask-all", "gmask-rf", "gmask-all", "pmask-rf", "pmask-all")

esl_alimask_cmd = shell.create(
    name="esl_alimask",
    image=HMMER_IMAGE,
    resources=SINGLE_RESOURCES,
    inputs={"alignment": File, "mask": File | None, "masks": str, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        ALN=({inputs.alignment})
        shopt -s nullglob; MASK=({inputs.mask}); shopt -u nullglob
        MASKS={inputs.masks}
        P={inputs.prefix}
        ARGS={inputs.args}
        OUT={outputs.results}
        cd "$OUT"
        OPTS=()
        for M in $MASKS; do OPTS+=("--$M" "$P.$M"); done
        esl-alimask "${OPTS[@]}" -o "$P.masked.sthlm" $ARGS "${ALN[0]}" "${MASK[@]}" || exit 1
        gzip "$P".*mask*
    """,
)

# The sed script is passed as one argument, so no shell reaches it.
esl_reformat_cmd = shell.create(
    name="esl_reformat",
    image=HMMER_IMAGE,
    resources=SINGLE_RESOURCES,
    inputs={"seqfile": File, "filename": str, "sed_script": str, "args": str},
    defaults={"sed_script": "", "args": ""},
    outputs={"results": Dir},
    script=r"""
        set -o pipefail
        SEQ=({inputs.seqfile})
        NAME={inputs.filename}
        SED={inputs.sed_script}
        ARGS={inputs.args}
        OUT={outputs.results}
        if [ -n "$SED" ]; then
            esl-reformat $ARGS "${SEQ[0]}" | sed -e "$SED" | gzip -c > "$OUT/$NAME" || exit 1
        else
            esl-reformat $ARGS "${SEQ[0]}" | gzip -c > "$OUT/$NAME" || exit 1
        fi
    """,
)

FORMAT_TSV_IMAGE = "community.wave.seqera.io/library/gawk_gzip:fe0c202000bc8ed4"
# Upstream's process_low; the formatter streams.
FORMAT_TSV_RESOURCES = flyte.Resources(cpu=1, memory="2Gi")
TABLE_HEADERS = {
    "tblout": "profile target_name target_accession query_name query_accession full_evalue full_score full_bias "
    "best_domain_evalue best_domain_score best_domain_bias exp reg clu ov env dom rep inc description",
    "domtblout": "profile target_name target_accession target_length query_name query_accession query_length "
    "full_evalue full_score full_bias domain_number domain_count c_evalue i_evalue domain_score domain_bias "
    "hmm_from hmm_to ali_from ali_to env_from env_to acc description",
}
TABLE_COLUMNS = {"tblout": 19, "domtblout": 23}

# Upstream's formatter templates (one per format, otherwise the same), with the
# format's header, column count and name passed in. ``names`` and
# ``labels`` carry the tables (in the caller's order, as a list[File] is staged
# under its files' names and globbed alphabetically) and their labels, one per
# line.
FORMAT_TSV_AWK = r"""
            { sub(/\r$/, "") }
            N == 0 && /^#/ {
                n = 0
                ok = 1
                first = $1
                sub(/^#/, "", first)
                if (length(first) > 0) {
                    if (first ~ /^-+$/) { n++ } else { ok = 0 }
                }
                if (ok) {
                    for (i = 2; i <= NF; i++) {
                        if ($i !~ /^-+$/) { ok = 0; break }
                        n++
                    }
                }
                if (ok && n > 0) {
                    if (n != expected) {
                        print "hmmer/formattsv: expected " expected " " kind " columns, found " n " -- unsupported layout" (kind == "tblout" ? " (e.g. nhmmer)" : "") > "/dev/stderr"
                        exit 1
                    }
                    N = n
                }
                next
            }
            N > 0 && $0 == "#" { intrailer = 1; next }
            intrailer {
                if ($0 != "" && $0 !~ /^#/) {
                    print "hmmer/formattsv: data row seen after the run-metadata footer -- looks like more than one HMMER run concatenated into one file" > "/dev/stderr"
                    exit 1
                }
                next
            }
            {
                if (N == 0) {
                    print "hmmer/formattsv: data row seen before a column-count header was found" > "/dev/stderr"
                    exit 1
                }
                printf "%s", label
                for (i = 1; i < N; i++) printf "\t%s", $i
                rest = $N
                for (i = N + 1; i <= NF; i++) rest = rest " " $i
                printf "\t%s\n", rest
            }
"""  # noqa: E501

format_tsv_cmd = shell.create(
    name="hmmer_format_tsv",
    image=FORMAT_TSV_IMAGE,
    resources=FORMAT_TSV_RESOURCES,
    inputs={
        "tables": list[File],
        "names": str,
        "labels": str,
        "header": str,
        "columns": str,
        "table_format": str,
        "filename": str,
    },
    outputs={"results": Dir},
    script=r"""
        set -euo pipefail
        IN=({inputs.tables})
        NAMES={inputs.names}
        LABELS={inputs.labels}
        HEADER={inputs.header}
        EXPECTED={inputs.columns}
        KIND={inputs.table_format}
        NAME={inputs.filename}
        OUT={outputs.results}
        D=$(dirname "${IN[0]}")
        mapfile -t files <<< "$NAMES"
        mapfile -t labels <<< "$LABELS"
        {
            printf '%s\n' "$HEADER" | tr ' ' '\t'
            for i in "${!files[@]}"; do
                zcat -f "$D/${files[i]}" | awk -v label="${labels[i]}" -v expected="$EXPECTED" -v kind="$KIND" '"""
    + FORMAT_TSV_AWK
    + r"""'
            done
        } | gzip -n > "$OUT/$NAME"
    """,
)

hmmer_env = flyte.TaskEnvironment.from_task(
    "hmmer_tools",
    hmmsearch_cmd.as_task(),
    hmmalign_cmd.as_task(),
    hmmbuild_cmd.as_task(),
    esl_alimask_cmd.as_task(),
    esl_reformat_cmd.as_task(),
)
hmmrank_env = flyte.TaskEnvironment.from_task("hmmer_rank", hmmrank_cmd.as_task())
sfetch_env = flyte.TaskEnvironment.from_task("hmmer_sfetch", esl_sfetch_index_cmd.as_task())
format_tsv_env = flyte.TaskEnvironment.from_task("hmmer_format_tsv", format_tsv_cmd.as_task())
env = flyte.TaskEnvironment(name="hmmer", depends_on=[hmmer_env, hmmrank_env, sfetch_env, format_tsv_env])


@dataclass
class HmmsearchResult:
    output: File  # <prefix>.txt.gz, the human-readable report
    alignments: File | None  # <prefix>.sto.gz, with write_align
    target_summary: File | None  # <prefix>.tbl.gz, with write_target
    domain_summary: File | None  # <prefix>.domtbl.gz, with write_domain
    results: Dir  # everything hmmsearch wrote


async def hmmsearch(
    hmm: File,
    seqdb: File,
    prefix: str = "sample",
    write_align: bool = False,
    write_target: bool = False,
    write_domain: bool = False,
    args: str = "",
) -> HmmsearchResult:
    """Search the profile(s) in ``hmm`` against ``seqdb``; the ``write_*`` flags add optional outputs."""
    results = await hmmsearch_cmd(
        hmm=hmm,
        seqdb=seqdb,
        prefix=prefix,
        write_align=write_align,
        write_target=write_target,
        write_domain=write_domain,
        args=args,
    )

    async def pick(name: str) -> File:
        found = await results.get_file(name)
        if found is None:
            raise FileNotFoundError(f"hmmsearch wrote no {name}")
        return found

    return HmmsearchResult(
        output=await pick(f"{prefix}.txt.gz"),
        alignments=await pick(f"{prefix}.sto.gz") if write_align else None,
        target_summary=await pick(f"{prefix}.tbl.gz") if write_target else None,
        domain_summary=await pick(f"{prefix}.domtbl.gz") if write_domain else None,
        results=results,
    )


async def hmmalign(fasta: File, hmm: File, prefix: str = "sample", args: str = "") -> File:
    """The sequences in ``fasta`` aligned to the profile in ``hmm``: ``<prefix>.sto.gz``."""
    results = await hmmalign_cmd(fasta=fasta, hmm=hmm, prefix=prefix, args=args)
    found = await results.get_file(f"{prefix}.sto.gz")
    if found is None:
        raise FileNotFoundError(f"hmmalign wrote no {prefix}.sto.gz")
    return found


@dataclass
class HmmbuildResult:
    hmm: File  # <prefix>.hmm.gz, the profile, named <prefix>
    report: File  # <prefix>.hmmbuild.txt


async def hmmbuild(
    alignment: File, mxfile: File | None = None, prefix: str = "sample", args: str = ""
) -> HmmbuildResult:
    """A profile HMM named ``prefix`` built from ``alignment``, with hmmbuild's report."""
    results = await hmmbuild_cmd(alignment=alignment, mxfile=mxfile, prefix=prefix, args=args)

    async def pick(name: str) -> File:
        found = await results.get_file(name)
        if found is None:
            raise FileNotFoundError(f"hmmbuild wrote no {name}")
        return found

    return HmmbuildResult(hmm=await pick(f"{prefix}.hmm.gz"), report=await pick(f"{prefix}.hmmbuild.txt"))


async def hmmrank(tblout: File, domtblout: File | None = None, prefix: str = "sample") -> File:
    """Each target's hits ranked across profiles: ``<prefix>.hmmrank.tsv.gz``.

    ``tblout`` and ``domtblout`` are Parquet forms of hmmsearch's per-target
    and per-domain tables (with a ``profile`` column naming the searched
    profile); with ``domtblout``, each hit also gets its domains' coordinates.
    """
    if "\n" in prefix:
        raise ValueError("hmmrank: the prefix can't contain a newline (object stores can't hold the name)")
    results = await hmmrank_cmd(tblout=tblout, domtblout=domtblout, prefix=prefix)
    found = await results.get_file(f"{prefix}.hmmrank.tsv.gz")
    if found is None:
        raise FileNotFoundError(f"hmmrank wrote no {prefix}.hmmrank.tsv.gz")
    return found


async def esl_sfetch_index(seqfile: File) -> File:
    """Easel's SSI index of ``seqfile``: ``<seqfile name>.ssi``."""
    results = await esl_sfetch_index_cmd(seqfile=seqfile)
    indexes = [f for f in await results.list_files() if f.name.endswith(".ssi")]
    if len(indexes) != 1:
        raise FileNotFoundError(f"esl-sfetch wrote {len(indexes)} .ssi indexes, expected 1")
    return indexes[0]


async def pick_file(results: Dir, name: str, tool: str) -> File:
    found = await results.get_file(name)
    if found is None:
        raise FileNotFoundError(f"{tool} wrote no {name}")
    return found


@dataclass
class EslAlimaskResult:
    masked: File  # <prefix>.masked.sthlm.gz
    masks: dict[str, File]  # each mask asked for ("fmask-rf", ...): <prefix>.<mask>.gz


async def esl_alimask(
    alignment: File,
    mask: File | None = None,
    prefix: str = "sample",
    write_fmask_rf: bool = False,
    write_fmask_all: bool = False,
    write_gmask_rf: bool = False,
    write_gmask_all: bool = False,
    write_pmask_rf: bool = False,
    write_pmask_all: bool = False,
    args: str = "",
) -> EslAlimaskResult:
    """``alignment`` with columns removed by esl-alimask: ``<prefix>.masked.sthlm.gz``, plus any masks asked for."""
    wanted = [
        m
        for m, on in zip(
            MASKS,
            (write_fmask_rf, write_fmask_all, write_gmask_rf, write_gmask_all, write_pmask_rf, write_pmask_all),
        )
        if on
    ]
    results = await esl_alimask_cmd(alignment=alignment, mask=mask, masks=" ".join(wanted), prefix=prefix, args=args)
    return EslAlimaskResult(
        masked=await pick_file(results, f"{prefix}.masked.sthlm.gz", "esl-alimask"),
        masks={m: await pick_file(results, f"{prefix}.{m}.gz", "esl-alimask") for m in wanted},
    )


async def esl_reformat(seqfile: File, prefix: str = "sample", sed_script: str = "", args: str = "") -> File:
    """``seqfile`` in the format ``args`` ends with (optionally through ``sed_script``): ``<prefix>.<format>.gz``."""
    suffix = args.rsplit(maxsplit=1)[-1] if args.strip() else "sequences"
    filename = f"{prefix}.{suffix}.gz"
    results = await esl_reformat_cmd(seqfile=seqfile, filename=filename, sed_script=sed_script, args=args)
    return await pick_file(results, filename, "esl-reformat")


async def format_tsv(tables: list[File], labels: list[str], table_format: str, prefix: str = "sample") -> File:
    """HMMER ``tables`` as one TSV with a ``profile`` column from ``labels``: ``<prefix>.<table_format>.tsv.gz``."""
    if len(labels) != len(tables):
        raise ValueError(
            f"HMMER_FORMATTSV: labels and files must be the same length (got {len(labels)} labels, {len(tables)} files)"
        )
    if table_format not in TABLE_HEADERS:
        raise ValueError(f"HMMER_FORMATTSV: format must be 'tblout' or 'domtblout', got '{table_format}'")
    if not tables:
        raise ValueError("format_tsv needs at least one table")
    if any("\n" in label for label in labels):
        raise ValueError("format_tsv: a label can't contain a newline")
    names = [t.path.rstrip("/").rsplit("/", 1)[-1] for t in tables]
    if len(set(names)) != len(names):
        raise ValueError(f"format_tsv inputs need distinct names, got {names}")
    filename = f"{prefix}.{table_format}.tsv.gz"
    results = await format_tsv_cmd(
        tables=tables,
        names="\n".join(names),
        labels="\n".join(labels),
        header=TABLE_HEADERS[table_format],
        columns=str(TABLE_COLUMNS[table_format]),
        table_format=table_format,
        filename=filename,
    )
    return await pick_file(results, filename, "the HMMER table formatter")
