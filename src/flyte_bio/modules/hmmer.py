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

hmmer_env = flyte.TaskEnvironment.from_task(
    "hmmer_tools", hmmsearch_cmd.as_task(), hmmalign_cmd.as_task(), hmmbuild_cmd.as_task()
)
hmmrank_env = flyte.TaskEnvironment.from_task("hmmer_rank", hmmrank_cmd.as_task())
sfetch_env = flyte.TaskEnvironment.from_task("hmmer_sfetch", esl_sfetch_index_cmd.as_task())
env = flyte.TaskEnvironment(name="hmmer", depends_on=[hmmer_env, hmmrank_env, sfetch_env])


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
