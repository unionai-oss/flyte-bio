"""humann3 — HUMAnN 3, functional profiling of metagenomes: gene families and metabolic pathways.

Exposes:

- :func:`humann3_humann`, which profiles one sample's reads against a
  nucleotide (ChocoPhlAn pangenome) and a protein (UniRef, DIAMOND) database,
  using a MetaPhlAn taxonomic profile to pick the pangenomes
  (``--taxonomic-profile``). Its outputs keep upstream's names, gzipped:
  ``<prefix>_genefamilies.tsv.gz``, ``<prefix>_pathabundance.tsv.gz``,
  ``<prefix>_pathcoverage.tsv.gz`` (when HUMAnN writes it),
  ``<prefix>_reactions.tsv.gz`` (when HUMAnN writes it) and the log
  ``<prefix>.log``. It runs on the task's CPUs (``--threads``).
- :func:`humann3_regroup`, which regroups a gene-family table (plain or
  gzipped) into other functional categories (``groups``, e.g.
  ``uniref90_level4ec``) with a utility mapping database:
  ``<prefix>_regroup.tsv.gz``. As upstream, it runs HUMAnN's regroup tool with
  the utility database set in HUMAnN's config, so no ``humann_config`` call is
  needed.
- :func:`humann3_renorm`, which renormalises a table (plain or gzipped), by
  default to copies per million: ``<prefix>_renorm.tsv.gz``.

Each database is a tar archive (any compression), extracted in the task. As
upstream, HUMAnN's databases are found inside it: the protein database is the
directory holding a ``*.dmnd`` file, the nucleotide database the directory
holding a ``*.ffn.gz`` file. A utility database archive whose entries share one
top-level directory is used from inside that directory, as upstream's untar
strips it.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

HUMANN_IMAGE = "quay.io/biocontainers/humann:3.6.1--pyh7cba7a3_0"
# Upstream labels humann process_high, regroup and renorm process_low.
HIGH_RESOURCES = flyte.Resources(cpu=8, memory="32Gi")
LOW_RESOURCES = flyte.Resources(cpu=2, memory="8Gi")

# Extracts the archive $ARCHIVE into ./$DEST and sets DB_DIR to it, or to its
# only entry when that is a directory (as upstream's untar strips it).
EXTRACT_DATABASE = r"""
        mkdir "$DEST"
        python -c 'import sys, tarfile; tarfile.open(sys.argv[1]).extractall(sys.argv[2])' "$ARCHIVE" "$DEST" || exit 1
        ENTRIES=("$DEST"/*)
        DB_DIR="$DEST"
        if [ ${#ENTRIES[@]} -eq 1 ] && [ -d "${ENTRIES[0]}" ]; then DB_DIR="${ENTRIES[0]}"; fi
"""

# Copies a table, gunzipping it if needed, to input.tsv.
STAGE_TABLE = r"""
        if [[ "$TABLE" == *.gz ]]; then gzip -c -d "$TABLE" > input.tsv; else cp "$TABLE" input.tsv; fi
"""

humann3_humann_cmd = shell.create(
    name="humann3_humann",
    image=HUMANN_IMAGE,
    resources=HIGH_RESOURCES,
    inputs={
        "input": File,
        "profile": File,
        "nucleotide_db": File,
        "protein_db": File,
        "prefix": str,
        "args": str,
    },
    defaults={"args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        INPUT=({{inputs.input}})
        PROFILE=({{inputs.profile}})
        NUC=({{inputs.nucleotide_db}})
        PROT=({{inputs.protein_db}})
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        mkdir nucleotide protein
        python -c 'import sys, tarfile; tarfile.open(sys.argv[1]).extractall(sys.argv[2])' \
            "${{NUC[0]}}" nucleotide || exit 1
        python -c 'import sys, tarfile; tarfile.open(sys.argv[1]).extractall(sys.argv[2])' \
            "${{PROT[0]}}" protein || exit 1
        PROTS=$(find -L protein -name '*.dmnd' | head -n 1)
        NUCS=$(find -L nucleotide -name '*.ffn.gz' | head -n 1)
        if [ -z "$PROTS" ] || [ -z "$NUCS" ]; then echo "no *.dmnd or *.ffn.gz in the databases" >&2; exit 1; fi
        ln -s "${{INPUT[0]}}" .
        ln -s "${{PROFILE[0]}}" .
        mkdir run; cd run

        humann \
            $ARGS \
            --threads $TASK_CPUS \
            --input "../$(basename "${{INPUT[0]}}")" \
            --protein-database "$W/$(dirname "$PROTS")" \
            --nucleotide-database "$W/$(dirname "$NUCS")" \
            --output-basename "$P" \
            --taxonomic-profile "../$(basename "${{PROFILE[0]}}")" \
            --o-log "$P.log" \
            --output . || exit 1

        gzip -n ./*.tsv
        mv ./*.tsv.gz ./*.log "$OUT"/
    """,
)

humann3_regroup_cmd = shell.create(
    name="humann3_regroup",
    image=HUMANN_IMAGE,
    resources=LOW_RESOURCES,
    inputs={"input": File, "groups": str, "utility_db": File, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        INPUT=({inputs.input})
        UTILITY=({inputs.utility_db})
        P={inputs.prefix}
        ARGS={inputs.args}
        OUT={outputs.results}
        W=$(mktemp -d); cd "$W"
        TABLE="${INPUT[0]}"
"""
    + STAGE_TABLE
    + r"""
        ARCHIVE="${UTILITY[0]}"; DEST=utility
"""
    + EXTRACT_DATABASE
    + r"""
        printf '%s\n' \
            'import os, sys' \
            'import humann.config as config' \
            'config.utility_mapping_database = os.environ["HUMANN_UTILITY_DB"]' \
            'from humann.tools.regroup_table import main' \
            'sys.exit(main())' \
            > run_regroup.py
        export HUMANN_UTILITY_DB="$W/$DB_DIR"

        python run_regroup.py \
            --input input.tsv \
            --output "${P}_regroup.tsv" \
            --groups {inputs.groups} \
            $ARGS || exit 1

        gzip -n "${P}_regroup.tsv"
        mv "${P}_regroup.tsv.gz" "$OUT"/
    """,
)

humann3_renorm_cmd = shell.create(
    name="humann3_renorm",
    image=HUMANN_IMAGE,
    resources=LOW_RESOURCES,
    inputs={"input": File, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        INPUT=({inputs.input})
        P={inputs.prefix}
        ARGS={inputs.args}
        OUT={outputs.results}
        W=$(mktemp -d); cd "$W"
        TABLE="${INPUT[0]}"
"""
    + STAGE_TABLE
    + r"""
        humann_renorm_table \
            --input input.tsv \
            --output "${P}_renorm.tsv" \
            $ARGS || exit 1

        gzip -n "${P}_renorm.tsv"
        mv "${P}_renorm.tsv.gz" "$OUT"/
    """,
)

env = flyte.TaskEnvironment.from_task(
    "humann3", humann3_humann_cmd.as_task(), humann3_regroup_cmd.as_task(), humann3_renorm_cmd.as_task()
)


@dataclass
class Humann3Result:
    genefamilies: File  # <prefix>_genefamilies.tsv.gz
    pathabundance: File  # <prefix>_pathabundance.tsv.gz
    pathcoverage: File | None  # <prefix>_pathcoverage.tsv.gz, when written
    reactions: File | None  # <prefix>_reactions.tsv.gz, when written
    log: File  # <prefix>.log
    results: Dir  # everything HUMAnN wrote


async def humann3_humann(
    input: File, profile: File, nucleotide_db: File, protein_db: File, prefix: str = "sample", args: str = ""
) -> Humann3Result:
    """Gene families and pathways of one sample's ``input`` reads, given its MetaPhlAn ``profile``."""
    results = await humann3_humann_cmd(
        input=input, profile=profile, nucleotide_db=nucleotide_db, protein_db=protein_db, prefix=prefix, args=args
    )
    return await humann3_result(results, prefix)


async def humann3_result(results: Dir, prefix: str) -> Humann3Result:
    """Pick HUMAnN's outputs out of ``results``."""
    genefamilies = await results.get_file(f"{prefix}_genefamilies.tsv.gz")
    pathabundance = await results.get_file(f"{prefix}_pathabundance.tsv.gz")
    log = await results.get_file(f"{prefix}.log")
    if genefamilies is None or pathabundance is None or log is None:
        raise FileNotFoundError(f"HUMAnN wrote no {prefix}_genefamilies.tsv.gz, _pathabundance.tsv.gz or .log")
    return Humann3Result(
        genefamilies=genefamilies,
        pathabundance=pathabundance,
        pathcoverage=await results.get_file(f"{prefix}_pathcoverage.tsv.gz"),
        reactions=await results.get_file(f"{prefix}_reactions.tsv.gz"),
        log=log,
        results=results,
    )


async def humann3_regroup(input: File, groups: str, utility_db: File, prefix: str = "sample", args: str = "") -> File:
    """``input``'s gene families regrouped into ``groups``: ``<prefix>_regroup.tsv.gz``."""
    results = await humann3_regroup_cmd(input=input, groups=groups, utility_db=utility_db, prefix=prefix, args=args)
    table = await results.get_file(f"{prefix}_regroup.tsv.gz")
    if table is None:
        raise FileNotFoundError(f"HUMAnN regroup wrote no {prefix}_regroup.tsv.gz")
    return table


async def humann3_renorm(input: File, prefix: str = "sample", args: str = "") -> File:
    """``input`` renormalised (copies per million by default): ``<prefix>_renorm.tsv.gz``."""
    results = await humann3_renorm_cmd(input=input, prefix=prefix, args=args)
    table = await results.get_file(f"{prefix}_renorm.tsv.gz")
    if table is None:
        raise FileNotFoundError(f"HUMAnN renorm wrote no {prefix}_renorm.tsv.gz")
    return table
