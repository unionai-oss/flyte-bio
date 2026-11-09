"""suppa — alternative splicing events and their inclusion levels with SUPPA2.

Exposes:

- :func:`suppa_generateevents` — the splicing events of a GTF. With
  ``format="ioe"``: one ``<prefix>_<event>_<strict|variable_<threshold>>.ioe``
  per event type (``SE``, ``A3``/``A5`` for ``SS``, ``MX``, ``RI``,
  ``AF``/``AL`` for ``FL``) and a ``.gtf`` of each; with ``format="ioi"``:
  ``<prefix>.ioi``, the transcripts of each gene.
- :func:`suppa_psiperevent` — the PSI of each event of an ``.ioe`` (or
  ``.ioi``) from a transcript expression table: ``<prefix>.psi``.
- :func:`suppa_psiperisoform` — the PSI of each transcript within its gene:
  ``<prefix>.psi_isoform.psi`` (SUPPA appends ``_isoform.psi`` to the output
  name ``<prefix>.psi``).
- :func:`suppa_diffsplice` — differential splicing between two conditions:
  ``<prefix>.dpsi``, ``<prefix>.psivec`` and, with ``save_tpm_events``,
  ``<prefix>_avglogtpm.tab``.
- :func:`suppa_clusterevents` — cluster events by their PSI across
  conditions: ``<prefix>.clustvec``.

As upstream, an option given a falsy value (``None``, ``0``, ``False``,
``""``) is left out, so SUPPA uses its default. ``args`` takes further
options for each tool.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

SUPPA_IMAGE = "community.wave.seqera.io/library/suppa:2.4--2612fcca3884f6bc"
# Upstream: process_low (generateEvents), process_medium (PSI) and process_high (diffSplice, clusterEvents).
LOW_RESOURCES = flyte.Resources(cpu=2, memory="4Gi")
MEDIUM_RESOURCES = flyte.Resources(cpu=2, memory="8Gi")
HIGH_RESOURCES = flyte.Resources(cpu=4, memory="16Gi")

suppa_generateevents_cmd = shell.create(
    name="suppa_generateevents",
    image=SUPPA_IMAGE,
    resources=LOW_RESOURCES,
    inputs={"gtf": File, "format": str, "prefix": str, "options": str, "args": str},
    defaults={"options": "", "args": ""},
    outputs={"results": Dir},
    script=r"""
        GTF=({inputs.gtf})
        FORMAT={inputs.format}
        P={inputs.prefix}
        OPTIONS={inputs.options}
        ARGS={inputs.args}
        cd {outputs.results}
        suppa.py generateEvents --input-file "${GTF[0]}" --format "$FORMAT" --output-file "$P" $OPTIONS $ARGS \
            || exit 1
    """,
)

suppa_psiperevent_cmd = shell.create(
    name="suppa_psiperevent",
    image=SUPPA_IMAGE,
    resources=MEDIUM_RESOURCES,
    inputs={"expression": File, "events": File, "total_filter": int, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        EXPRESSION=({inputs.expression}); EVENTS=({inputs.events})
        P={inputs.prefix}
        ARGS={inputs.args}
        cd {outputs.results}
        suppa.py psiPerEvent --ioe-file "${EVENTS[0]}" --expression-file "${EXPRESSION[0]}" \
            --total-filter {inputs.total_filter} --output-file "$P" $ARGS || exit 1
    """,
)

suppa_psiperisoform_cmd = shell.create(
    name="suppa_psiperisoform",
    image=SUPPA_IMAGE,
    resources=MEDIUM_RESOURCES,
    inputs={"expression": File, "gtf": File, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        EXPRESSION=({inputs.expression}); GTF=({inputs.gtf})
        P={inputs.prefix}
        ARGS={inputs.args}
        cd {outputs.results}
        suppa.py psiPerIsoform --expression "${EXPRESSION[0]}" --gtf "${GTF[0]}" --output-file "$P.psi" $ARGS \
            || exit 1
    """,
)

suppa_diffsplice_cmd = shell.create(
    name="suppa_diffsplice",
    image=SUPPA_IMAGE,
    resources=HIGH_RESOURCES,
    inputs={
        "events": File,
        "psi_1": File,
        "psi_2": File,
        "expression_1": File,
        "expression_2": File,
        "method": str,
        "prefix": str,
        "options": str,
        "args": str,
    },
    defaults={"options": "", "args": ""},
    outputs={"results": Dir},
    script=r"""
        EVENTS=({inputs.events}); PSI1=({inputs.psi_1}); PSI2=({inputs.psi_2})
        TPM1=({inputs.expression_1}); TPM2=({inputs.expression_2})
        METHOD={inputs.method}
        P={inputs.prefix}
        OPTIONS={inputs.options}
        ARGS={inputs.args}
        cd {outputs.results}
        suppa.py diffSplice --method "$METHOD" --psi "${PSI1[0]}" "${PSI2[0]}" --tpm "${TPM1[0]}" "${TPM2[0]}" \
            --input "${EVENTS[0]}" --output "$P" $OPTIONS $ARGS || exit 1
    """,
)

suppa_clusterevents_cmd = shell.create(
    name="suppa_clusterevents",
    image=SUPPA_IMAGE,
    resources=HIGH_RESOURCES,
    inputs={"dpsi": File, "psivec": File, "groups": str, "prefix": str, "options": str, "args": str},
    defaults={"options": "", "args": ""},
    outputs={"results": Dir},
    script=r"""
        DPSI=({inputs.dpsi}); PSIVEC=({inputs.psivec})
        GROUPS_ARG={inputs.groups}
        P={inputs.prefix}
        OPTIONS={inputs.options}
        ARGS={inputs.args}
        OUT={outputs.results}
        W=$(mktemp -d); cd "$W"
        DPSI_NAME=$(basename "${DPSI[0]}"); PSIVEC_NAME=$(basename "${PSIVEC[0]}")
        cp -L "${DPSI[0]}" "local_$DPSI_NAME"
        cp -L "${PSIVEC[0]}" "local_$PSIVEC_NAME"
        suppa.py clusterEvents --dpsi "local_$DPSI_NAME" --psivec "local_$PSIVEC_NAME" $OPTIONS \
            --groups "$GROUPS_ARG" $ARGS --output "$P" || exit 1
        mv ./*.clustvec "$OUT"/
    """,
)

env = flyte.TaskEnvironment.from_task(
    "suppa",
    suppa_generateevents_cmd.as_task(),
    suppa_psiperevent_cmd.as_task(),
    suppa_psiperisoform_cmd.as_task(),
    suppa_diffsplice_cmd.as_task(),
    suppa_clusterevents_cmd.as_task(),
)


def option(flag: str, value: object) -> str:
    """``flag value`` if ``value`` is truthy (upstream's rule), else nothing."""
    return f"{flag} {value}" if value else ""


def options(*parts: str) -> str:
    return " ".join(p for p in parts if p)


async def files_with_suffix(results: Dir, suffixes: tuple[str, ...]) -> list[File]:
    found = [f async for f in results.walk() if f.path.endswith(suffixes)]
    return sorted(found, key=lambda f: f.name)


async def pick(results: Dir, name: str, tool: str) -> File:
    found = await results.get_file(name)
    if found is None:
        raise FileNotFoundError(f"SUPPA {tool} wrote no {name}")
    return found


@dataclass
class SuppaEvents:
    events: list[File]  # the .ioe files (one per event type) or the .ioi, sorted by name
    gtfs: list[File]  # with format="ioe": a .gtf per event type, sorted by name
    results: Dir


async def suppa_generateevents(
    gtf: File,
    format: str,
    prefix: str = "events",
    pool_genes: bool = False,
    event_type: str = "",
    boundary: str = "",
    threshold: int | None = None,
    exon_length: int | None = None,
    args: str = "",
) -> SuppaEvents:
    """The splicing events of ``gtf`` (``format`` ``ioe`` or ``ioi``).

    For ``ioe``: ``event_type`` lists event types (e.g. ``"SE SS MX RI FL"``),
    ``boundary`` is ``S`` (strict) or ``V`` (variable, within ``threshold``),
    and ``exon_length`` sets the length of the exons in the event GTFs. They
    are ignored for ``ioi``, as upstream.
    """
    ioe = ""
    if format == "ioe":
        ioe = options(
            option("--event-type", event_type),
            option("--boundary", boundary),
            option("--threshold", threshold),
            option("--exon-length", exon_length),
        )
    results = await suppa_generateevents_cmd(
        gtf=gtf,
        format=format,
        prefix=prefix,
        options=options("--pool-genes" if pool_genes else "", ioe),
        args=args,
    )
    return SuppaEvents(
        events=await files_with_suffix(results, (".ioe", ".ioi")),
        gtfs=await files_with_suffix(results, (".gtf",)),
        results=results,
    )


async def suppa_psiperevent(expression: File, events: File, prefix: str, total_filter: int = 0, args: str = "") -> File:
    """The PSI of each event in ``events`` from the transcript ``expression`` table: ``<prefix>.psi``."""
    results = await suppa_psiperevent_cmd(
        expression=expression, events=events, total_filter=total_filter, prefix=prefix, args=args
    )
    return await pick(results, f"{prefix}.psi", "psiPerEvent")


async def suppa_psiperisoform(expression: File, gtf: File, prefix: str, args: str = "") -> File:
    """The PSI of each transcript of ``gtf`` within its gene: ``<prefix>.psi_isoform.psi``."""
    results = await suppa_psiperisoform_cmd(expression=expression, gtf=gtf, prefix=prefix, args=args)
    return await pick(results, f"{prefix}.psi_isoform.psi", "psiPerIsoform")


@dataclass
class SuppaDiffSplice:
    dpsi: File  # <prefix>.dpsi
    psivec: File  # <prefix>.psivec
    tpm_events: File | None  # <prefix>_avglogtpm.tab, with save_tpm_events
    results: Dir


async def suppa_diffsplice(
    events: File,
    psi_1: File,
    expression_1: File,
    psi_2: File,
    expression_2: File,
    method: str,
    prefix: str,
    local_area: int | None = None,
    lower_bound_delta_psi: float | None = None,
    is_paired: bool = False,
    gene_correction: bool = False,
    multi_testing_correction_alpha: float | None = None,
    save_tpm_events: bool = False,
    combination_analysis: bool = False,
    use_median_delta_psi: bool = False,
    tpm_threshold: float | None = None,
    nan_tpm_threshold: float | None = None,
    args: str = "",
) -> SuppaDiffSplice:
    """Differential splicing of ``events`` between condition 1 and condition 2.

    Each condition has its PSI table (:func:`suppa_psiperevent`) and its
    transcript expression table. ``method`` is ``empirical`` or ``classical``.
    """
    if method not in ("empirical", "classical"):
        raise ValueError(f"Invalid method: {method}. Must be one of 'empirical' or 'classical'.")
    flags = options(
        option("--area", local_area),
        option("--lower-bound", lower_bound_delta_psi),
        "--paired" if is_paired else "",
        "--gene-correction" if gene_correction else "",
        "--combination" if combination_analysis else "",
        option("--alpha", multi_testing_correction_alpha),
        "--save_tpm_events" if save_tpm_events else "",
        "--median" if use_median_delta_psi else "",
        option("--tpm-threshold", tpm_threshold),
        option("--nan-tpm-threshold", nan_tpm_threshold),
    )
    results = await suppa_diffsplice_cmd(
        events=events,
        psi_1=psi_1,
        psi_2=psi_2,
        expression_1=expression_1,
        expression_2=expression_2,
        method=method,
        prefix=prefix,
        options=flags,
        args=args,
    )
    return SuppaDiffSplice(
        dpsi=await pick(results, f"{prefix}.dpsi", "diffSplice"),
        psivec=await pick(results, f"{prefix}.psivec", "diffSplice"),
        tpm_events=await results.get_file(f"{prefix}_avglogtpm.tab"),
        results=results,
    )


async def suppa_clusterevents(
    dpsi: File,
    psivec: File,
    groups: str,
    prefix: str,
    significance_threshold: float | None = None,
    dpsi_threshold: float | None = None,
    maximum_distance: float | None = None,
    metric: str = "",
    separation: float | None = None,
    min_number_events: int | None = None,
    clustering_method: str = "",
    args: str = "",
) -> File:
    """Cluster the events of a :func:`suppa_diffsplice` run by their PSI: ``<prefix>.clustvec``.

    ``groups`` gives the psivec columns of each condition, e.g. ``"1-2,3-4"``;
    ``metric`` is ``euclidean``, ``manhattan`` or ``cosine``, and
    ``clustering_method`` ``DBSCAN`` or ``OPTICS``.
    """
    if metric and metric not in ("euclidean", "manhattan", "cosine"):
        raise ValueError(f"Invalid metric: {metric}. Must be one of: euclidean, manhattan, cosine")
    if clustering_method and clustering_method not in ("DBSCAN", "OPTICS"):
        raise ValueError(f"Invalid clustering method: {clustering_method}. Must be one of: DBSCAN, OPTICS")
    flags = options(
        option("--sig-threshold", significance_threshold),
        option("--dpsi-threshold", dpsi_threshold),
        option("--eps", maximum_distance),
        option("--metric", metric),
        option("--separation", separation),
        option("--min-pts", min_number_events),
        option("--clustering", clustering_method),
    )
    results = await suppa_clusterevents_cmd(
        dpsi=dpsi, psivec=psivec, groups=groups, prefix=prefix, options=flags, args=args
    )
    return await pick(results, f"{prefix}.clustvec", "clusterEvents")
