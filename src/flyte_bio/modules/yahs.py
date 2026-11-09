"""yahs — YaHS, scaffolding an assembly with Hi-C contacts.

Exposes :func:`yahs`, which scaffolds a contig assembly (FASTA with its
``.fai``) with Hi-C reads aligned to it (BAM, BED or PA5), optionally
starting from an existing AGP (``-a``). As upstream, the outputs keep YaHS'
names under ``-o <prefix>``:

- ``<prefix>_scaffolds_final.fa`` and ``<prefix>_scaffolds_final.agp``;
- the initial-break AGP ``<prefix>_inital_break_*.agp`` (YaHS' spelling) or
  ``<prefix>_no_break*.agp``, and the per-round AGPs ``<prefix>_r*_*.agp``;
- the contact map ``<prefix>.bin`` and YaHS' stderr, ``<prefix>.log``.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

YAHS_IMAGE = "community.wave.seqera.io/library/yahs:1.2.2--5a61538ab6dbd63d"
DEFAULT_RESOURCES = flyte.Resources(cpu=2, memory="16Gi")

yahs_cmd = shell.create(
    name="yahs",
    image=YAHS_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"fasta": File, "fai": File, "hic_map": File, "agp": File | None, "prefix": str, "args": str},
    defaults={"args": ""},
    outputs={"results": Dir},
    script=r"""
        FA=({inputs.fasta}); FAI=({inputs.fai}); MAP=({inputs.hic_map})
        shopt -s nullglob; AGP=({inputs.agp}); shopt -u nullglob
        P={inputs.prefix}
        ARGS={inputs.args}
        OUT={outputs.results}
        W=$(mktemp -d); cd "$W"
        # The FASTA needs its index beside it; inputs go in their own directory, apart from the outputs.
        mkdir in
        FASTA=in/$(basename "${FA[0]}")
        ln -s "${FA[0]}" "$FASTA"
        ln -s "${FAI[0]}" "$FASTA.fai"
        AGP_INPUT=()
        if [ ${#AGP[@]} -gt 0 ]; then AGP_INPUT=(-a "${AGP[0]}"); fi
        yahs -o "$P" "${AGP_INPUT[@]}" $ARGS "$FASTA" "${MAP[0]}" 2> "$P.log" || { tail -n 20 "$P.log" >&2; exit 1; }
        cat "$P.log" >&2
        shopt -s nullglob
        mv "${P}"_scaffolds_final.fa "${P}"_scaffolds_final.agp "${P}"_inital_break*.agp "${P}"_no_break*.agp \
            "${P}"_r*_*.agp "$P.bin" "$P.log" "$OUT"/
    """,
)

env = flyte.TaskEnvironment.from_task("yahs", yahs_cmd.as_task())


@dataclass
class YahsResult:
    scaffolds_fasta: File | None  # <prefix>_scaffolds_final.fa
    scaffolds_agp: File | None  # <prefix>_scaffolds_final.agp
    initial_break_agp: list[File]  # <prefix>_{inital,no}_break*.agp
    round_agp: list[File]  # <prefix>_r*_*.agp
    binary: File  # <prefix>.bin
    log: File  # <prefix>.log
    results: Dir  # everything above


async def yahs(
    fasta: File, fai: File, hic_map: File, agp: File | None = None, prefix: str = "sample", args: str = ""
) -> YahsResult:
    """Scaffold ``fasta`` with the Hi-C alignments in ``hic_map``."""
    results = await yahs_cmd(fasta=fasta, fai=fai, hic_map=hic_map, agp=agp, prefix=prefix, args=args)
    files = {f.name: f async for f in results.walk()}
    for name in (f"{prefix}.bin", f"{prefix}.log"):
        if name not in files:
            raise FileNotFoundError(f"YaHS wrote no {name}")
    final = f"{prefix}_scaffolds_final"
    breaks = (f"{prefix}_inital_break", f"{prefix}_no_break")
    return YahsResult(
        scaffolds_fasta=files.get(f"{final}.fa"),
        scaffolds_agp=files.get(f"{final}.agp"),
        initial_break_agp=[files[n] for n in sorted(files) if n.startswith(breaks) and n.endswith(".agp")],
        round_agp=[files[n] for n in sorted(files) if n.startswith(f"{prefix}_r") and n.endswith(".agp")],
        binary=files[f"{prefix}.bin"],
        log=files[f"{prefix}.log"],
        results=results,
    )
