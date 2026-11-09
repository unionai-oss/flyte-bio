"""resolvetaxonomy — one taxonomy table and bare-id sequences from either of two taxonomy sources.

Exposes :func:`resolve_taxonomy`, which takes reference ``sequences`` (FASTA,
Clustal or PHYLIP, sniffed from the content; may be gzipped) and an optional
``taxonomy`` file (a two-column TSV of id and lineage; may be gzipped). As
upstream:

- an explicit ``taxonomy`` file wins and is copied as it is, with a warning if
  the sequence headers also carry taxonomy text (anything after the id);
- otherwise the taxonomy is read from the headers, ``<id> <lineage>``. Records
  with none are an error with ``taxonomy_required``, and otherwise are left
  out of the table with a warning;
- the sequences are rewritten with headers cut down to the bare id, in their
  own format.

Outputs: ``<prefix>.resolved.tax``, ``<prefix>.resolved.<ext>`` (``ext`` is
the sequences' extension after any ``.gz``, or ``fasta`` if there is none) and
``<prefix>.warnings.txt`` (empty without warnings). Upstream's inline Python
runs in its Biopython image, unchanged but for spaces inside one f-string
field (``{ record.id }``), as the shell template reserves ``{name.name}``.
"""

from dataclasses import dataclass

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

# Upstream's Docker image: Biopython 1.84.
RESOLVE_IMAGE = "community.wave.seqera.io/library/biopython_python:1ed2958705acf8e6"
# Upstream's process_low asks for 2 CPUs; the script is single-threaded.
DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory="2Gi")

# Upstream's script, unescaped from Groovy, in a quoted heredoc.
resolve_taxonomy_cmd = shell.create(
    name="resolve_taxonomy",
    image=RESOLVE_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={"taxonomy": File | None, "sequences": File, "taxonomy_required": bool, "prefix": str, "extension": str},
    outputs={"results": Dir},
    script=r"""
        shopt -s nullglob; TAX=({inputs.taxonomy}); shopt -u nullglob
        SEQ=({inputs.sequences})
        REQUIRED={inputs.taxonomy_required}
        P={inputs.prefix}
        EXT={inputs.extension}
        OUT={outputs.results}
        cd "$OUT"
        python3 - "${TAX[0]:-}" "${SEQ[0]}" "$REQUIRED" \
            "$P.resolved.tax" "$P.resolved.$EXT" "$P.warnings.txt" << 'PYEOF' || exit 1
import gzip
import sys
from Bio import SeqIO

taxonomy_in, sequences_in, taxonomy_required, taxonomy_out, sequences_out, warnings_out = sys.argv[1:7]
taxonomy_required = taxonomy_required == 'true'

def opentext(path):
    # Sniffed from the magic bytes rather than the suffix: an input can arrive with any
    # name, and reading a gzip member as text yields mojibake instead of an error.
    with open(path, 'rb') as fh:
        compressed = fh.read(2) == b'\x1f\x8b'
    return gzip.open(path, 'rt') if compressed else open(path)

# Format sniffed from content (FASTA, Clustal or PHYLIP -- sequences can be any of
# these, depending on the caller).
with opentext(sequences_in) as fh:
    first_line = next((l.strip() for l in fh if l.strip()), '')
if first_line.startswith('>'):
    sequences_format = 'fasta'
elif first_line.upper().startswith('CLUSTAL'):
    sequences_format = 'clustal'
else:
    sequences_format = 'phylip-relaxed'

with opentext(sequences_in) as fh:
    records = list(SeqIO.parse(fh, sequences_format))
warnings = []

def embedded_taxonomy(record):
    # record.description is the *whole* header line (id + any trailing text);
    # record.id is just its first token -- GTDB's own single-file convention puts
    # the taxonomy string right after the id, space-separated.
    return record.description[len(record.id):].strip()

if taxonomy_in:
    # An explicit taxonomy file always wins. Warn (not fail) rather than silently
    # using the wrong source if the sequences also happen to carry embedded text.
    if any(embedded_taxonomy(record) for record in records):
        warnings.append(
            'An explicit taxonomy file was provided; ignoring embedded taxonomy '
            'text found in sequence record headers.'
        )
    with opentext(taxonomy_in) as fh_in, open(taxonomy_out, 'w') as fh_out:
        fh_out.write(fh_in.read())
else:
    missing = [record.id for record in records if not embedded_taxonomy(record)]
    if missing and taxonomy_required:
        sys.exit(
            'No taxonomy file was provided, and these records have no embedded '
            'taxonomy in their header either: ' + ', '.join(missing)
        )
    if missing:
        warnings.append(
            'No taxonomy file was provided; records with no embedded taxonomy in '
            'their header are omitted from the resolved taxonomy: ' + ', '.join(missing)
        )
    with open(taxonomy_out, 'w') as fh:
        for record in records:
            tax = embedded_taxonomy(record)
            if tax:
                print(f"{ record.id }\t{tax}", file=fh)

# Always strip headers down to a bare id -- some downstream tools keep the whole
# header line as the sequence/leaf name rather than just the first token.
for record in records:
    record.description = record.id
SeqIO.write(records, sequences_out, sequences_format)

with open(warnings_out, 'w') as fh:
    for warning in warnings:
        print(warning, file=fh)
PYEOF
    """,
)

env = flyte.TaskEnvironment.from_task("resolvetaxonomy", resolve_taxonomy_cmd.as_task())


def resolved_extension(sequences: File) -> str:
    """The sequences' extension after any ``.gz``, or ``fasta`` with none (upstream's naming)."""
    name = sequences.path.rstrip("/").rsplit("/", 1)[-1]
    stem = name.removesuffix(".gz") if name.endswith(".gz") else name
    return stem.rsplit(".", 1)[-1] if "." in stem else "fasta"


async def pick(results: Dir, name: str) -> File:
    found = await results.get_file(name)
    if found is None:
        raise FileNotFoundError(f"resolve_taxonomy wrote no {name}")
    return found


@dataclass
class ResolvedTaxonomy:
    taxonomy: File  # <prefix>.resolved.tax
    sequences: File  # <prefix>.resolved.<ext>, headers cut to the bare id
    warnings: File  # <prefix>.warnings.txt


async def resolve_taxonomy(
    sequences: File, taxonomy: File | None = None, taxonomy_required: bool = True, prefix: str = "sample"
) -> ResolvedTaxonomy:
    """The taxonomy of ``sequences``, from ``taxonomy`` or their headers, and the sequences with bare-id headers."""
    extension = resolved_extension(sequences)
    results = await resolve_taxonomy_cmd(
        taxonomy=taxonomy, sequences=sequences, taxonomy_required=taxonomy_required, prefix=prefix, extension=extension
    )
    return ResolvedTaxonomy(
        taxonomy=await pick(results, f"{prefix}.resolved.tax"),
        sequences=await pick(results, f"{prefix}.resolved.{extension}"),
        warnings=await pick(results, f"{prefix}.warnings.txt"),
    )
