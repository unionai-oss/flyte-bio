"""Tests for flyte_bio.modules.resolvetaxonomy.

The cases are upstream's custom/resolvetaxonomy tests (modules commit
9a1cc9f, Biopython 1.84), with upstream's assertions: content checks on the
resolved taxonomy, sequences and warnings. The extensionless and gzipped
inputs are built in the test, as upstream builds them. Upstream's failing
case checks the task's output for "taxonomy file"; here it checks that the
task fails and, when the error carries the script's message, that it's
upstream's.
"""

import gzip
import tempfile
from pathlib import Path

from flyte.io import File

from flyte_bio.modules.resolvetaxonomy import ResolvedTaxonomy, resolve_taxonomy
from tests.framework import env, fixture

EMBEDDED = "generic/fasta/resolvetaxonomy_embedded.fasta"
PLAIN = "generic/fasta/resolvetaxonomy_plain.fasta"
PARTIAL = "generic/fasta/resolvetaxonomy_partial.fasta"
OVERRIDE_TAX = "generic/tax/resolvetaxonomy_override.tax"
PLAIN_TAX = "generic/tax/resolvetaxonomy_plain.tax"


async def file_bytes(file: File) -> bytes:
    async with file.open("rb") as fh:
        return bytes(await fh.read())


async def text_of(file: File) -> str:
    return (await file_bytes(file)).decode()


async def local(name: str, data: bytes) -> File:
    path = Path(tempfile.mkdtemp()) / name
    path.write_bytes(data)
    return await File.from_local(str(path))


async def texts(result: ResolvedTaxonomy) -> tuple[str, str, str]:
    return await text_of(result.taxonomy), await text_of(result.sequences), await text_of(result.warnings)


@env.task
async def test_resolvetaxonomy_embedded_required() -> None:
    result = await resolve_taxonomy(await fixture(EMBEDDED), taxonomy_required=True, prefix="test")
    assert result.taxonomy.name == "test.resolved.tax" and result.sequences.name == "test.resolved.fasta"
    assert result.warnings.name == "test.warnings.txt", result.warnings.name
    taxonomy, sequences, warnings = await texts(result)
    assert "SeqA\tBacteria;PhylumA;ClassA;OrderA;FamilyA;GenusA;SpeciesA" in taxonomy, taxonomy
    assert "SeqB\tBacteria;PhylumA;ClassB;OrderB;FamilyB;GenusB;SpeciesB" in taxonomy, taxonomy
    assert ">SeqA\n" in sequences and ">SeqB\n" in sequences, sequences
    assert "Bacteria" not in sequences, sequences
    assert warnings.strip() == "", warnings


@env.task
async def test_resolvetaxonomy_external_file_wins() -> None:
    result = await resolve_taxonomy(
        await fixture(EMBEDDED), await fixture(OVERRIDE_TAX), taxonomy_required=True, prefix="test"
    )
    taxonomy, sequences, warnings = await texts(result)
    assert "SeqA\tArchaea;OverridePhylum;OverrideClass" in taxonomy, taxonomy
    assert "SeqB\tArchaea;OverridePhylum;OverrideClass2" in taxonomy, taxonomy
    assert "Bacteria" not in taxonomy, taxonomy
    assert ">SeqA\n" in sequences and "Bacteria" not in sequences, sequences
    assert "ignoring" in warnings.lower(), warnings


@env.task
async def test_resolvetaxonomy_external_file_no_warning() -> None:
    result = await resolve_taxonomy(
        await fixture(PLAIN), await fixture(PLAIN_TAX), taxonomy_required=True, prefix="test"
    )
    taxonomy, _, warnings = await texts(result)
    assert "SeqGood1\tBacteria;PhylumA;ClassA" in taxonomy, taxonomy
    assert warnings.strip() == "", warnings


@env.task
async def test_resolvetaxonomy_no_taxonomy_required_fails() -> None:
    try:
        await resolve_taxonomy(await fixture(PLAIN), taxonomy_required=True, prefix="test")
    except Exception as e:  # the script exits with its message on stderr
        print(f"resolve_taxonomy failed as expected: {e}")
        if "embedded" in str(e):
            assert "taxonomy file" in str(e), str(e)
    else:
        raise AssertionError("resolve_taxonomy succeeded with no taxonomy anywhere and taxonomy required")


@env.task
async def test_resolvetaxonomy_no_taxonomy_not_required_warns() -> None:
    result = await resolve_taxonomy(await fixture(PLAIN), taxonomy_required=False, prefix="test")
    taxonomy, _, warnings = await texts(result)
    assert taxonomy.strip() == "", taxonomy
    assert "SeqGood1" in warnings and "SeqGood2" in warnings, warnings


@env.task
async def test_resolvetaxonomy_partial_drops_records() -> None:
    result = await resolve_taxonomy(await fixture(PARTIAL), taxonomy_required=False, prefix="test")
    taxonomy, _, warnings = await texts(result)
    assert "SeqC\tBacteria;PhylumA;ClassA;OrderA;FamilyA;GenusA;SpeciesC" in taxonomy, taxonomy
    assert "SeqD" not in taxonomy, taxonomy
    assert "SeqD" in warnings, warnings


@env.task
async def test_resolvetaxonomy_no_extension() -> None:
    sequences = await local("resolvetaxonomy_embedded", await file_bytes(await fixture(EMBEDDED)))
    result = await resolve_taxonomy(sequences, taxonomy_required=True, prefix="test")
    assert result.sequences.name.endswith(".resolved.fasta"), result.sequences.name
    assert ">SeqA\n" in await text_of(result.sequences)


@env.task
async def test_resolvetaxonomy_gzipped_inputs() -> None:
    taxonomy_gz = await local("resolvetaxonomy_plain.tax.gz", gzip.compress(await file_bytes(await fixture(PLAIN_TAX))))
    fasta_gz = await local("resolvetaxonomy_plain.fasta.gz", gzip.compress(await file_bytes(await fixture(PLAIN))))
    result = await resolve_taxonomy(fasta_gz, taxonomy_gz, taxonomy_required=True, prefix="test")
    taxonomy, sequences, warnings = await texts(result)
    assert "SeqGood1\tBacteria;PhylumA;ClassA" in taxonomy, taxonomy
    assert warnings.strip() == "", warnings
    assert result.sequences.name.endswith(".resolved.fasta"), result.sequences.name
    assert sequences.startswith(">"), sequences[:40]


tests = [
    test_resolvetaxonomy_embedded_required,
    test_resolvetaxonomy_external_file_wins,
    test_resolvetaxonomy_external_file_no_warning,
    test_resolvetaxonomy_no_taxonomy_required_fails,
    test_resolvetaxonomy_no_taxonomy_not_required_warns,
    test_resolvetaxonomy_partial_drops_records,
    test_resolvetaxonomy_no_extension,
    test_resolvetaxonomy_gzipped_inputs,
]
