"""Tests for flyte_bio.modules.bioawk.

The cases and expected md5s are upstream's bioawk tests (modules commit
9a1cc9f, bioawk 1.0). Upstream's ``ext.args2`` programs are shell-quoted; here
they are what the shell hands bioawk. In the "tsv" case upstream's quoting
turns ``OFS = '\\t'`` into ``OFS = t`` (an unset variable), which is kept.
nf-test hashes ``.gz`` outputs decompressed.
"""

import tempfile
from pathlib import Path

from flyte.io import File

from flyte_bio.modules.bioawk import BioawkResult, bioawk
from tests.framework import assert_gunzipped_md5, assert_md5, env, fixture

GENOME = "genomics/sarscov2/genome/genome.fasta"
TELOMERES = "generic/telomere/telomeric_motifs.fasta"
TELOMERE_MD5 = "45a26fc19df7fafe489fb66ee057b069"


async def program_file(text: str) -> File:
    path = Path(tempfile.mkdtemp()) / "program.awk"
    path.write_text(text)
    return await File.from_local(str(path))


def single_output(result: BioawkResult, name: str) -> File:
    names = [f.path.rsplit("/", 1)[-1] for f in result.outputs]
    assert names == [name], f"expected outputs [{name}], got {names}"
    return result.outputs[0]


@env.task
async def test_bioawk_fasta_to_fasta_gz() -> None:
    result = await bioawk(
        [await fixture(GENOME)],
        "fa.gz",
        program='{print ">" $name ORS length($seq)}',
        prefix="sample_1",
        args="-c fastx",
    )
    await assert_gunzipped_md5(single_output(result, "sample_1.fa.gz"), "b558dd15d8940373a032a827d490e693")


@env.task
async def test_bioawk_fasta_to_tsv() -> None:
    result = await bioawk(
        [await fixture(TELOMERES)],
        "tsv",
        program='BEGIN {OFS = t} {print toupper($seq), gsub(/G/, "", $seq)} ',
        prefix="telomere_summary",
        args="-c fastx",
    )
    await assert_md5(single_output(result, "telomere_summary.tsv"), TELOMERE_MD5, label="telomere_summary.tsv")


@env.task
async def test_bioawk_fasta_to_tsv_no_redirect() -> None:
    program = await program_file('BEGIN {OFS = "\t"} {print toupper($seq), gsub(/G/, "", $seq) > "output.tsv"}')
    result = await bioawk(
        [await fixture(TELOMERES)],
        "tsv",
        program_file=program,
        prefix="telomere_summary",
        redirect_output=False,
        args="-c fastx",
    )
    assert result.output is None, "no-redirect run wrote telomere_summary.tsv"
    await assert_md5(single_output(result, "output.tsv"), TELOMERE_MD5, label="output.tsv")


@env.task
async def test_bioawk_fasta_to_tsv_program_file() -> None:
    program = await program_file('BEGIN {OFS = "\t"} {print toupper($seq), gsub(/G/, "", $seq) }')
    result = await bioawk(
        [await fixture(TELOMERES)], "tsv", program_file=program, prefix="telomere_summary", args="-c fastx"
    )
    await assert_md5(single_output(result, "telomere_summary.tsv"), TELOMERE_MD5, label="telomere_summary.tsv")


tests = [
    test_bioawk_fasta_to_fasta_gz,
    test_bioawk_fasta_to_tsv,
    test_bioawk_fasta_to_tsv_no_redirect,
    test_bioawk_fasta_to_tsv_program_file,
]
