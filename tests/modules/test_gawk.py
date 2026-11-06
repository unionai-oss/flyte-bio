"""Tests for flyte_bio.modules.gawk.

The cases and expected md5s are upstream's gawk tests (modules commit
9a1cc9f, gawk 5.3.1). Upstream's ``ext.args2`` programs are shell-quoted; here
they are the bare awk program. nf-test hashes ``.gz`` outputs decompressed.
"""

import tempfile
from pathlib import Path

from flyte.io import File

from flyte_bio.modules.gawk import GawkResult, gawk
from tests.framework import assert_gunzipped_md5, assert_md5, env, fixture

FAI = "genomics/homo_sapiens/genome/genome.fasta.fai"
TO_BED = r'BEGIN { FS = OFS = "\t"}; { print $1, "0", $2 }'
FIRST_COLUMN = r'BEGIN {FS=" "}; {print $1}'
HEADER_ONLY = r'BEGIN {OFS="\t"} END {print "chrom", "start", "end"}'
BED_MD5 = "87a15eb9c2ff20ccd5cd8735a28708f7"


async def program_file(text: str, name: str = "program.awk") -> File:
    path = Path(tempfile.mkdtemp()) / name
    path.write_text(text + "\n")
    return await File.from_local(str(path))


async def text_files() -> list[File]:
    return [await fixture("generic/txt/hello.txt"), await fixture("generic/txt/species_names.txt")]


def single_output(result: GawkResult, name: str) -> File:
    names = [f.path.rsplit("/", 1)[-1] for f in result.outputs]
    assert names == [name], f"expected outputs [{name}], got {names}"
    assert result.output is not None
    return result.output


async def expect_error(message: str, **kwargs) -> None:
    try:
        await gawk(**kwargs)
    except ValueError as e:
        assert message in str(e), f"expected {message!r}, got {e!r}"
    else:
        raise AssertionError(f"gawk did not fail with {message!r}")


@env.task
async def test_gawk_fasta_to_bed() -> None:
    result = await gawk([await fixture(FAI)], program=TO_BED, prefix="test", suffix="bed")
    await assert_md5(single_output(result, "test.bed"), BED_MD5, label="test.bed")


@env.task
async def test_gawk_fasta_to_bed_program_file() -> None:
    result = await gawk([await fixture(FAI)], program_file=await program_file(TO_BED), prefix="test", suffix="bed")
    await assert_md5(single_output(result, "test.bed"), BED_MD5, label="test.bed")


@env.task
async def test_gawk_fasta_to_bed_awk_redirect() -> None:
    program = r'BEGIN { FS = OFS = "\t"}; { print $1, "0", $2 > "test.bed" }'
    result = await gawk([await fixture(FAI)], program=program, prefix="test", suffix="bed", redirect_output=False)
    await assert_md5(single_output(result, "test.bed"), BED_MD5, label="test.bed")


@env.task
async def test_gawk_first_column_multiple_files() -> None:
    files = await text_files()
    result = await gawk(files, program_file=await program_file(FIRST_COLUMN), prefix="test", suffix="bed")
    await assert_md5(single_output(result, "test.bed"), "566c51674bd643227bb2d83e0963376d", label="test.bed")


@env.task
async def test_gawk_unzip_inputs() -> None:
    vcfs = "genomics/homo_sapiens/illumina/vcf/"
    files = [await fixture(vcfs + "NA12878_chrM.vcf.gz"), await fixture(vcfs + "NA24385_sv.vcf.gz")]
    program = await program_file("/^#CHROM/ { print $1, $10 }", "column_header.awk")
    result = await gawk(files, program_file=program, prefix="test", suffix="bed")
    await assert_md5(single_output(result, "test.bed"), "1e31ebd4a060aab5433bbbd9ab24e403", label="test.bed")


@env.task
async def test_gawk_compress_output() -> None:
    result = await gawk([await fixture(FAI)], program=TO_BED, prefix="test", suffix="txt.gz")
    await assert_gunzipped_md5(single_output(result, "test.txt.gz"), BED_MD5, label="test.txt.gz")


@env.task
async def test_gawk_input_and_output_names_clash() -> None:
    await expect_error(
        "Input and output names are the same, set prefix in module configuration to disambiguate!",
        files=await text_files(),
        program_file=await program_file(FIRST_COLUMN),
        prefix="hello",
        suffix="txt",
    )


@env.task
async def test_gawk_suffix_from_input_extension() -> None:
    result = await gawk([await fixture(FAI)], program=TO_BED, prefix="test")
    await assert_md5(single_output(result, "test.fai"), BED_MD5, label="test.fai")


@env.task
async def test_gawk_no_input_with_suffix() -> None:
    result = await gawk([], program_file=await program_file(HEADER_ONLY), prefix="hello", suffix="tsv")
    await assert_md5(single_output(result, "hello.tsv"), "08153dd97bae3a3cc0f1433d01940d7d", label="hello.tsv")


@env.task
async def test_gawk_no_input_no_suffix() -> None:
    await expect_error(
        "suffix not set and no input provided",
        files=[],
        program_file=await program_file(HEADER_ONLY),
        prefix="hello",
    )


tests = [
    test_gawk_fasta_to_bed,
    test_gawk_fasta_to_bed_program_file,
    test_gawk_fasta_to_bed_awk_redirect,
    test_gawk_first_column_multiple_files,
    test_gawk_unzip_inputs,
    test_gawk_compress_output,
    test_gawk_input_and_output_names_clash,
    test_gawk_suffix_from_input_extension,
    test_gawk_no_input_with_suffix,
    test_gawk_no_input_no_suffix,
]
