"""Tests for flyte_bio.modules.krakenuniq.

The cases and expected md5s are upstream's krakenuniq/build tests (modules
commit 9a1cc9f, KrakenUniq 1.0.4), with the test config's options (its
``--jellyfish-bin`` is passed by the module). Upstream builds over the NCBI
taxonomy that ``krakenuniq-download`` fetches; these tests use the test
data's SARS-CoV-2 protein taxonomy (``prot_names.dmp``, ``prot_nodes.dmp``:
the genome's lineage) instead, renamed ``names.dmp`` and ``nodes.dmp``. As
upstream, the build log, report and ``taxDB`` are checked by name only (they
hold timestamps, paths and the taxonomy's labels), every other file by md5.

The classification cases are upstream's krakenuniq/preloadedkrakenuniq tests
(same commit), against upstream's test database archive
(``genomics/sarscov2/genome/db/krakenuniq.tar.gz``), passed as the archive:
two samples per batch, named ``sample_1`` and ``sample.2``. As upstream, the
reports are checked by name (they hold a timestamp) and the per-read
classifications by md5. Upstream's snapshots have no read outputs: its
compression step (``find . -maxdepth 0``) never matches, so its reads stay
uncompressed and its ``*.gz`` outputs are empty. Here the reads are gzipped
and returned, so the tests check that each sample's are there. Upstream also
greps its ``.command.sh`` for where the options went; that has no
counterpart here, so those cases check the outputs only.
"""

import shutil
import tempfile
from pathlib import Path

from flyte.io import Dir, File

from flyte_bio.modules.krakenuniq import KrakenuniqResult, krakenuniq_build, krakenuniq_preloadedkrakenuniq
from tests.framework import assert_md5, env, fixture

METAGENOME = "genomics/sarscov2/metagenome/"
ARGS = "--work-on-disk --max-db-size 12 --kmer-len 15 --minimizer-len 13"
UNSTABLE = ["database-build.log", "database.report.tsv", "taxDB"]
STABLE_MD5S = {
    "database.idx": "f5298aa4215c956d2f737b5988721a24",
    "database.kdb": "46542be19979e867075f9e455a314e00",
    "database.kdb.counts": "a82149480a435210ec1e870f06b6cdb3",
    "database0.kdb": "46542be19979e867075f9e455a314e00",
}
INTERMEDIATE_MD5S = {
    "database.jdb": "48db5f6bb9b9e9db7f5bbdd0d8bec4a8",
    "database.kraken.tsv": "be9bd0cb3a8bea30989e65aa42a73727",
    "library-files.txt": "8a9817d2de8465eb436aa3aa8696a717",
}


async def krakenuniq_renamed(path: str, name: str) -> File:
    """The fixture at ``path``, uploaded again under ``name``."""
    local = Path(tempfile.mkdtemp()) / name
    shutil.copy(await (await fixture(path)).download(), local)
    return await File.from_local(str(local))


async def krakenuniq_test_db(keep_intermediate: bool) -> Dir:
    return await krakenuniq_build(
        await fixture("genomics/sarscov2/genome/genome.fasta"),
        taxonomy=[
            await krakenuniq_renamed(METAGENOME + "prot_names.dmp", "names.dmp"),
            await krakenuniq_renamed(METAGENOME + "prot_nodes.dmp", "nodes.dmp"),
        ],
        seqid2taxid=await fixture(METAGENOME + "seqid2taxid.map"),
        keep_intermediate=keep_intermediate,
        prefix="test",
        args=ARGS,
    )


async def check_krakenuniq_db(db: Dir, stable: dict[str, str], unstable: list[str]) -> None:
    files = {f.name: f async for f in db.walk(recursive=True)}
    names = sorted(files)
    assert names == sorted([*stable, *unstable]), f"krakenuniq db: {names}"
    for name, md5 in stable.items():
        await assert_md5(files[name], md5, label=name)


@env.task
async def test_krakenuniq_build() -> None:
    # upstream case: krakenuniq/build "sarscov2 - fasta"
    await check_krakenuniq_db(await krakenuniq_test_db(False), STABLE_MD5S, ["taxDB"])


@env.task
async def test_krakenuniq_build_nocleanup() -> None:
    # upstream case: krakenuniq/build "sarscov2 - fasta - nocleanup"
    await check_krakenuniq_db(await krakenuniq_test_db(True), STABLE_MD5S | INTERMEDIATE_MD5S, UNSTABLE)


SARSCOV2 = "genomics/sarscov2/"
PREFIXES = ["sample_1", "sample.2"]
FASTA_MD5S = {"sample_1": "2bea6c2195c400a909a2d4cca2e3045e", "sample.2": "8aafacd89a6aac98aaf512df0a7493d1"}


async def krakenuniq_classify(
    sequences: list[str], sequence_type: str, paired: bool = False, save: bool = True, args: str = "", args2: str = ""
) -> KrakenuniqResult:
    return await krakenuniq_preloadedkrakenuniq(
        [await fixture(SARSCOV2 + path) for path in sequences],
        PREFIXES,
        await fixture(SARSCOV2 + "genome/db/krakenuniq.tar.gz"),
        sequence_type=sequence_type,
        paired=paired,
        save_output_reads=save,
        report_file=save,
        save_output=save,
        args=args,
        args2=args2,
    )


async def check_krakenuniq_classified(
    result: KrakenuniqResult, md5s: dict[str, str], sequence_type: str, stem: str = ""
) -> None:
    reports = sorted(f.name for f in result.report)
    assert reports == sorted(f"{p}.krakenuniq.report.txt" for p in PREFIXES), reports
    assignments = {f.name: f for f in result.classified_assignment}
    assert sorted(assignments) == sorted(f"{p}.krakenuniq.classified.txt" for p in md5s), sorted(assignments)
    for prefix, md5 in md5s.items():
        name = f"{prefix}.krakenuniq.classified.txt"
        await assert_md5(assignments[name], md5, label=name)
    for kind, files in (("classified", result.classified_reads), ("unclassified", result.unclassified_reads)):
        names = sorted(f.name for f in files)
        assert names == sorted(f"{p}{stem}.{kind}.{sequence_type}.gz" for p in PREFIXES), f"{kind} reads: {names}"


FASTA = ["illumina/fasta/contigs.fasta", "genome/genome.fasta.gz"]


@env.task
async def test_krakenuniq_preloadedkrakenuniq_fasta() -> None:
    # upstream case: krakenuniq/preloadedkrakenuniq "sarscov2 - FASTA"
    await check_krakenuniq_classified(await krakenuniq_classify(FASTA, "fasta"), FASTA_MD5S, "fasta")


@env.task
async def test_krakenuniq_preloadedkrakenuniq_fasta_ramchunksize() -> None:
    # upstream case: krakenuniq/preloadedkrakenuniq "sarscov2 - FASTA - ramchunksize"
    result = await krakenuniq_classify(FASTA, "fasta", args2="--preload-size 1GB")
    await check_krakenuniq_classified(result, FASTA_MD5S, "fasta")


@env.task
async def test_krakenuniq_preloadedkrakenuniq_fasta_extra_args() -> None:
    # upstream case: krakenuniq/preloadedkrakenuniq "sarscov2 - FASTA - extra args"
    result = await krakenuniq_classify(FASTA, "fasta", args="--hll-precision 14", args2="--exact")
    await check_krakenuniq_classified(result, FASTA_MD5S, "fasta")


@env.task
async def test_krakenuniq_preloadedkrakenuniq_fasta_no_optional_outputs() -> None:
    # upstream case: krakenuniq/preloadedkrakenuniq "sarscov2 - FASTA - no optional outputs"
    result = await krakenuniq_classify(FASTA, "fasta", save=False)
    assert result.classified_reads == [] and result.unclassified_reads == []
    assert result.classified_assignment == [] and result.report == []


@env.task
async def test_krakenuniq_preloadedkrakenuniq_fastq_single() -> None:
    # upstream case: krakenuniq/preloadedkrakenuniq "sarscov2 - Illumina FASTQ single"
    result = await krakenuniq_classify(
        ["illumina/fastq/test_interleaved.fastq.gz", "illumina/fastq/test2_1.fastq.gz"], "fastq"
    )
    md5s = {"sample_1": "88a734a9a9216cb0770a77f36c9f4e78", "sample.2": "f885fa1cdbfc5460af0772219991bf6d"}
    await check_krakenuniq_classified(result, md5s, "fastq")


@env.task
async def test_krakenuniq_preloadedkrakenuniq_fastq_paired() -> None:
    # upstream case: krakenuniq/preloadedkrakenuniq "sarscov2 - Illumina FASTQ paired-end"
    fastqs = ["test_1.fastq.gz", "test_2.fastq.gz", "test2_1.fastq.gz", "test2_2.fastq.gz"]
    result = await krakenuniq_classify([f"illumina/fastq/{f}" for f in fastqs], "fastq", paired=True)
    md5 = "ed5e19c7a88312cc04e483ac5f2579cd"
    await check_krakenuniq_classified(result, {"sample_1": md5, "sample.2": md5}, "fastq", ".merged")


tests = [
    test_krakenuniq_build,
    test_krakenuniq_build_nocleanup,
    test_krakenuniq_preloadedkrakenuniq_fasta,
    test_krakenuniq_preloadedkrakenuniq_fasta_ramchunksize,
    test_krakenuniq_preloadedkrakenuniq_fasta_extra_args,
    test_krakenuniq_preloadedkrakenuniq_fasta_no_optional_outputs,
    test_krakenuniq_preloadedkrakenuniq_fastq_single,
    test_krakenuniq_preloadedkrakenuniq_fastq_paired,
]
