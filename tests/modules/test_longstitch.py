"""Tests for flyte_bio.modules.longstitch.

The cases and expected md5s are upstream's longstitch tests (modules commit
9a1cc9f, LongStitch 1.0.5). As upstream, outputs it lists as unstable (the
scaffold graph, LINKS' graph and log, ARKS' graph and, for test data #2, the
tig-pair checkpoints) are checked by name. ``run`` and ``tigmint-ntLink``
write the same outputs, as upstream's snapshots show.

The thread count changes Tigmint's corrections; upstream's outputs match 2
threads (not 1, 3, 4, 6 or 8), so the tests run on 2 CPUs.

Upstream emits no final ``<prefix>.tigmint-ntLink.fa`` or
``<prefix>.ntLink-arks.fa`` (its shell tests for them are malformed; see the
module); this port writes them, so beside upstream's outputs the tests allow
those and check the pipeline's final scaffolds exist.
"""

import flyte

from flyte_bio.modules.longstitch import longstitch_cmd, read_longstitch_results
from tests.framework import assert_md5, assert_nonempty, env, fixture

DATA = "genomics/scaffolding/longstitch/"
T = "assembly.cut250.tigmint.fa.k32.w100"
N = "test.k32.w100"
A = "ntLink.scaffolds_c4_m8-10000_cut250_k20_r0.05_e30000_z1000"
FINAL = {"test.tigmint-ntLink.fa", "test.ntLink-arks.fa", "test.tigmint-ntLink-arks.fa"}

TIGMINT_1 = {
    f"{T}.tsv": "d168051e4b208350f1dd4a497a2e59e3",
    f"{T}.z1000.n2.scaffold.dot": None,
    f"{T}.z1000.ntLink.scaffolds.fa": "16be0be0da4a0a506972f9d5e59a66aa",
    f"{T}.z1000.stitch.abyss-scaffold.fa": "16be0be0da4a0a506972f9d5e59a66aa",
    f"{T}.z1000.stitch.path": "c54af479a9a1ce52a99eca6bd3306295",
    f"{T}.z1000.trimmed_scafs.agp": "4972d12e493703a4c157abcb6b8d2793",
    f"{T}.z1000.trimmed_scafs.fa": "addcad088a26382fb25282de508a1662",
    f"{T}.z1000.trimmed_scafs.path": "0e351d3339f764ff7a03cb5e62e86c92",
    f"{T}.z1000.trimmed_scafs.tsv": "727434abb839c06e6cf3889ba77e4f4b",
    f"{T}.z1000.verbose_mapping.tsv": "58475b43cd8dfc255024f275ad719f5f",
}
TIGMINT_ARKS_1 = {
    f"{T}.tsv": "d168051e4b208350f1dd4a497a2e59e3",
    f"{T}.z1000.n2.scaffold.dot": None,
    f"{T}.z1000.ntLink.scaffolds.fa": "16be0be0da4a0a506972f9d5e59a66aa",
    f"{T}.z1000.ntLink.scaffolds.renamed.fa": "be2718cc1762de765f001757a341825c",
    f"{T}.z1000.{A}.dist.gv": "fc21a3b342ed3ba968304ae62fe914f5",
    f"{T}.z1000.{A}.tigpair_checkpoint.tsv": "d69bada77f0306911aa2de9f1dc247a5",
    f"{T}.z1000.{A}_l4_a0.3.assembly_correspondence.tsv": "df66f8c6e47f3a0a0d748e3c73de7f89",
    f"{T}.z1000.{A}_l4_a0.3.gv": None,
    f"{T}.z1000.{A}_l4_a0.3.log": None,
    f"{T}.z1000.{A}_l4_a0.3.scaffolds": "e03cc0724ef0b27aab91e46fd2ff5086",
    f"{T}.z1000.{A}_l4_a0.3.scaffolds.fa": "b118dbe22d0f77890426dedd6cb62eee",
    f"{T}.z1000.{A}_l4_a0.3.tigpair_checkpoint.tsv": "d69bada77f0306911aa2de9f1dc247a5",
    f"{T}.z1000.{A}_main.tsv": "1de8b67483640b84c42ca1056e103b2e",
    f"{T}.z1000.{A}_original.gv": None,
    f"{T}.z1000.stitch.abyss-scaffold.fa": "16be0be0da4a0a506972f9d5e59a66aa",
    f"{T}.z1000.stitch.path": "c54af479a9a1ce52a99eca6bd3306295",
    f"{T}.z1000.trimmed_scafs.agp": "4972d12e493703a4c157abcb6b8d2793",
    f"{T}.z1000.trimmed_scafs.fa": "addcad088a26382fb25282de508a1662",
    f"{T}.z1000.trimmed_scafs.path": "0e351d3339f764ff7a03cb5e62e86c92",
    f"{T}.z1000.trimmed_scafs.tsv": "727434abb839c06e6cf3889ba77e4f4b",
    f"{T}.z1000.verbose_mapping.tsv": "58475b43cd8dfc255024f275ad719f5f",
    "test.tigmint-ntLink-arks.fa": "5190838b0d62f2a2c93aec4a17545841",
}
NTLINK_ARKS_1 = {
    f"{N}.tsv": "3df675c11d2e74134e7c7fb7d1e0cb30",
    f"{N}.z1000.n2.scaffold.dot": None,
    f"{N}.z1000.ntLink.scaffolds.renamed.fa": "62dc7eeee6276c1aca40844af0f2e92e",
    f"{N}.z1000.{A}.dist.gv": "bcc8de2496ddff7d1eef0bac1afff6fc",
    f"{N}.z1000.{A}.tigpair_checkpoint.tsv": "15b34804e0f465ca0cafe0f7e67b6543",
    f"{N}.z1000.{A}_l4_a0.3.assembly_correspondence.tsv": "d9e0e50914ff20eb311939d77e987f06",
    f"{N}.z1000.{A}_l4_a0.3.gv": None,
    f"{N}.z1000.{A}_l4_a0.3.log": None,
    f"{N}.z1000.{A}_l4_a0.3.scaffolds": "c049358e40a982f54e2a6dc2d37bab30",
    f"{N}.z1000.{A}_l4_a0.3.scaffolds.fa": "ff4779050147fcf5e9e7bf07097b2376",
    f"{N}.z1000.{A}_main.tsv": "f21e5fd6be86a75483f60c15b70b793a",
    f"{N}.z1000.{A}_original.gv": None,
    f"{N}.z1000.stitch.abyss-scaffold.fa": "a42d63d43b054a0e980f76a5335904d0",
    f"{N}.z1000.stitch.path": "6087209d2321cc7e1e4a664c98f8b7cc",
    f"{N}.z1000.trimmed_scafs.agp": "266d8c22189e32ca840e32f15ea253a7",
    f"{N}.z1000.trimmed_scafs.fa": "71a2fd8e043fd7661e899ac02815438d",
    f"{N}.z1000.trimmed_scafs.path": "bb18307ab27d705b8df58e2d0c22548d",
    f"{N}.z1000.trimmed_scafs.tsv": "ec53bc9a919cf452508d40d724b8f381",
    f"{N}.z1000.verbose_mapping.tsv": "a91f681a7bedfc89a73e1ea3f0367f28",
}
TIGMINT_2 = {
    f"{T}.tsv": "7ccd5619d889ae6ce018efc11008f613",
    f"{T}.z1000.n2.scaffold.dot": None,
    f"{T}.z1000.ntLink.scaffolds.fa": "fe8b860019ab275a923033ffddec7ecf",
    f"{T}.z1000.stitch.abyss-scaffold.fa": "fe8b860019ab275a923033ffddec7ecf",
    f"{T}.z1000.stitch.path": "1287e5876d465e16e445557bff5015af",
    f"{T}.z1000.trimmed_scafs.agp": "ab8e9c15b4b9c155752ba28015f245ce",
    f"{T}.z1000.trimmed_scafs.fa": "157af72f0cc0c4bdfbd07f77a0a49bdf",
    f"{T}.z1000.trimmed_scafs.path": "7fb5b5d5c07dc70d6bfaa7747e4f1887",
    f"{T}.z1000.trimmed_scafs.tsv": "523925e6dbde0d7c48dcbe2d42303b01",
    f"{T}.z1000.verbose_mapping.tsv": "9ef46b00e351e589e4bb3b728bca36e0",
}
TIGMINT_ARKS_2 = {
    f"{T}.tsv": "7ccd5619d889ae6ce018efc11008f613",
    f"{T}.z1000.n2.scaffold.dot": None,
    f"{T}.z1000.ntLink.scaffolds.fa": "fe8b860019ab275a923033ffddec7ecf",
    f"{T}.z1000.ntLink.scaffolds.renamed.fa": "a31035b82753802aa1b207540df17f78",
    f"{T}.z1000.{A}.dist.gv": "b1a6eb5399ec9f06ac732927574cdd4b",
    f"{T}.z1000.{A}.tigpair_checkpoint.tsv": None,
    f"{T}.z1000.{A}_l4_a0.3.assembly_correspondence.tsv": "8bd8dcbd2f8159a65a4bdb6637a41e14",
    f"{T}.z1000.{A}_l4_a0.3.gv": None,
    f"{T}.z1000.{A}_l4_a0.3.log": None,
    f"{T}.z1000.{A}_l4_a0.3.scaffolds": "b3bbf43b5314c01fd7f78dd52aa24f76",
    f"{T}.z1000.{A}_l4_a0.3.scaffolds.fa": "b5507c20bd520e0a8ee31e63d388623d",
    f"{T}.z1000.{A}_l4_a0.3.tigpair_checkpoint.tsv": None,
    f"{T}.z1000.{A}_main.tsv": "dba0904ccd4241cac377c7ba55ffa168",
    f"{T}.z1000.{A}_original.gv": None,
    f"{T}.z1000.stitch.abyss-scaffold.fa": "fe8b860019ab275a923033ffddec7ecf",
    f"{T}.z1000.stitch.path": "1287e5876d465e16e445557bff5015af",
    f"{T}.z1000.trimmed_scafs.agp": "ab8e9c15b4b9c155752ba28015f245ce",
    f"{T}.z1000.trimmed_scafs.fa": "157af72f0cc0c4bdfbd07f77a0a49bdf",
    f"{T}.z1000.trimmed_scafs.path": "7fb5b5d5c07dc70d6bfaa7747e4f1887",
    f"{T}.z1000.trimmed_scafs.tsv": "523925e6dbde0d7c48dcbe2d42303b01",
    f"{T}.z1000.verbose_mapping.tsv": "9ef46b00e351e589e4bb3b728bca36e0",
    "test.tigmint-ntLink-arks.fa": "c1f95d0be81facab4b7f4ff9615caf99",
}
NTLINK_ARKS_2 = {
    f"{N}.tsv": "a503dbcb68d6496f41c4fe192581a22d",
    f"{N}.z1000.n2.scaffold.dot": None,
    f"{N}.z1000.ntLink.scaffolds.renamed.fa": "19a351ba9b2d210e9c85e51a564914c7",
    f"{N}.z1000.{A}.dist.gv": "475a4a507efa20e07bed9c5a692096e8",
    f"{N}.z1000.{A}.tigpair_checkpoint.tsv": None,
    f"{N}.z1000.{A}_l4_a0.3.assembly_correspondence.tsv": "4eaa69384dda42a7c53b8c0af17ae5dd",
    f"{N}.z1000.{A}_l4_a0.3.gv": None,
    f"{N}.z1000.{A}_l4_a0.3.log": None,
    f"{N}.z1000.{A}_l4_a0.3.scaffolds": "2fc807fa033afb034e590f11285997d5",
    f"{N}.z1000.{A}_l4_a0.3.scaffolds.fa": "d91e1396cdb094ed451d9bf28639e76b",
    f"{N}.z1000.{A}_main.tsv": "dba0904ccd4241cac377c7ba55ffa168",
    f"{N}.z1000.{A}_original.gv": None,
    f"{N}.z1000.stitch.abyss-scaffold.fa": "e8337c0d74d0443a110b79ef1bed5b0f",
    f"{N}.z1000.stitch.path": "a72e013fa0d2b4f543ded5886e3b8f29",
    f"{N}.z1000.trimmed_scafs.agp": "f35e5b8993e00d18582418147bc42d16",
    f"{N}.z1000.trimmed_scafs.fa": "a817a1bbfd335f371a30426c34c4b2fe",
    f"{N}.z1000.trimmed_scafs.path": "07a0b5fd42d15f5108ac3b2453ffb5da",
    f"{N}.z1000.trimmed_scafs.tsv": "8d8fe50235b4954091e23d7b596ce87f",
    f"{N}.z1000.verbose_mapping.tsv": "eb36b8971fdcbf89405ed05a235d1bde",
}


async def check(dataset: int, genomesize: str, command: str, expected: dict[str, str | None]) -> None:
    assembly = await fixture(f"{DATA}test_scaffolds{dataset}.fa")
    reads = await fixture(f"{DATA}test_reads{dataset}.fa.gz")
    # The cache key ignores resources, so salt it to keep the 2-CPU run apart from a default one.
    task = longstitch_cmd.as_task().override(
        resources=flyte.Resources(cpu=2, memory="16Gi"), cache=flyte.Cache(behavior="auto", salt="cpu=2")
    )
    results = await task(
        assembly=assembly,
        reads=reads,
        command=command,
        span="",
        genomesize=genomesize,
        longmap="ont",
        prefix="test",
        args="",
    )
    result = await read_longstitch_results(results, command, "test")
    names = set(result.outputs)
    missing = sorted(set(expected) - names)
    assert not missing, f"longstitch {command} wrote no {missing}"
    extra = sorted(names - set(expected) - FINAL)
    assert not extra, f"longstitch {command} wrote unexpected {extra}"
    errors = []
    for name, md5 in expected.items():
        if md5 is not None:
            try:
                await assert_md5(result.outputs[name], md5, label=name)
            except AssertionError as e:
                errors.append(str(e))
    assert not errors, "; ".join(errors)
    await assert_nonempty(result.scaffolds, label=result.scaffolds.name)


@env.task
async def test_longstitch_run_1() -> None:
    await check(1, "482700", "run", TIGMINT_1)


@env.task
async def test_longstitch_tigmint_ntlink_arks_1() -> None:
    await check(1, "482700", "tigmint-ntLink-arks", TIGMINT_ARKS_1)


@env.task
async def test_longstitch_tigmint_ntlink_1() -> None:
    await check(1, "482700", "tigmint-ntLink", TIGMINT_1)


@env.task
async def test_longstitch_ntlink_arks_1() -> None:
    await check(1, "482700", "ntLink-arks", NTLINK_ARKS_1)


@env.task
async def test_longstitch_run_2() -> None:
    await check(2, "250000", "run", TIGMINT_2)


@env.task
async def test_longstitch_tigmint_ntlink_arks_2() -> None:
    await check(2, "250000", "tigmint-ntLink-arks", TIGMINT_ARKS_2)


@env.task
async def test_longstitch_tigmint_ntlink_2() -> None:
    await check(2, "250000", "tigmint-ntLink", TIGMINT_2)


@env.task
async def test_longstitch_ntlink_arks_2() -> None:
    await check(2, "250000", "ntLink-arks", NTLINK_ARKS_2)


tests = [
    test_longstitch_run_1,
    test_longstitch_tigmint_ntlink_arks_1,
    test_longstitch_tigmint_ntlink_1,
    test_longstitch_ntlink_arks_1,
    test_longstitch_run_2,
    test_longstitch_tigmint_ntlink_arks_2,
    test_longstitch_tigmint_ntlink_2,
    test_longstitch_ntlink_arks_2,
]
