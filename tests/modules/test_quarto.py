"""Tests for flyte_bio.modules.quarto.

The quarto/notebook cases and expected md5s are upstream's (modules commit
9a1cc9f, Quarto 1.7.31). Upstream pins them to 1 CPU (``cpus`` is written into
``params.yml``), so these run the task at 1 CPU; ``versions.yml`` is keyed by
the process name, ``QUARTO_NOTEBOOK``, as upstream. The quartonotebook cases
are that module's last working version's tests (commit 6d46786; it is
deprecated at 9a1cc9f in favour of quarto/notebook): the same notebooks at
upstream's default 2 CPUs, with no ``versions.csv`` required, and the Jupyter
notebook rendered without ``--execute`` (so it writes no artifact). As
upstream, the HTML is checked for the notebook's output, not by md5.
"""

import flyte
from flyte.io import File

from flyte_bio.modules.quarto import (
    QuartoNotebookResult,
    notebook_parameters,
    notebook_parameters_file,
    quarto_notebook,
    quarto_notebook_cmd,
    quarto_notebook_result,
)
from tests.framework import assert_md5, env, fixture

NOTEBOOKS = "generic/notebooks/"
QMD_R = NOTEBOOKS + "quarto/quarto_r.qmd"
QMD_PYTHON = NOTEBOOKS + "quarto/quarto_python.qmd"
RMD = NOTEBOOKS + "rmarkdown/rmarkdown_notebook.Rmd"
IPYNB = NOTEBOOKS + "jupyter/ipython_notebook.ipynb"
PARAMETERS = {"input_filename": "hello.txt", "n_iter": 12}

R_ARTIFACT_MD5 = "b10a8db164e0754105b7a99be72e3fe5"  # "Hello World", from R's cat()
PYTHON_ARTIFACT_MD5 = "8ddd8be4b179a529afa5f2ffae4b9858"  # "Hello World!\n"
R_VERSIONS_MD5 = "ffbc211e9dcc1252c411bdd2f072b92e"  # R 4.4.3
PYTHON_VERSIONS_MD5 = "68fd39542815eb5537cdb4b234e2d30d"  # python 3.12.11 (csv rows end in \r\n)
PARAMS_MD5 = {  # by CPUs and whether parameters were given
    (1, False): "6d2ec3cf3e2f15958ab8ab440d651f8d",
    (1, True): "4fd353560b2c3366515af2cc4d49957d",
    (2, False): "8b9e438c0eb4850e20035a1222cd151e",
    (2, True): "0dbcd58276179f0700e144e45d8d3c76",
}


async def text(file: File) -> str:
    async with file.open("rb") as fh:
        return bytes(await fh.read()).decode()


async def render_at_one_cpu(path: str, parameters: dict | None, args: str = "") -> QuartoNotebookResult:
    """quarto/notebook at 1 CPU, as upstream's tests pin it."""
    notebook = await fixture(path)
    input_files = [await fixture("generic/txt/hello.txt")] if parameters else []
    task = quarto_notebook_cmd.as_task().override(
        resources=flyte.Resources(cpu=1, memory="8Gi"), cache=flyte.Cache(behavior="auto", salt="cpu=1")
    )
    results = await task(
        notebook=notebook,
        params_yaml=await notebook_parameters_file(notebook_parameters({"id": "test"}, parameters)),
        input_files=input_files,
        extensions=None,
        artifact_dir="artifacts",
        prefix="test",
        notebook_versions=True,
        versions_key="QUARTO_NOTEBOOK",
        args=args,
    )
    return await quarto_notebook_result(results, "test", path.rsplit("/", 1)[-1])


async def render_legacy(path: str, parameters: dict | None) -> QuartoNotebookResult:
    """quartonotebook: no versions.csv required, at the default 2 CPUs."""
    input_files = [await fixture("generic/txt/hello.txt")] if parameters else []
    return await quarto_notebook(
        await fixture(path),
        parameters=parameters,
        input_files=input_files,
        prefix="test",
        notebook_versions=False,
    )


async def check(
    result: QuartoNotebookResult,
    cpus: int,
    parametrized: bool,
    artifact_md5: str | None,
    versions_md5: str | None,
    html: str = "Hello World 1",
) -> None:
    assert html in await text(result.html), f"test.html lacks {html!r}"
    await assert_md5(result.params_yaml, PARAMS_MD5[cpus, parametrized], label="params.yml")
    if parametrized:
        assert "meta:" in await text(result.params_yaml), "params.yml lacks meta:"
    names = [f.path.rsplit("/", 1)[-1] for f in result.artifacts]
    assert names == (["artifact.txt"] if artifact_md5 else []), f"artifacts are {names}"
    if artifact_md5:
        await assert_md5(result.artifacts[0], artifact_md5, label="artifact.txt")
    if versions_md5:
        assert result.versions is not None, "quarto wrote no versions.yml"
        await assert_md5(result.versions, versions_md5, label="versions.yml")
    else:
        assert result.versions is None, "quarto wrote versions.yml without notebook_versions"


@env.task
async def test_quarto_notebook_qmd_r() -> None:
    # upstream case: quarto/notebook "test notebook - [qmd:r]"
    await check(await render_at_one_cpu(QMD_R, None), 1, False, R_ARTIFACT_MD5, R_VERSIONS_MD5)


@env.task
async def test_quarto_notebook_qmd_python() -> None:
    # upstream case: quarto/notebook "test notebook - [qmd:python]"
    await check(await render_at_one_cpu(QMD_PYTHON, None), 1, False, PYTHON_ARTIFACT_MD5, PYTHON_VERSIONS_MD5)


@env.task
async def test_quarto_notebook_parametrized_qmd_r() -> None:
    # upstream case: quarto/notebook "test notebook - parametrized - [qmd:r]"
    await check(await render_at_one_cpu(QMD_R, PARAMETERS), 1, True, R_ARTIFACT_MD5, R_VERSIONS_MD5)


@env.task
async def test_quarto_notebook_parametrized_qmd_python() -> None:
    # upstream case: quarto/notebook "test notebook - parametrized - [qmd:python]"
    result = await render_at_one_cpu(QMD_PYTHON, PARAMETERS)
    await check(result, 1, True, PYTHON_ARTIFACT_MD5, PYTHON_VERSIONS_MD5)


@env.task
async def test_quarto_notebook_parametrized_rmd() -> None:
    # upstream case: quarto/notebook "test notebook - parametrized - [rmd]"
    await check(await render_at_one_cpu(RMD, PARAMETERS), 1, True, R_ARTIFACT_MD5, R_VERSIONS_MD5)


@env.task
async def test_quarto_notebook_parametrized_ipynb() -> None:
    # upstream case: quarto/notebook "test notebook - parametrized - [ipynb]" (with --execute)
    result = await render_at_one_cpu(IPYNB, PARAMETERS, args="--execute")
    await check(result, 1, True, PYTHON_ARTIFACT_MD5, PYTHON_VERSIONS_MD5, html="Hello World")


@env.task
async def test_quarto_legacy_qmd_r() -> None:
    # upstream case: quartonotebook "test notebook - [qmd:r]"
    await check(await render_legacy(QMD_R, None), 2, False, R_ARTIFACT_MD5, None)


@env.task
async def test_quarto_legacy_qmd_python() -> None:
    # upstream case: quartonotebook "test notebook - [qmd:python]"
    await check(await render_legacy(QMD_PYTHON, None), 2, False, PYTHON_ARTIFACT_MD5, None)


@env.task
async def test_quarto_legacy_parametrized_qmd_r() -> None:
    # upstream case: quartonotebook "test notebook - parametrized - [qmd:r]"
    await check(await render_legacy(QMD_R, PARAMETERS), 2, True, R_ARTIFACT_MD5, None)


@env.task
async def test_quarto_legacy_parametrized_qmd_python() -> None:
    # upstream case: quartonotebook "test notebook - parametrized - [qmd:python]"
    await check(await render_legacy(QMD_PYTHON, PARAMETERS), 2, True, PYTHON_ARTIFACT_MD5, None)


@env.task
async def test_quarto_legacy_parametrized_rmd() -> None:
    # upstream case: quartonotebook "test notebook - parametrized - [rmd]"
    await check(await render_legacy(RMD, PARAMETERS), 2, True, R_ARTIFACT_MD5, None)


@env.task
async def test_quarto_legacy_parametrized_ipynb() -> None:
    # upstream case: quartonotebook "test notebook - parametrized - [ipynb]" (not executed: no artifact)
    await check(await render_legacy(IPYNB, PARAMETERS), 2, True, None, None, html="Hello World")


tests = [
    test_quarto_notebook_qmd_r,
    test_quarto_notebook_qmd_python,
    test_quarto_notebook_parametrized_qmd_r,
    test_quarto_notebook_parametrized_qmd_python,
    test_quarto_notebook_parametrized_rmd,
    test_quarto_notebook_parametrized_ipynb,
    test_quarto_legacy_qmd_r,
    test_quarto_legacy_qmd_python,
    test_quarto_legacy_parametrized_qmd_r,
    test_quarto_legacy_parametrized_qmd_python,
    test_quarto_legacy_parametrized_rmd,
    test_quarto_legacy_parametrized_ipynb,
]
