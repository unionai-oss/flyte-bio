"""Flyte-native test framework for flyte-bio.

Tests are themselves ``@env.task`` functions; the driver fans them out
via ``asyncio.gather``. Each test ends up as a Flyte action — own pod,
own logs, cached fixtures — so a suite of hundreds scales horizontally
instead of contending for one machine.

Tests are run with the native ``flyte run`` CLI — there is no bespoke
entry point. ``flyte run`` reads the standard Flyte config (``--config``
flag, ``FLYTE_CONFIG``, or ``~/.flyte/config.yaml``) to pick the cluster,
project, domain, and auth, then imports this package and submits the
driver task::

    # against whichever cluster your default config points at
    flyte run tests/run_all.py run_all

    # against a specific cluster
    flyte --config ~/.flyte/prod.yaml run tests/run_all.py run_all

    # locally, without a cluster
    flyte run --local tests/run_all.py run_all

The same code runs everywhere; the only thing that changes is the config
``flyte run`` resolves. ``flyte run`` calls ``flyte.init`` *before*
importing this module, so the image registry below auto-resolves to the
cluster you're targeting (see ``IMAGE_REGISTRY``). Run from the project
root so the ``tests`` package is importable (``flyte run`` puts the cwd
on ``sys.path``).

Conventions
-----------

A test is an async function decorated with ``@env.task`` that takes no
arguments and returns ``None``. It awaits the tool wrappers being tested
and raises ``AssertionError`` on failure::

    from tests.framework import env, fixture, assert_md5
    from flyte_bio.modules.bedtools import bedtools_intersect

    @env.task
    async def test_intersect_bed_bed() -> None:
        a = await fixture("genomics/sarscov2/.../a.bed")
        b = await fixture("genomics/sarscov2/.../b.bed")
        out = await bedtools_intersect(a=a, b=[b])
        await assert_md5(out, "<md5 from upstream snapshot>")

A driver — typically ``tests/run_all.py`` or a per-module ``run.py`` —
collects tests into a list and hands them to :func:`gather_tests`, which
runs them in parallel, prints a summary, and raises if any failed.

Fixtures
--------

:func:`fixture` downloads a test-data file (resolved against
:data:`TEST_DATA_BASE`) and stages it as a Flyte :class:`File`. It's an
ordinary ``@env.task(cache="auto")`` — Flyte handles the caching on the
task's input path, so the same fixture is fetched at most once per
deployment. Nothing is cached by hand.
"""



import asyncio
import hashlib
import os
import tempfile
from collections.abc import Awaitable, Callable
from pathlib import Path
from urllib.parse import urlsplit

import flyte
import flyte.storage as storage
import requests
from flyte.io import File

from flyte_bio.modules import env as modules_env

# The package under test must be *installed* into the test image, not just
# code-bundled. flyte's default bundle copies source files preserving their
# path relative to the project root, so this `src/`-layout package would land
# at `/root/src/flyte_bio/...` — but only `/root` is on `sys.path`, so
# `import flyte_bio` would fail in the pod. `with_uv_project` with
# `install_project` builds the wheel and `uv`-installs `flyte_bio` + its
# locked deps (incl. flyte) into the image, so the import resolves regardless
# of layout. The lockfile pins flyte, so we don't also `install_flyte`.
# (Test code under `tests/` is *not* installed — it rides in via the code
# bundle, which works because it sits at the project root, not `src/`.)
#
# The push registry is left to flyte: it reads the active `flyte.init` config
# (set up by `flyte run` before this module imports) and resolves to match
# the target cluster.
PYPROJECT = Path(__file__).resolve().parent.parent / "pyproject.toml"

TEST_IMAGE = flyte.Image.from_debian_base(
    install_flyte=False
).with_uv_project(PYPROJECT, project_install_mode="install_project")

env = flyte.TaskEnvironment(
    name="flyte_bio_tests",
    image=TEST_IMAGE,
    depends_on=[modules_env],
)

# Base URL that relative fixture paths are resolved against.
TEST_DATA_BASE = "https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/"


def download_http(url: str) -> str:
    """Stream an ``http(s)`` URL to a local temp file and return its path.

    Deliberately uses ``requests`` rather than :func:`flyte.storage.get`:
    that routes ``https`` through fsspec's ``HTTPFileSystem``, which binds
    its aiohttp ``ClientSession`` to fsspec's own background event loop.
    flyte then ``await``\\ s the filesystem coroutine on the task's loop
    instead, and aiohttp's timeout machinery raises "Timeout context
    manager should be used inside a task" across the loop mismatch. A plain
    synchronous ``requests`` download has no event loop to mismatch; call
    it from a worker thread via :func:`asyncio.to_thread`.
    """
    suffix = Path(urlsplit(url).path).suffix
    fd, path = tempfile.mkstemp(suffix=suffix)
    with os.fdopen(fd, "wb") as out, requests.get(url, stream=True, timeout=60) as resp:
        resp.raise_for_status()
        for chunk in resp.iter_content(1 << 20):
            out.write(chunk)
    return path


@env.task(cache="auto")
async def fixture(rel_path: str) -> File:
    """Download a test-data fixture and stage it as a Flyte :class:`File`.

    ``rel_path`` is resolved against :data:`TEST_DATA_BASE`
    (e.g. ``genomics/sarscov2/illumina/bam/test.paired_end.bam``).

    Dispatches on the URL scheme: ``http(s)`` is streamed directly with
    :func:`download_http` (see there for why we bypass fsspec), while other
    schemes go through :func:`flyte.storage.get`, which is obstore-native
    and loop-safe for ``s3://``. So moving the fixtures into S3 later stays
    a one-line :data:`TEST_DATA_BASE` edit. ``cache="auto"`` does the real
    caching, so the same path is fetched at most once per deployment.
    """
    url = TEST_DATA_BASE + rel_path.lstrip("/")
    if url.startswith(("http://", "https://")):
        local = await asyncio.to_thread(download_http, url)
    else:
        local = await storage.get(url)
    return await File.from_local(local)


async def assert_nonempty(file: File, *, label: str = "") -> None:
    """Assert ``file`` exists and has nonzero length (a "run-to-green" check).

    Used for outputs whose bytes aren't reproducible enough to md5 — e.g.
    tool reports that embed the command line, version banner, or a
    timestamp, and compressed outputs whose block layout depends on the
    thread count. We still want to prove the wrapper ran and produced the
    file, just not pin its exact content.
    """
    size = 0
    async with file.open("rb") as fh:
        while True:
            chunk = await fh.read(1 << 20)
            if not chunk:
                break
            size += len(chunk)
            if size:
                break
    if size == 0:
        prefix = f"{label}: " if label else ""
        raise AssertionError(f"{prefix}expected a nonempty file, got 0 bytes")


async def assert_dir_nonempty(directory, *, label: str = "") -> None:
    """Assert ``directory`` (a :class:`Dir`) holds at least one file.

    A run-to-green check for outputs like a built index, whose individual
    files aren't reproducible enough to md5 but whose presence proves the
    tool ran.
    """
    files = await directory.list_files()
    if not files:
        prefix = f"{label}: " if label else ""
        raise AssertionError(f"{prefix}expected a nonempty directory, got 0 files")


async def assert_md5(file: File, expected: str, *, label: str = "") -> None:
    """Assert that ``file``'s md5 matches ``expected``.

    The expected MD5 typically comes from an upstream test snapshot,
    where each entry is shaped ``"<filename>:md5,<hex>"``. Pass just the
    hex.

    Streams the file so it stays memory-bounded even for large BAMs.
    """
    h = hashlib.md5()
    async with file.open("rb") as fh:
        while True:
            chunk = await fh.read(1 << 20)
            if not chunk:
                break
            h.update(bytes(chunk))
    actual = h.hexdigest()
    if actual != expected:
        prefix = f"{label}: " if label else ""
        raise AssertionError(
            f"{prefix}md5 mismatch (expected {expected}, got {actual})"
        )


# A "test" is anything callable with no args returning Awaitable[None] —
# typically a `@env.task` function. Tests carry a `.name` (the qualified
# task name); plain async functions fall back to `__name__`.
Test = Callable[[], Awaitable[None]]


def test_name(t: Test) -> str:
    return getattr(t, "name", None) or getattr(t, "__name__", repr(t))


async def gather_tests(tests: list[Test]) -> str:
    """Run ``tests`` concurrently, print a summary, raise if any failed.

    Designed to be the body of a top-level driver task::

        @env.task
        async def run_all() -> str:
            return await gather_tests([
                test_intersect_bed_bed,
                test_bedtools_sort,
                ...,
            ])

    Exceptions raised inside individual tests are caught (each test's
    failure is visible per-task in the Flyte UI). After all have settled,
    a summary is printed and an aggregated ``AssertionError`` is raised
    if any failed.
    """
    results = await asyncio.gather(
        *(t() for t in tests),
        return_exceptions=True,
    )
    passed: list[str] = []
    failed: list[tuple[str, BaseException]] = []
    for t, r in zip(tests, results):
        name = test_name(t)
        if isinstance(r, BaseException):
            failed.append((name, r))
        else:
            passed.append(name)

    bar = "=" * 72
    lines = [
        "",
        bar,
        f"{len(passed)} passed, {len(failed)} failed (of {len(tests)})",
        bar,
    ]
    for name in passed:
        lines.append(f"  PASS  {name}")
    for name, err in failed:
        lines.append(f"  FAIL  {name}")
        lines.append(f"        {type(err).__name__}: {err}")
    summary = "\n".join(lines)
    print(summary)

    if failed:
        raise AssertionError(f"{len(failed)}/{len(tests)} tests failed")
    return summary


__all__ = ["assert_dir_nonempty", "assert_md5", "assert_nonempty", "env", "fixture", "gather_tests"]
