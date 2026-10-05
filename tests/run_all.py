"""Top-level test driver — runs every test module's suite as one action.

Run it with the native ``flyte run`` CLI from the project root. The
cluster, project, domain, and auth all come from your Flyte config
(``--config`` flag, ``FLYTE_CONFIG``, or ``~/.flyte/config.yaml``)::

    # whichever cluster your default config points at (e.g. a devbox)
    flyte run tests/run_all.py run_all

    # a specific cluster
    flyte --config ~/.flyte/prod.yaml run tests/run_all.py run_all

    # locally, without a cluster
    flyte run --local tests/run_all.py run_all

Each test in :data:`ALL_TESTS` becomes its own action, so hundreds of
tests scale horizontally on whichever cluster the config points at.

To iterate on a subset, write a sibling driver — e.g.
``tests/modules/run.py`` that imports only the module tests it wants and
passes them to :func:`gather_tests` — and run it the same way:
``flyte run tests/modules/run.py <task>``.
"""



from tests.framework import env, gather_tests
from tests.modules import (
    test_bedtools,
    test_cat,
    test_catadditionalfasta,
    test_gffread,
    test_gtf2bed,
    test_gtffilter,
    test_salmon,
    test_samtools,
    test_star,
)
from tests.pipelines import test_rnaseq

ALL_TESTS = [
    *test_bedtools.tests,
    *test_cat.tests,
    *test_catadditionalfasta.tests,
    *test_gffread.tests,
    *test_gtf2bed.tests,
    *test_gtffilter.tests,
    *test_salmon.tests,
    *test_samtools.tests,
    *test_star.tests,
    *test_rnaseq.tests,
]


@env.task
async def run_all() -> str:
    return await gather_tests(ALL_TESTS)
