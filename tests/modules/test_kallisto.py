"""Tests for flyte_bio.modules.kallisto.

The case is upstream's kallisto/index test (modules commit 9a1cc9f, kallisto
0.52.0). Upstream marks the index unstable and snapshots only its name, so the
test checks the name and that the index isn't empty.
"""

from flyte_bio.modules.kallisto import kallisto_index
from tests.framework import assert_nonempty, env, fixture


@env.task
async def test_kallisto_index_transcriptome() -> None:
    # upstream case: kallisto/index "sarscov2 transcriptome.fasta"
    index = await kallisto_index(await fixture("genomics/sarscov2/genome/transcriptome.fasta"))
    assert index.name == "kallisto", index.name
    await assert_nonempty(index, label="kallisto index")


tests = [test_kallisto_index_transcriptome]
