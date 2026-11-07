"""Tests for flyte_bio.modules.rsemmergecounts.

The case is upstream's custom/rsemmergecounts test (modules commit 9a1cc9f,
sed 4.7 image), with its snapshot md5s for all six tables: three yeast
samples' RSEM results, unpacked from upstream's tarball.
"""

from flyte_bio.modules.rsemmergecounts import rsem_merge_counts
from tests.framework import assert_md5, env, fixture_dir

MD5S = {
    "counts_gene": "0a100470da54fc96e19ae70f0afdfbf0",
    "tpm_gene": "338689117908cdd08fe10d6fd4b38f17",
    "counts_transcript": "d1c7e727f1d1da94def08a0d5e735fb1",
    "tpm_transcript": "bf07c9cba378995b3626ff60c0a2d6b9",
    "genes_long": "93c7d2e9d1f2d49bc367939f4f15ff1d",
    "isoforms_long": "17c5217b25bc42db6f9fce802d1dbe73",
}


@env.task
async def test_rsemmergecounts_saccharomyces_cerevisiae() -> None:
    # upstream case: custom/rsemmergecounts "saccharomyces_cerevisiae"
    unpacked = await fixture_dir("genomics/eukaryotes/saccharomyces_cerevisiae/rsem_results.tar.gz")
    files = [f async for f in unpacked.walk()]
    genes = [f for f in files if f.path.endswith(".genes.results")]
    isoforms = [f for f in files if f.path.endswith(".isoforms.results")]
    assert len(genes) == 3 and len(isoforms) == 3, f"rsem_results.tar.gz holds {[f.path for f in files]}"
    merged = await rsem_merge_counts(genes, isoforms, prefix="rsem_merged")
    for field, md5 in MD5S.items():
        await assert_md5(getattr(merged, field), md5, label=field)


tests = [test_rsemmergecounts_saccharomyces_cerevisiae]
