"""Tests for flyte_bio.modules.sigprofiler.

The cases are upstream's sigprofiler tests (modules commit 9a1cc9f,
SigProfilerExtractor 1.2.2, SigProfilerMatrixGenerator 1.3.5), with their
``ext.args``:

- a single-sample run against GRCh38 with upstream's fixed seed: as upstream,
  every file under ``results/`` is compared by name, and by md5 except for the
  PDFs, ``JOB_METADATA.txt`` and the ``*_seqinfo.txt`` files (names only) and
  the timestamped ``*.err`` / ``*.out`` logs (ignored). Upstream's snapshot
  nests the files by directory without naming the directories, so the check
  compares the multiset of ``name:md5`` entries. Upstream's test downloads
  GRCh38 in the run; since 2026-09-25 that download is a newer build the
  image's SigProfilerMatrixGenerator rejects by checksum, so the test installs
  the reference upstream's snapshot was made with (now ``GRCh38_Legacy``) as
  GRCh38 first. The template's install step then finds it installed.
- two runs that must fail on bad ``args``. Upstream checks the error report
  for a message, which Nextflow's report satisfies by quoting the template;
  here the task's failure is checked. Upstream's second case
  (``download_genome_sigprofiler:False`` without ``genome_installed_path``)
  also passes ``make_decomposition_plots:maybe``, which the template rejects
  first.
"""

import asyncio
import hashlib
import tarfile
import tempfile
from collections import Counter
from pathlib import Path

from flyte.io import File

from flyte_bio.modules.sigprofiler import sigprofiler, sigprofiler_install_genome
from tests.framework import assert_md5, env, fixture

TABLE = "delete_me/sigprofiler/test_GRCh38/HCC1395_WXS_HCC1395_joint_table.tsv"
ARGS = (
    "[download_genome_sigprofiler:True, input_type:matrix, context_type:'96,DINUC,ID', minimum_signatures:1, "
    "maximum_signatures:2, nmf_replicates:10, min_nmf_iterations:10, max_nmf_iterations:100, nmf_test_conv:100, "
    "make_decomposition_plots:False ]"
)
BAD_ARGS = (
    "[{download}, volume:./, input_type:matrix, context_type:'96,DINUC,ID', minimum_signatures:1, "
    "maximum_signatures:2, nmf_replicates:10, min_nmf_iterations:10, max_nmf_iterations:100, nmf_test_conv:100, "
    "make_decomposition_plots:maybe ]"
)
# Upstream's seed channel, collected into seed.txt with a newline after each line.
SEEDS = "\tSeed\n0\t315836613344651455423314970566641622150\n"
VERSIONS_MD5 = "480098d2a5d6c34f7d3eb468201d123c"

# Upstream's snapshot of results/, flattened and sorted: "name:md5,<hex>", or the bare name of an unstable file.
EXPECTED_RESULTS = (
    "10_seqinfo.txt",
    "10_seqinfo.txt",
    "10_seqinfo.txt",
    "10_seqinfo.txt",
    "11_seqinfo.txt",
    "11_seqinfo.txt",
    "11_seqinfo.txt",
    "11_seqinfo.txt",
    "12_seqinfo.txt",
    "12_seqinfo.txt",
    "12_seqinfo.txt",
    "12_seqinfo.txt",
    "13_seqinfo.txt",
    "13_seqinfo.txt",
    "13_seqinfo.txt",
    "13_seqinfo.txt",
    "14_seqinfo.txt",
    "14_seqinfo.txt",
    "14_seqinfo.txt",
    "14_seqinfo.txt",
    "15_seqinfo.txt",
    "15_seqinfo.txt",
    "15_seqinfo.txt",
    "15_seqinfo.txt",
    "16_seqinfo.txt",
    "16_seqinfo.txt",
    "16_seqinfo.txt",
    "16_seqinfo.txt",
    "17_seqinfo.txt",
    "17_seqinfo.txt",
    "17_seqinfo.txt",
    "17_seqinfo.txt",
    "18_seqinfo.txt",
    "18_seqinfo.txt",
    "18_seqinfo.txt",
    "18_seqinfo.txt",
    "19_seqinfo.txt",
    "19_seqinfo.txt",
    "19_seqinfo.txt",
    "19_seqinfo.txt",
    "1_seqinfo.txt",
    "1_seqinfo.txt",
    "1_seqinfo.txt",
    "1_seqinfo.txt",
    "20_seqinfo.txt",
    "20_seqinfo.txt",
    "20_seqinfo.txt",
    "20_seqinfo.txt",
    "21_seqinfo.txt",
    "21_seqinfo.txt",
    "21_seqinfo.txt",
    "21_seqinfo.txt",
    "22_seqinfo.txt",
    "22_seqinfo.txt",
    "22_seqinfo.txt",
    "22_seqinfo.txt",
    "2_seqinfo.txt",
    "2_seqinfo.txt",
    "2_seqinfo.txt",
    "2_seqinfo.txt",
    "3_seqinfo.txt",
    "3_seqinfo.txt",
    "3_seqinfo.txt",
    "3_seqinfo.txt",
    "4_seqinfo.txt",
    "4_seqinfo.txt",
    "4_seqinfo.txt",
    "4_seqinfo.txt",
    "5_seqinfo.txt",
    "5_seqinfo.txt",
    "5_seqinfo.txt",
    "5_seqinfo.txt",
    "6_seqinfo.txt",
    "6_seqinfo.txt",
    "6_seqinfo.txt",
    "6_seqinfo.txt",
    "7_seqinfo.txt",
    "7_seqinfo.txt",
    "7_seqinfo.txt",
    "7_seqinfo.txt",
    "8_seqinfo.txt",
    "8_seqinfo.txt",
    "8_seqinfo.txt",
    "8_seqinfo.txt",
    "9_seqinfo.txt",
    "9_seqinfo.txt",
    "9_seqinfo.txt",
    "9_seqinfo.txt",
    "All_solutions_stat.csv:md5,b2e621090867708c4b3b113c02017982",
    "All_solutions_stat.csv:md5,c9e60f5b13219244a5c295c5076b3f48",
    "All_solutions_stat.csv:md5,d26a62df9f7668246b1c8b26a4222af4",
    "COSMIC_DBS78_Activities.txt:md5,c20aab9a5b0e72e3f451c6b81e2326cc",
    "COSMIC_DBS78_Activity_Plots.pdf",
    "COSMIC_DBS78_Samples_Stats.txt:md5,2ca5d1d54782a18429566e1b688eb644",
    "COSMIC_DBS78_Signature_Assignment_log.txt:md5,2204c302eb0a691e175c582557847bb8",
    "COSMIC_DBS78_Signatures.txt:md5,fe1675e6e6f26f06550cb2cecfd24828",
    "COSMIC_DBS78_TMB_plot.pdf",
    "COSMIC_ID83_Activities.txt:md5,06115e584d4c6552ccc7363b112327e7",
    "COSMIC_ID83_Activity_Plots.pdf",
    "COSMIC_ID83_Samples_Stats.txt:md5,81a95dc25aaaa9a245b85fab2a3b1878",
    "COSMIC_ID83_Signature_Assignment_log.txt:md5,afa42172d82a4eb61f5b0650360d48af",
    "COSMIC_ID83_Signatures.txt:md5,351b21dd0837c0dbb575bc8ad15ed208",
    "COSMIC_ID83_TMB_plot.pdf",
    "COSMIC_SBS96_Activities.txt:md5,1e7c6e4fcd6db9aa4e7f6f4dcdfeb165",
    "COSMIC_SBS96_Activity_Plots.pdf",
    "COSMIC_SBS96_Samples_Stats.txt:md5,568ed4a39c1b47505c8a7cee6810e6df",
    "COSMIC_SBS96_Signature_Assignment_log.txt:md5,5d79aab4ba27961177859e026e8fdabb",
    "COSMIC_SBS96_Signatures.txt:md5,833454d46f359bbd5edd21f33204eee3",
    "COSMIC_SBS96_TMB_plot.pdf",
    "Cosmic_DBS78_Decomposition_Log.txt:md5,949615a327a4231c94d46b4c881c02ff",
    "Cosmic_ID83_Decomposition_Log.txt:md5,d900442dd5497ee8d6506d882ab849ff",
    "Cosmic_SBS96_Decomposition_Log.txt:md5,c160005e4e25d9709961569c258e39dd",
    "DBS78_De-Novo_Activities_refit.txt:md5,f8ff32f6b7d13150d55c03288e524a84",
    "DBS78_De-Novo_Activity_Plots_refit.pdf",
    "DBS78_De-Novo_Samples_Stats_refit.txt:md5,d2a44d214a6ce187001c4ce5707937e9",
    "DBS78_De-Novo_Signature_Assignment_log.txt:md5,0f06e52a7c4a04a6a3882ca77c246b41",
    "DBS78_De-Novo_Signatures.txt:md5,092f589e7e1c993357a534a747fbadc5",
    "DBS78_De-Novo_Signatures_SEM_Error.txt:md5,cc904e8d4cb9999b2aadd8362e00acd8",
    "DBS78_De-Novo_TMB_plot_refit.pdf",
    "DBS78_De-Novo_refit_Signature_Assignment_log.txt:md5,4185dbec591763f2e1d4233dbf328efb",
    "DBS78_S1_NMF_Activities.txt:md5,f8ff32f6b7d13150d55c03288e524a84",
    "DBS78_S1_NMF_Activities_SEM_Error.txt:md5,c7ce4ba94cc043709d94858cd5808ec5",
    "DBS78_S1_NMF_Activity_Plots.pdf",
    "DBS78_S1_NMF_Convergence_Information.txt:md5,e85c5ed874d700c99fad20cb1ce54f9b",
    "DBS78_S1_Samples_stats.txt:md5,d2a44d214a6ce187001c4ce5707937e9",
    "DBS78_S1_Signatures.txt:md5,2ddaca19400389a6b5038fd39a92d185",
    "DBS78_S1_Signatures_SEM_Error.txt:md5,cc904e8d4cb9999b2aadd8362e00acd8",
    "DBS78_S1_Signatures_stats.txt:md5,c781a450d458dba2c7268e6a83855aec",
    "DBS78_S1_TMB_NMF_plot.pdf",
    "DBS78_selection_plot.pdf",
    "DBS_78_plots_COSMIC_DBS78.pdf",
    "DBS_78_plots_DBS78_De-Novo.pdf",
    "DBS_78_plots_S1.pdf",
    "De_Novo_MutationType_Probabilities.txt:md5,63f2d244a39e503ebc78b759fb0d916e",
    "De_Novo_MutationType_Probabilities.txt:md5,9e2c10da35e3e4a3552d42ace8d127dd",
    "De_Novo_MutationType_Probabilities.txt:md5,f6bf0e921fec88a5530260b71cedba29",
    "De_Novo_MutationType_Probabilities_refit.txt:md5,93e74401506989182928195b7268d999",
    "De_Novo_MutationType_Probabilities_refit.txt:md5,d2099d19ac1034d624443ae88b43c0be",
    "De_Novo_MutationType_Probabilities_refit.txt:md5,fbb27313f0c8704b40dcff7e509beb5a",
    "De_Novo_map_to_COSMIC_DBS78.csv:md5,437090cdf67c7772f82b8c8f1a2a5b85",
    "De_Novo_map_to_COSMIC_ID83.csv:md5,14e9613461e18a51eb5b170d9e361e45",
    "De_Novo_map_to_COSMIC_SBS96.csv:md5,440a3e9b78505dfa14cbb75a64113a6b",
    "ID83_De-Novo_Activities_refit.txt:md5,f5c960d38792bae21b357f79328d443d",
    "ID83_De-Novo_Activity_Plots_refit.pdf",
    "ID83_De-Novo_Samples_Stats_refit.txt:md5,a14ff9e2c1d060fdf80b0a429e7072a3",
    "ID83_De-Novo_Signature_Assignment_log.txt:md5,4dfcb629c95d3a8ad9a6a9cbaaab42ac",
    "ID83_De-Novo_Signatures.txt:md5,d34d0ead07eed437de64e692b922a723",
    "ID83_De-Novo_Signatures_SEM_Error.txt:md5,55db03cf759e983cfb6b0d3504540fa3",
    "ID83_De-Novo_TMB_plot_refit.pdf",
    "ID83_De-Novo_refit_Signature_Assignment_log.txt:md5,4185dbec591763f2e1d4233dbf328efb",
    "ID83_S1_NMF_Activities.txt:md5,f5c960d38792bae21b357f79328d443d",
    "ID83_S1_NMF_Activities_SEM_Error.txt:md5,aef7288ee354681c638000c67c18f266",
    "ID83_S1_NMF_Activity_Plots.pdf",
    "ID83_S1_NMF_Convergence_Information.txt:md5,e85c5ed874d700c99fad20cb1ce54f9b",
    "ID83_S1_Samples_stats.txt:md5,a14ff9e2c1d060fdf80b0a429e7072a3",
    "ID83_S1_Signatures.txt:md5,83712a37739e06671f199e7caffa4603",
    "ID83_S1_Signatures_SEM_Error.txt:md5,55db03cf759e983cfb6b0d3504540fa3",
    "ID83_S1_Signatures_stats.txt:md5,907cceb66b0d028e406247a1b9ff87ac",
    "ID83_S1_TMB_NMF_plot.pdf",
    "ID83_selection_plot.pdf",
    "ID_83_plots_COSMIC_ID83.pdf",
    "ID_83_plots_ID83_De-Novo.pdf",
    "ID_83_plots_S1.pdf",
    "JOB_METADATA.txt",
    "JOB_METADATA.txt",
    "JOB_METADATA.txt",
    "MT_seqinfo.txt",
    "MT_seqinfo.txt",
    "MT_seqinfo.txt",
    "MT_seqinfo.txt",
    "SBS96_De-Novo_Activities_refit.txt:md5,37fb52947411c1ce69b62c656cb539db",
    "SBS96_De-Novo_Activity_Plots_refit.pdf",
    "SBS96_De-Novo_Samples_Stats_refit.txt:md5,7a3f12db508724491597e0b279c3d9d2",
    "SBS96_De-Novo_Signature_Assignment_log.txt:md5,5b2aa2df5c7b95985f495aba88b5f960",
    "SBS96_De-Novo_Signatures.txt:md5,f81b6baf37222e0ba8abcad45b5ae883",
    "SBS96_De-Novo_Signatures_SEM_Error.txt:md5,6a726117def78216d0f2ec4d043afc5a",
    "SBS96_De-Novo_TMB_plot_refit.pdf",
    "SBS96_De-Novo_refit_Signature_Assignment_log.txt:md5,4185dbec591763f2e1d4233dbf328efb",
    "SBS96_S1_NMF_Activities.txt:md5,37fb52947411c1ce69b62c656cb539db",
    "SBS96_S1_NMF_Activities_SEM_Error.txt:md5,8214ca37159418e18f6fb852b6f41a72",
    "SBS96_S1_NMF_Activity_Plots.pdf",
    "SBS96_S1_NMF_Convergence_Information.txt:md5,e85c5ed874d700c99fad20cb1ce54f9b",
    "SBS96_S1_Samples_stats.txt:md5,7a3f12db508724491597e0b279c3d9d2",
    "SBS96_S1_Signatures.txt:md5,5aa9d63d406dc474fd774d8447202573",
    "SBS96_S1_Signatures_SEM_Error.txt:md5,6a726117def78216d0f2ec4d043afc5a",
    "SBS96_S1_Signatures_stats.txt:md5,8dc7655067609a51e7a374a69a63313e",
    "SBS96_S1_TMB_NMF_plot.pdf",
    "SBS96_selection_plot.pdf",
    "SBS_96_plots_COSMIC_SBS96.pdf",
    "SBS_96_plots_S1.pdf",
    "SBS_96_plots_SBS96_De-Novo.pdf",
    "Samples.txt:md5,6aa8d6af573c4e66a295f92c52e34b9d",
    "Samples.txt:md5,c3dee44caf26c936de8fef73ab167b39",
    "Samples.txt:md5,d0eaf04de705941e78102ef0b52546f5",
    "Seeds.txt:md5,b294ba598152bc73d9e8e54db56459f0",
    "Seeds.txt:md5,b294ba598152bc73d9e8e54db56459f0",
    "Seeds.txt:md5,b294ba598152bc73d9e8e54db56459f0",
    "X_seqinfo.txt",
    "X_seqinfo.txt",
    "X_seqinfo.txt",
    "X_seqinfo.txt",
    "Y_seqinfo.txt",
    "Y_seqinfo.txt",
    "Y_seqinfo.txt",
    "Y_seqinfo.txt",
    "input_data.txt:md5,f8b97a9d1b25dbfb380211b452a506b5",
    "input_data.txt:md5,f8b97a9d1b25dbfb380211b452a506b5",
    "test.DBS1248.all:md5,38184798b4af6ac52c8776ec0716772f",
    "test.DBS150.all:md5,4b29cc951d00fb2180664f7d51ec302e",
    "test.DBS186.all:md5,c4e1c5223084da1f2210034cb09da0db",
    "test.DBS2400.all:md5,4fc604f283263e9a97a1b3c6013e21ce",
    "test.DBS2976.all:md5,275285fb0ee750f284d446b43d520ed2",
    "test.DBS78.all:md5,d0eaf04de705941e78102ef0b52546f5",
    "test.ID28.all:md5,597ee1dfdc82778c3ceee0b1f6afae3f",
    "test.ID332.all:md5,4730dda0c4c7f32afc6d7b65b1685823",
    "test.ID415.all:md5,f83d142841518c8b4051986a6f6eb8b9",
    "test.ID83.all:md5,6aa8d6af573c4e66a295f92c52e34b9d",
    "test.ID8628.all:md5,57ea0750faac3e828bea0c68d0c4ee2c",
    "test.ID96.all:md5,877e611537075e90f05716c6eec7d929",
    "test.SBS1536.all:md5,8be962560a5900145ac808d0b2fcf7e3",
    "test.SBS18.all:md5,a71d51d117195b78f0b947b19e2ceb70",
    "test.SBS24.all:md5,acbdd71c53cde171a39b038424e026b7",
    "test.SBS288.all:md5,80755093742f9e3531804bdf40f75e14",
    "test.SBS384.all:md5,8c79472deb4e78ac1ed78b37c9658fdd",
    "test.SBS4608.all:md5,19ecb85e364aa17917c8ae568c51a88c",
    "test.SBS6.all:md5,7ee16357da652f5e4e1c46c43fbea998",
    "test.SBS6144.all:md5,0e301200b308c463e1dfd0e055d64760",
    "test.SBS96.all:md5,c3dee44caf26c936de8fef73ab167b39",
)


def results_entries(tarball: str) -> list[str]:
    """``name:md5,<hex>`` (or the bare name, for unstable files) of every file under the tarball's results/."""
    unstable = (".pdf", "JOB_METADATA.txt", "seqinfo.txt")
    entries = []
    with tarfile.open(tarball) as tar:
        for member in tar.getmembers():
            name = Path(member.name).name
            if not member.isfile() or name.endswith((".err", ".out")):
                continue
            if name.endswith(unstable):
                entries.append(name)
                continue
            fh = tar.extractfile(member)
            assert fh is not None
            entries.append(f"{name}:md5,{hashlib.md5(fh.read()).hexdigest()}")
    return sorted(entries)


@env.task
async def test_sigprofiler_single_sample() -> None:
    # upstream case: "test_sigprofiler_single_sample"
    seeds_path = Path(tempfile.mkdtemp()) / "seed.txt"
    seeds_path.write_text(SEEDS)
    seeds = await File.from_local(str(seeds_path))
    volume = await sigprofiler_install_genome("GRCh38_Legacy", "GRCh38")
    result = await sigprofiler(
        [await fixture(TABLE)], "GRCh38", genome_installed_path=volume, seeds=seeds, prefix="test", args=ARGS
    )
    entries = await asyncio.to_thread(results_entries, str(await result.results_tar.download()))
    missing = Counter(EXPECTED_RESULTS) - Counter(entries)
    extra = Counter(entries) - Counter(EXPECTED_RESULTS)
    assert not missing and not extra, f"results/ differs: missing {sorted(missing)}, unexpected {sorted(extra)}"
    await assert_md5(result.versions, VERSIONS_MD5, label="versions.yml")


async def expect_failure(args: str, what: str) -> None:
    try:
        await sigprofiler([await fixture(TABLE)], "GRCh38", prefix="test", args=args)
    except Exception as e:  # the task fails: the template rejects its arguments
        print(f"sigprofiler failed as expected ({what}): {e}")
    else:
        raise AssertionError(f"sigprofiler accepted {what}")


@env.task
async def test_sigprofiler_invalid_boolean_argument() -> None:
    # upstream case: "test_sigprofiler - invalid boolean argument"
    await expect_failure(BAD_ARGS.format(download="download_genome_sigprofiler:True"), "make_decomposition_plots:maybe")


@env.task
async def test_sigprofiler_no_download_without_genome_path() -> None:
    # upstream case: "test_sigprofiler - download False but no genome_installed_path"
    await expect_failure(
        BAD_ARGS.format(download="download_genome_sigprofiler:False"), "no download and no genome_installed_path"
    )


tests = [
    test_sigprofiler_single_sample,
    test_sigprofiler_invalid_boolean_argument,
    test_sigprofiler_no_download_without_genome_path,
]
