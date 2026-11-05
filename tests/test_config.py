"""Tests for flyte_bio.config and the pipelines' Hydra configs in flyte_bio/conf."""

from flyte.io import Dir, File
from hydra import compose, initialize_config_dir
from omegaconf import OmegaConf

from flyte_bio.config import CONFIG_DIR, bind_config, pipeline_defaults, run_pipeline
from flyte_bio.pipelines.rnaseq import rnaseq
from flyte_bio.pipelines.scrnaseq import scrnaseq
from flyte_bio.pipelines.variant_calling import variant_calling
from tests.framework import env

PIPELINES = (rnaseq, scrnaseq, variant_calling)


def compose_config(name: str, overrides: list[str]):
    with initialize_config_dir(config_dir=str(CONFIG_DIR / name), version_base=None):
        return compose(config_name=name, overrides=overrides)


def expect_error(cfg, text: str) -> None:
    """bind_config(scrnaseq, cfg) raises a ValueError that mentions ``text``."""
    try:
        bind_config(scrnaseq, cfg)
    except ValueError as e:
        assert text in str(e), f"expected {text!r} in: {e}"
    else:
        raise AssertionError(f"expected ValueError mentioning {text!r}")


@env.task
async def test_config_files() -> None:
    # Each pipeline's config holds exactly its options and defaults, and its
    # test profile only sets options that exist, with the right types.
    for pipeline in PIPELINES:
        name = pipeline.__name__
        options = OmegaConf.load(CONFIG_DIR / name / f"{name}.yaml")
        options.pop("defaults")
        assert options == pipeline_defaults(pipeline), f"{name}.yaml is out of date with {name}()"

        cfg = compose_config(name, ["profile=test"])
        kwargs = bind_config(pipeline, cfg, samplesheet=File.from_existing_remote("s3://bucket/samplesheet.csv"))
        assert isinstance(kwargs["fasta"], File), kwargs["fasta"]


@env.task
async def test_bind_config() -> None:
    required = {"samplesheet": "s3://b/s.csv", "fasta": "s3://b/g.fa", "gtf": "s3://b/g.gtf", "protocol": "10XV2"}
    kwargs = bind_config(
        scrnaseq,
        OmegaConf.create({**required, "simpleaf_index": "s3://b/index", "qcatch_n_partitions": 4, "hydra": {}}),
    )
    assert isinstance(kwargs["fasta"], File) and kwargs["fasta"].path == "s3://b/g.fa"
    assert isinstance(kwargs["simpleaf_index"], Dir)
    assert kwargs["qcatch_n_partitions"] == 4
    assert "skip_cellbender" not in kwargs  # left to the function's default

    kwargs = bind_config(
        rnaseq,
        {
            "samplesheet": "s3://b/s.csv",
            "fasta": "s3://b/g.fa",
            "rseqc_modules": ["bam_stat"],
            "min_mapped_reads": 5,
            "bbsplit_fasta_list": {"human": "s3://b/human.fa"},
        },
    )
    assert kwargs["rseqc_modules"] == ("bam_stat",)
    assert kwargs["min_mapped_reads"] == 5.0
    assert kwargs["bbsplit_fasta_list"]["human"].path == "s3://b/human.fa"

    kwargs = bind_config(
        variant_calling, {"samplesheet": "s3://b/s.csv", "fasta": "s3://b/g.fa", "known_indels": ["s3://b/k.vcf.gz"]}
    )
    assert [f.path for f in kwargs["known_indels"]] == ["s3://b/k.vcf.gz"]

    # Given arguments win over the config, and None ones are ignored.
    sheet = File.from_existing_remote("s3://b/given.csv")
    assert bind_config(scrnaseq, required, samplesheet=sheet)["samplesheet"] is sheet
    assert bind_config(scrnaseq, required, samplesheet=None)["samplesheet"].path == "s3://b/s.csv"

    expect_error({**required, "skip_cellbnder": True}, "'skip_cellbnder' (did you mean 'skip_cellbender'?)")
    expect_error({**required, "skip_cellbender": "yes"}, "skip_cellbender: expected bool, got 'yes'")
    expect_error({**required, "qcatch_n_partitions": True}, "qcatch_n_partitions: expected int, got True")
    expect_error({**required, "fasta": None}, "fasta: can't be null")
    expect_error({k: v for k, v in required.items() if k != "gtf"}, "gtf: required")
    expect_error(compose_config("scrnaseq", []), "protocol: required")
    # Tests run in the cluster, which can't read a local path.
    expect_error({**required, "fasta": "genome.fa"}, "is a local path")


@env.task
async def test_run_pipeline() -> None:
    async def pipeline(samplesheet: File, n: int = 1, flag: bool = False) -> tuple[str, int, bool]:
        return samplesheet.path, n, flag

    assert await run_pipeline(pipeline, {"samplesheet": "s3://b/s.csv", "n": 3}, flag=True) == ("s3://b/s.csv", 3, True)


tests = [test_config_files, test_bind_config, test_run_pipeline]
