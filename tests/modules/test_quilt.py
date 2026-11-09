"""Tests for flyte_bio.modules.quilt.

The cases and expected values are upstream's quilt/quilt tests (modules commit
9a1cc9f, QUILT 1.0.5), run on 1 CPU as upstream's test config pins. Upstream
snapshots the imputed VCF's ``variantsMD5`` (the md5 of its record lines, as
nft-vcf computes it), the ``.tbi``'s md5, and the names of the RData and plot
files. As upstream, the phasefile and genfile cases first make NA12878's truth
file from its genotypes: ``bcftools query -f '[%GT]\\n' -r chr22``, then an awk
program that adds the sample name (and, for the genfile, codes genotypes as
0, 1, 2 or NA). The genfile case also checks that QUILT reported the
imputation accuracy (upstream reads it from ``.command.err``; here it's
QUILT's log).
"""

import asyncio
import gzip
import hashlib
import tempfile
from pathlib import Path

from flyte.io import File

from flyte_bio.modules.bcftools import bcftools_query
from flyte_bio.modules.gawk import gawk
from flyte_bio.modules.quilt import QuiltResult, quilt_quilt
from tests.framework import assert_md5, env, fixture

HS = "genomics/homo_sapiens/"
BAM = HS + "illumina/bam/NA12878.chr21_22.1X.bam"
SAMPLE_VCF = HS + "illumina/vcf/NA12878.chr21_22.1X.glimpse2.vcf.gz"
PANEL = HS + "popgen/1000GP.chr22"
REGION = {"chromosome": "chr22", "region_start": 16570000, "region_end": 16610000, "ngen": 100, "buffer": 10000}
PREPARED_REFERENCE = "QUILT_prepared_reference.chr22.16570000.16610000.RData"
# Upstream's genfile program: the sample name, then each genotype coded 0, 1, 2 or NA.
GENFILE_PROGRAM = r"""BEGIN { OFS="\t"; print "NA12878" }
{
    for (i=1; i<=NF; i++) {
        if ($i ~ /\./) {
            $i = "NA"
        } else {
            n = split($i, a, /[\/|]/)
            if (n != 2)
                $i = "NA"
            else if (a[1] == 0 && a[2] == 0)
                $i = 0
            else if (a[1] == a[2])
                $i = 2
            else
                $i = 1
        }
    }
    print
}"""


async def quilt_write(name: str, text: str) -> File:
    path = Path(tempfile.mkdtemp()) / name
    path.write_text(text)
    return await File.from_local(str(path))


async def quilt_vcf_lines(file: File) -> list[str]:
    async with file.open("rb") as fh:
        return gzip.decompress(bytes(await fh.read())).decode().splitlines()


async def quilt_check(result: QuiltResult, prefix: str, variants: str, tbi: str, rdata: list[str]) -> None:
    """Check the VCF's name and variantsMD5, the .tbi's md5 and the RData files' names."""
    assert result.vcf.name == f"{prefix}.vcf.gz", result.vcf.name
    records = [line for line in await quilt_vcf_lines(result.vcf) if not line.startswith("#")]
    actual = format(int(hashlib.md5("".join(records).encode()).hexdigest(), 16), "x")
    assert actual == variants, f"variantsMD5 {actual}, expected {variants}"
    assert result.tbi is not None and result.tbi.name == f"{prefix}.vcf.gz.tbi", "QUILT wrote no .tbi"
    await assert_md5(result.tbi, tbi, label="tbi")
    names = [f.name for f in result.rdata]
    assert names == rdata, f"RData: {names}"


async def quilt_inputs() -> dict:
    bam, bai, hap, legend = await asyncio.gather(
        fixture(BAM), fixture(BAM + ".bai"), fixture(PANEL + ".hap.gz"), fixture(PANEL + ".legend.gz")
    )
    return {"reads": [bam], "indexes": [bai], "haplotypes": hap, "legend": legend, **REGION}


async def quilt_truth_inputs(program: str, suffix: str) -> dict:
    """The inputs of upstream's truth cases: posfile, genetic map, reference, and the truth file from ``program``."""
    vcf, csi, posfile, genetic_map, fasta, fai = await asyncio.gather(
        fixture(SAMPLE_VCF),
        fixture(SAMPLE_VCF + ".csi"),
        fixture(PANEL + ".posfile"),
        fixture(HS + "genome/genetic_map/genome.GRCh38.chr22.stitch.map"),
        fixture(HS + "genome/genome.fasta"),
        fixture(HS + "genome/genome.fasta.fai"),
    )
    genotypes = await bcftools_query(vcf, "NA12878", index=csi, format=r"[%GT]\n", args="-r chr22")
    truth = await gawk([genotypes], program=program, prefix="NA12878", suffix=suffix)
    assert truth.output is not None, "gawk wrote no truth file"
    return {
        "posfile": posfile,
        "genetic_map": genetic_map,
        "fasta": fasta,
        "fasta_index": fai,
        "truth": truth.output,
        **await quilt_inputs(),
    }


@env.task
async def test_quilt_phasefile() -> None:
    # upstream case: quilt/quilt "homo_sapiens - bam, hap, legend, posfile, phasefile, map - fasta"
    inputs = await quilt_truth_inputs('BEGIN{print "NA12878"} {print}', "phasefile.name.txt")
    phasefile = inputs.pop("truth")
    result = await quilt_quilt(
        prefix="NA12878",
        phasefile=phasefile,
        args="--save_prepared_reference=TRUE --make_plots=TRUE --seed=1",
        **inputs,
    )
    await quilt_check(
        result, "NA12878", "df1391e233c9a8b237e0f77e188ee2b7", "f4e800b09116b8ccac8e62c0f835eeed", [PREPARED_REFERENCE]
    )
    expected = sorted(
        f"haps.NA12878.chr22.16570000.16610000_igs.{i}.{stage}.png"
        for i in range(1, 9)
        for stage in ("0.truth", "it1.gibbs", "it2.gibbs", "it3.gibbs")
    )
    names = [f.name for f in result.plots]
    assert names == expected, f"plots: {names}"


@env.task
async def test_quilt_genfile() -> None:
    # upstream case: quilt/quilt "homo_sapiens - bam, hap, legend, posfile, genfile, map - fasta"
    inputs = await quilt_truth_inputs(GENFILE_PROGRAM, "genfile.name.txt")
    genfile = inputs.pop("truth")
    result = await quilt_quilt(
        prefix="NA12878", genfile=genfile, args="--save_prepared_reference=TRUE --seed=1", **inputs
    )
    await quilt_check(
        result, "NA12878", "288ad019a68e19a83e797ff2a48fcece", "4e2bcdc3d5d3b113f66f3a1b87e8f0b9", [PREPARED_REFERENCE]
    )
    assert result.plots == [], "QUILT made plots without --make_plots"
    async with result.log.open("rb") as fh:
        log = bytes(await fh.read()).decode()
    assert "Final imputation dosage accuracy for sample NA12878" in log, "QUILT reported no accuracy"


@env.task
async def test_quilt_no_optional_files() -> None:
    # upstream case: quilt/quilt "QUILT no optional files"
    result = await quilt_quilt(prefix="test", args="--seed=1", **await quilt_inputs())
    await quilt_check(result, "test", "bb640c38eea9b8df18d1debad11a21e0", "4fed28161b3da64b2f8302e3cbb4c0d9", [])
    assert result.plots == [], "QUILT made plots without --make_plots"


@env.task
async def test_quilt_bamlist_renaming() -> None:
    # upstream case: quilt/quilt "QUILT with bamlist and renaming"
    bamlist, names = await asyncio.gather(
        quilt_write("bamlist.txt", "NA12878.chr21_22.1X.bam\n"), quilt_write("bamnames.txt", "Mysample1\n")
    )
    result = await quilt_quilt(
        prefix="test", bamlist=bamlist, sample_names=names, args="--seed=1", **await quilt_inputs()
    )
    await quilt_check(result, "test", "bb640c38eea9b8df18d1debad11a21e0", "487166bf9532223bbde1156b80e5cfa3", [])
    header = next(line for line in await quilt_vcf_lines(result.vcf) if line.startswith("#CHROM"))
    assert header.split("\t")[9:] == ["Mysample1"], f"samples: {header.split(chr(9))[9:]}"


tests = [test_quilt_phasefile, test_quilt_genfile, test_quilt_no_optional_files, test_quilt_bamlist_renaming]
