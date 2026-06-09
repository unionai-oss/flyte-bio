"""flyte_bio.modules — typed shell-task wrappers around bio CLI tools.

Each submodule wraps one tool family, sharing a single biocontainer image
across its tasks:

- :mod:`flyte_bio.modules.bedtools` — genome arithmetic (intersect, sort,
  merge).
- :mod:`flyte_bio.modules.cat` — file concatenation (``cat_fastq``).
- :mod:`flyte_bio.modules.catadditionalfasta` — append an extra FASTA + GTF.
- :mod:`flyte_bio.modules.gffread` — GFF/GTF conversion + transcript FASTA.
- :mod:`flyte_bio.modules.gtf2bed` — derive a BED12 gene model from a GTF.
- :mod:`flyte_bio.modules.gtffilter` — restrict a GTF to a genome's sequences.
- :mod:`flyte_bio.modules.gunzip` — single-file gzip decompression.
- :mod:`flyte_bio.modules.samtools` — alignment sort/index/stats/faidx.
- :mod:`flyte_bio.modules.untar` — tar archive extraction.

The module-level :data:`env` here is an aggregate
:class:`flyte.TaskEnvironment` depending on every submodule's env.
Pipelines depend on it once to gain access to every wrapped tool::

    from flyte_bio.modules import env as modules_env
"""



import flyte

from .bedtools import env as bedtools_env
from .cat import env as cat_env
from .catadditionalfasta import env as catadditionalfasta_env
from .gffread import env as gffread_env
from .gtf2bed import env as gtf2bed_env
from .gtffilter import env as gtffilter_env
from .gunzip import env as gunzip_env
from .samtools import env as samtools_env
from .untar import env as untar_env

env = flyte.TaskEnvironment(
    name="flyte_bio_modules",
    depends_on=[
        bedtools_env,
        cat_env,
        catadditionalfasta_env,
        gffread_env,
        gtf2bed_env,
        gtffilter_env,
        gunzip_env,
        samtools_env,
        untar_env,
    ],
)

__all__ = [
    "bedtools_env",
    "cat_env",
    "catadditionalfasta_env",
    "env",
    "gffread_env",
    "gtf2bed_env",
    "gtffilter_env",
    "gunzip_env",
    "samtools_env",
    "untar_env",
]
