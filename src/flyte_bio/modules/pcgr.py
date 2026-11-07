"""pcgr — download the Personal Cancer Genome Reporter's reference data bundle.

Exposes :func:`pcgr_getref`, which downloads PCGR's reference bundle for a
bundle version (such as ``20250314``) and genome (``grch37`` or ``grch38``)
from PCGR's server, ``pcgr_ref_data.<version>.<genome>.tgz``, and unpacks it
into ``<version>/data``, as upstream. That folder is what PCGR's
``--refdata_dir`` takes.

It's returned as ``<version>.tar.gz``, an archive of the ``<version>``
folder, rather than a Dir: the bundle holds hundreds of files, too many for a
Dir input to stage reliably. The archive starts with the bundle's version
file, ``<version>/data/<genome>/.PCGR_BUNDLE_VERSION``. The download is a
task of its own, cached on the version and genome, so each bundle (about
5 GB packed) downloads once.
"""

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

PCGR_GETREF_IMAGE = "community.wave.seqera.io/library/coreutils_curl_gzip_tar:17a6ea9a6766c02a"
PCGR_BUNDLE_URL = "https://insilico.hpc.uio.no/pcgr"
# Upstream's process_single. The task holds the bundle twice, unpacked and
# repacked (about 5 GB packed), so it needs the disk.
DOWNLOAD_RESOURCES = flyte.Resources(cpu=1, memory="2Gi", disk="40Gi")

pcgr_getref_cmd = shell.create(
    name="pcgr_getref",
    image=PCGR_GETREF_IMAGE,
    resources=DOWNLOAD_RESOURCES,
    inputs={"bundle_version": str, "genome": str},
    outputs={"results": Dir},
    script=rf"""
        set -o pipefail
        V={{inputs.bundle_version}}
        G={{inputs.genome}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"

        curl -fsSL "{PCGR_BUNDLE_URL}/pcgr_ref_data.$V.$G.tgz" | gzip -dc | tar xf - || exit 1
        if [ ! -d data ]; then echo "the PCGR bundle has no data/ folder" >&2; exit 1; fi
        mkdir "$V"
        mv data/ "$V"

        FIRST="$V/data/$G/.PCGR_BUNDLE_VERSION"
        if [ ! -f "$FIRST" ]; then echo "the PCGR bundle has no $FIRST" >&2; exit 1; fi
        shopt -s globstar dotglob
        {{
            echo "$FIRST"
            echo "$V"
            for F in "$V"/**; do
                if [ "$F" != "$FIRST" ] && [ "$F" != "$V/" ]; then echo "$F"; fi
            done
        }} > members.txt
        tar --no-recursion -czf "$OUT/$V.tar.gz" -T members.txt
    """,
)

env = flyte.TaskEnvironment.from_task("pcgr", pcgr_getref_cmd.as_task())


async def pcgr_getref(bundle_version: str, genome: str) -> File:
    """PCGR's reference bundle ``bundle_version`` for ``genome``: ``<bundle_version>.tar.gz`` of ``<version>/data``."""
    results = await pcgr_getref_cmd(bundle_version=bundle_version, genome=genome)
    archive = await results.get_file(f"{bundle_version}.tar.gz")
    if archive is None:
        raise FileNotFoundError(f"pcgr_getref wrote no {bundle_version}.tar.gz")
    return archive
