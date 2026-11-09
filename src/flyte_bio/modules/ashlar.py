"""ashlar — stitch and register multiplexed tissue images into one pyramidal OME-TIFF.

Exposes :func:`ashlar`, which aligns one or more cycles' tiled images (in
cycle order) into ``<prefix>.ome.tif``, optionally correcting illumination
with dark-field (``dfp``) and flat-field (``ffp``) profiles, e.g. from
:func:`flyte_bio.modules.basicpy.basicpy`. Give zero, one, or one profile per
image of each kind.

As upstream, each image is staged as ``image<n>/<name>`` (and the profiles as
``dfp<n>/`` and ``ffp<n>/``), and the random ``UUID`` attributes in the
output's OME-XML are blanked to spaces of the same length, so the output's
bytes are reproducible.
"""

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

ASHLAR_IMAGE = "community.wave.seqera.io/library/ashlar_tifffile:6466db6f351397ba"
DEFAULT_RESOURCES = flyte.Resources(cpu=1, memory="6Gi")

# A list[File] is staged under its basenames and globbed alphabetically, so the
# wrapper passes each list's order as "/"-separated names. Each file is linked
# into its own numbered directory, as upstream stages them.
ashlar_cmd = shell.create(
    name="ashlar",
    image=ASHLAR_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "images": list[File],
        "image_names": str,
        "dfp": list[File],
        "dfp_names": str,
        "ffp": list[File],
        "ffp_names": str,
        "prefix": str,
        "args": str,
    },
    defaults={"dfp": [], "dfp_names": "", "ffp": [], "ffp_names": "", "args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        # The image's entrypoint activates its environment; sourced in a function so that
        # the hook's closing `exec "$@"` sees no arguments.
        activate_environment() {{ . /shell-hook.sh; }}
        activate_environment
        shopt -s nullglob
        IMAGES=({{inputs.images}}); DFP=({{inputs.dfp}}); FFP=({{inputs.ffp}})
        shopt -u nullglob
        IMAGE_NAMES={{inputs.image_names}}
        DFP_NAMES={{inputs.dfp_names}}
        FFP_NAMES={{inputs.ffp_names}}
        P={{inputs.prefix}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"
        export JAVA_TOOL_OPTIONS='-XX:+PerfDisableSharedMem'
        export HOME="$W" XDG_CACHE_HOME="$W/.cache" MPLCONFIGDIR="$W/.matplotlib"

        # stage NAME_LIST STAGE_PREFIX FIRST_FILE: link the files in order, filling STAGED
        stage() {{
            STAGED=()
            local names="$1" dir_prefix="$2" source_dir n=1
            [ -n "$names" ] || return 0
            source_dir=$(dirname "$3")
            IFS=/ read -ra ORDERED <<< "$names"
            for NAME in "${{ORDERED[@]}}"; do
                mkdir "$dir_prefix$n"
                ln -s "$source_dir/$NAME" "$dir_prefix$n/$NAME"
                STAGED+=("$dir_prefix$n/$NAME")
                n=$(( n + 1 ))
            done
        }}
        stage "$IMAGE_NAMES" image "${{IMAGES[0]}}"; INPUTS=("${{STAGED[@]}}")
        PROFILES=()
        if [ ${{#DFP[@]}} -gt 0 ]; then stage "$DFP_NAMES" dfp "${{DFP[0]}}"; PROFILES+=(--dfp "${{STAGED[@]}}"); fi
        if [ ${{#FFP[@]}} -gt 0 ]; then stage "$FFP_NAMES" ffp "${{FFP[0]}}"; PROFILES+=(--ffp "${{STAGED[@]}}"); fi

        ashlar -o "$P.ome.tif" "${{INPUTS[@]}}" $ARGS "${{PROFILES[@]}}" || exit 1
        UUID='UUID="urn:uuid:[[:xdigit:]]{{8}}-[[:xdigit:]]{{4}}-[[:xdigit:]]{{4}}-[[:xdigit:]]{{4}}-[[:xdigit:]]{{12}}"'
        sed -i -E "s/$UUID/                                                    /g" "$P.ome.tif" || exit 1
        mv "$P.ome.tif" "$OUT"/
    """,
)

env = flyte.TaskEnvironment.from_task("ashlar", ashlar_cmd.as_task())


def ordered_names(files: list[File], label: str) -> str:
    names = [f.name for f in files]
    if len(set(names)) != len(names):
        raise ValueError(f"ashlar: the {label} must have distinct file names, got {names}")
    if any("/" in n for n in names):
        raise ValueError(f"ashlar: {label} names can't contain '/': {names}")
    return "/".join(names)


async def ashlar(
    images: list[File],
    dfp: list[File] | None = None,
    ffp: list[File] | None = None,
    prefix: str = "sample",
    args: str = "",
) -> File:
    """Stitch and register ``images`` (one per cycle, in order) into ``<prefix>.ome.tif``.

    ``dfp`` and ``ffp`` each hold zero, one, or ``len(images)`` illumination profiles.
    """
    dfp, ffp = dfp or [], ffp or []
    if not images:
        raise ValueError("ashlar: give at least one image")
    for profiles, label in ((dfp, "dfp"), (ffp, "ffp")):
        if len(profiles) not in (0, 1, len(images)):
            raise ValueError(f"ashlar: give zero, one, or {len(images)} {label} files, got {len(profiles)}")
    results = await ashlar_cmd(
        images=images,
        image_names=ordered_names(images, "images"),
        dfp=dfp,
        dfp_names=ordered_names(dfp, "dfp profiles"),
        ffp=ffp,
        ffp_names=ordered_names(ffp, "ffp profiles"),
        prefix=prefix,
        args=args,
    )
    found = await results.get_file(f"{prefix}.ome.tif")
    if found is None:
        raise FileNotFoundError(f"ashlar wrote no {prefix}.ome.tif")
    return found
