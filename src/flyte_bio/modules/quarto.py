"""quarto — render a parameterised notebook to HTML with Quarto.

Exposes :func:`quarto_notebook`, which executes and renders one notebook
(Quarto ``.qmd``, R Markdown ``.Rmd`` or Jupyter ``.ipynb``) to
``<prefix>.html``, as upstream's quarto/notebook module and its predecessor
quartonotebook do. The two differ only in the notebook's versions:
quarto/notebook (``notebook_versions=True``, the default) requires the
notebook to write a ``versions.csv`` of ``tool,version`` rows and turns it into
``versions.yml``; quartonotebook (``notebook_versions=False``) doesn't.

The notebook gets its parameters from ``params.yml`` (``quarto render
--execute-params``), as upstream: ``meta`` (default ``{"id": prefix}``),
``cpus`` (the task's CPUs, read when it runs), ``artifact_dir``
(``"artifacts"``, created empty for the notebook to write into), then
``parameters``, which may override any of them. The wrapper writes the YAML,
so parameter values never pass through the shell. ``input_files`` are linked
into the working directory, where the notebook runs, and ``extensions`` (a
Quarto ``_extensions`` directory) is copied there.

The results Dir holds, under upstream's names: ``<prefix>.html``, the
notebook, ``params.yml``, ``versions.yml`` (with ``notebook_versions``), the
artifact directory with whatever the notebook wrote into it, and
``_extensions`` when present. A notebook that writes hundreds of artifacts
should pack them into an archive: very large output directories fail to
upload.

The image is upstream's (Quarto, Papermill, Jupyter with matplotlib, R with
rmarkdown); a notebook that needs other packages needs an image with them
(``quarto_notebook_cmd`` with another ``image``). BLAS, OpenMP and numba
threads are the task's CPUs.
"""

import json
import tempfile
from dataclasses import dataclass
from pathlib import Path

import flyte
from flyte.extras import shell
from flyte.io import Dir, File

from flyte_bio.limits import TASK_LIMITS

QUARTO_IMAGE = "community.wave.seqera.io/library/jupyter_matplotlib_papermill_quarto_r-rmarkdown:6d15193ce3dfc665"
# Upstream's process_low.
DEFAULT_RESOURCES = flyte.Resources(cpu=2, memory="8Gi")

# Written as the value of the implicit ``cpus`` parameter; the script replaces it with the task's CPUs.
RUNTIME_CPUS = "RUNTIME_CPUS"

quarto_notebook_cmd = shell.create(
    name="quarto_notebook",
    image=QUARTO_IMAGE,
    resources=DEFAULT_RESOURCES,
    inputs={
        "notebook": File,
        "params_yaml": File,
        "input_files": list[File],
        "extensions": Dir | None,
        "artifact_dir": str,
        "prefix": str,
        "notebook_versions": bool,
        "versions_key": str,
        "args": str,
    },
    defaults={"notebook_versions": True, "versions_key": "quarto_notebook", "args": ""},
    outputs={"results": Dir},
    script=rf"""
        {TASK_LIMITS}
        NB=({{inputs.notebook}})
        PARAMS=({{inputs.params_yaml}})
        shopt -s nullglob; INPUTS=({{inputs.input_files}}); shopt -u nullglob
        EXTENSIONS={{inputs.extensions}}
        ARTIFACTS={{inputs.artifact_dir}}
        P={{inputs.prefix}}
        KEY={{inputs.versions_key}}
        ARGS={{inputs.args}}
        OUT={{outputs.results}}
        W=$(mktemp -d); cd "$W"

        for F in "${{INPUTS[@]}}"; do ln -s "$F" "$(basename "$F")"; done
        if [ -d "$EXTENSIONS" ] && [ -n "$(ls -A "$EXTENSIONS")" ]; then cp -rL "$EXTENSIONS" _extensions; fi
        NAME=$(basename "${{NB[0]}}")
        cp "${{NB[0]}}" "$NAME"
        sed "s/^cpus: {RUNTIME_CPUS}$/cpus: $TASK_CPUS/" "${{PARAMS[0]}}" > params.yml
        mkdir -p "$ARTIFACTS"

        export HOME="$W" XDG_CACHE_HOME="$W/.xdg_cache_home" XDG_DATA_HOME="$W/.xdg_data_home"
        export MPLCONFIGDIR="$W/.matplotlib"
        # Without a _quarto.yml here, Quarto searches parent directories for a project.
        [ -e _quarto.yml ] || printf 'project:\n  type: default\n' > _quarto.yml
        ENV_QUARTO=/opt/conda/etc/conda/activate.d/quarto.sh
        set +u
        if [ -z "${{QUARTO_DENO:-}}" ] && [ -f "$ENV_QUARTO" ]; then source "$ENV_QUARTO"; fi
        export MKL_NUM_THREADS=$TASK_CPUS OPENBLAS_NUM_THREADS=$TASK_CPUS
        export OMP_NUM_THREADS=$TASK_CPUS NUMBA_NUM_THREADS=$TASK_CPUS

        quarto render "$NAME" $ARGS --execute-params params.yml --output "$P.html" || exit 1

        if [ {{inputs.notebook_versions}} = true ]; then
            if [ ! -f versions.csv ]; then
                echo "ERROR: versions.csv not found; the notebook must write [tool,version] pairs." >&2
                exit 1
            fi
            printf '"%s":\n%s\n' "$KEY" "$(awk -F',' '{{printf "    %s: %s\n", $1, $2}}' versions.csv)" > versions.yml
        fi

        for F in ./*.html; do [ -L "$F" ] || mv "$F" "$OUT"/; done
        cp "$NAME" params.yml "$OUT"/
        if [ -f versions.yml ]; then mv versions.yml "$OUT"/; fi
        mkdir -p "$OUT/$ARTIFACTS"
        if [ -n "$(ls -A "$ARTIFACTS")" ]; then mv "$ARTIFACTS"/* "$OUT/$ARTIFACTS"/; fi
        if [ -d _extensions ]; then cp -rL _extensions "$OUT"/; fi
        true
    """,
)

env = flyte.TaskEnvironment.from_task("quarto", quarto_notebook_cmd.as_task())


def yaml_scalar(value: object) -> str:
    if value is None:
        return "null"
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, int | float):
        return str(value)
    if isinstance(value, dict):
        return "{}"
    if isinstance(value, list | tuple):
        return "[]"
    return json.dumps(str(value))


def yaml_lines(value: dict | list, indent: int) -> list[str]:
    """Block YAML for ``value``, laid out as upstream's YAML builder writes it (lists not indented)."""
    pad = " " * indent
    lines = []
    items = value.items() if isinstance(value, dict) else ((None, v) for v in value)
    for key, item in items:
        if not (isinstance(item, dict | list | tuple) and len(item) > 0):
            lines.append(f"{pad}{key}: {yaml_scalar(item)}" if key is not None else f"{pad}- {yaml_scalar(item)}")
            continue
        child = list(item) if isinstance(item, tuple) else item
        if key is None:
            inner = yaml_lines(child, indent + 2)
            lines.append(f"{pad}- {inner[0].lstrip()}")
            lines.extend(inner[1:])
        else:
            lines.append(f"{pad}{key}:")
            lines.extend(yaml_lines(child, indent + 2 if isinstance(child, dict) else indent))
    return lines


def notebook_parameters(meta: dict, parameters: dict | None) -> dict:
    """The notebook's parameters: upstream's implicit ones, overridden by ``parameters``."""
    return {"meta": meta, "cpus": RUNTIME_CPUS, "artifact_dir": "artifacts"} | (parameters or {})


def notebook_parameters_yaml(merged: dict) -> str:
    """``params.yml``'s text; an implicit ``cpus`` is filled in when the task runs."""
    lines = ["---"]
    for line in yaml_lines(merged, 0):
        lines.append("cpus: " + RUNTIME_CPUS if line == "cpus: " + json.dumps(RUNTIME_CPUS) else line)
    return "\n".join(lines) + "\n"


async def notebook_parameters_file(merged: dict) -> File:
    path = Path(tempfile.mkdtemp()) / "params.yml"
    path.write_text(notebook_parameters_yaml(merged))
    return await File.from_local(str(path))


@dataclass
class QuartoNotebookResult:
    html: File  # <prefix>.html
    notebook: File  # the notebook as rendered
    params_yaml: File  # params.yml
    versions: File | None  # versions.yml, with notebook_versions
    artifacts: list[File]  # what the notebook wrote into its artifact directory
    extensions: list[File]  # the _extensions directory's files, when there is one
    results: Dir  # all of the above


async def quarto_notebook_result(
    results: Dir, prefix: str, notebook_name: str, artifact_dir: str = "artifacts"
) -> QuartoNotebookResult:
    """Pick :func:`quarto_notebook`'s outputs out of ``quarto_notebook_cmd``'s results."""
    root = results.path.rstrip("/") + "/"
    files = {f.path.removeprefix(root): f async for f in results.walk()}

    def pick(name: str) -> File:
        if name not in files:
            raise FileNotFoundError(f"quarto wrote no {name}")
        return files[name]

    def under(directory: str) -> list[File]:
        return [files[name] for name in sorted(files) if name.startswith(directory.rstrip("/") + "/")]

    return QuartoNotebookResult(
        html=pick(f"{prefix}.html"),
        notebook=pick(notebook_name),
        params_yaml=pick("params.yml"),
        versions=files.get("versions.yml"),
        artifacts=under(artifact_dir),
        extensions=under("_extensions"),
        results=results,
    )


async def quarto_notebook(
    notebook: File,
    meta: dict | None = None,
    parameters: dict | None = None,
    input_files: list[File] | None = None,
    extensions: Dir | None = None,
    prefix: str = "sample",
    notebook_versions: bool = True,
    versions_key: str = "quarto_notebook",
    args: str = "",
) -> QuartoNotebookResult:
    """Render ``notebook`` with Quarto to ``<prefix>.html``; ``versions_key`` heads ``versions.yml``."""
    merged = notebook_parameters(meta if meta is not None else {"id": prefix}, parameters)
    results = await quarto_notebook_cmd(
        notebook=notebook,
        params_yaml=await notebook_parameters_file(merged),
        input_files=input_files or [],
        extensions=extensions,
        artifact_dir=str(merged["artifact_dir"]),
        prefix=prefix,
        notebook_versions=notebook_versions,
        versions_key=versions_key,
        args=args,
    )
    name = notebook.path.rstrip("/").rsplit("/", 1)[-1]
    return await quarto_notebook_result(results, prefix, name, str(merged["artifact_dir"]))
