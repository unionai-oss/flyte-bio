"""Vendored helper scripts that tool biocontainers don't ship themselves.

Some upstream tool wrappers aren't a single CLI invocation — they drive a
small helper that lives alongside the tool definition, not inside the
tool's container image. We vendor those helpers here as **real script
files** (core logic preserved; the original templated entry point swapped
for a plain argv/argparse interface), in two flavours:

- **Pure-Python, stdlib-only helpers** (e.g. :mod:`gtffilter`,
  :mod:`catadditionalfasta`) are imported and called directly from a native
  Flyte ``@env.task`` — no biocontainer involved.

- **Perl / R helpers** (e.g. ``gtf2bed.pl``) need their interpreter from the
  tool's biocontainer. A biocontainer never receives Flyte's code bundle
  (only native-Python task pods do), so we hand the script to the shell task
  as a normal **File input**, staged into the container like any input — no
  custom image build. The caller materializes that File from :func:`path`.

The catch: a loose ``.pl``/``.r`` isn't an imported module, so Flyte's
default ``loaded_modules`` code bundle would drop it. :func:`ship_in_bundle`
fixes that without forcing ``--copy-style all`` on every run: Flyte bundles
the files backing ``sys.modules`` entries (and applies no ``.py`` filter —
see ``flyte._code_bundle._utils.list_imported_modules_as_files``), so we
register each script as a stub module whose ``__file__`` points at it. The
bundler then copies the file at its path relative to the project root, so
:func:`path` resolves both locally and in-cluster.
"""



import sys
import types
from pathlib import Path

SCRIPTS_DIR = Path(__file__).resolve().parent


def path(name: str) -> Path:
    """Return the on-disk path of vendored script ``name``."""
    return SCRIPTS_DIR / name


def ship_in_bundle(*names: str) -> None:
    """Make Flyte's default code bundle ship non-Python vendored scripts.

    Registers each script as a stub ``sys.modules`` entry whose ``__file__``
    points at it, so the ``loaded_modules`` bundler copies it like any
    imported module — no ``--copy-style all`` needed. Idempotent.
    """
    for name in names:
        key = f"{__name__}._vendored_{name.replace('.', '_')}"
        if key in sys.modules:
            continue
        stub = types.ModuleType(key)
        stub.__file__ = str(SCRIPTS_DIR / name)
        sys.modules[key] = stub


# Vendored non-Python scripts that ride into a biocontainer as File inputs.
ship_in_bundle("gtf2bed.pl")


__all__ = ["path", "ship_in_bundle"]
