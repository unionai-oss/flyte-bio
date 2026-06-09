"""Vendored helper scripts that tool biocontainers don't ship themselves.

Some upstream tool wrappers aren't a single CLI invocation — they drive a
small helper (Python/Perl/R) that lives alongside the tool definition, not
inside the tool's container image. We vendor those helpers here as **real
script files** (core logic preserved; the original Nextflow ``template``
entry point swapped for a plain argparse CLI).

Every one runs the same way: a **shell task** on a stock interpreter
biocontainer, with the script handed in as a normal **File input** (staged
into the container like any input). That keeps tool pods self-contained —
no custom image build, and crucially no pod ever needs ``flyte_bio``
installed (a native Flyte task would, because the worker imports its
defining module; a shell task is just a serialized command). The wrapper in
each module materializes the script File from :func:`path`.

The catch: a script file used only as data isn't an imported module, so
Flyte's default ``loaded_modules`` code bundle would drop it.
:func:`ship_in_bundle` fixes that without forcing ``--copy-style all``:
Flyte bundles the files backing ``sys.modules`` entries (and applies no
``.py`` filter — see
``flyte._code_bundle._utils.list_imported_modules_as_files``), so we
register each script as a stub module whose ``__file__`` points at it. The
bundler then copies the file at its path relative to the project root, so
:func:`path` resolves in our dev runs. (For an external consumer the script
ships inside the installed wheel, so :func:`path` resolves in site-packages
and no bundling is needed.)
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


# Vendored scripts that ride into a biocontainer as File inputs. Listed here
# (not imported anywhere) so the default code bundle still ships them.
ship_in_bundle("gtffilter.py", "catadditionalfasta.py", "gtf2bed.pl")


__all__ = ["path", "ship_in_bundle"]
