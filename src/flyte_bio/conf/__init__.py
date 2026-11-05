"""Hydra configs for each pipeline.

``<pipeline>/<pipeline>.yaml`` holds every option with its default, and
``<pipeline>/profile/test.yaml`` upstream's test data.

Importing this package registers each file as a stub module, so Flyte's
default code bundle ships them in dev runs, as ``flyte_bio.scripts`` does
for its scripts. An installed flyte-bio has them in its wheel.
"""

import sys
import types
from pathlib import Path

CONFIG_DIR = Path(__file__).resolve().parent

for config_file in sorted(CONFIG_DIR.glob("*/**/*.yaml")):
    key = f"{__name__}.bundled_" + "_".join(config_file.relative_to(CONFIG_DIR).with_suffix("").parts)
    if key not in sys.modules:
        stub = types.ModuleType(key)
        stub.__file__ = str(config_file)
        sys.modules[key] = stub
