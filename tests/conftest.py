"""Shared fixtures: load a notebook or a script as a module, and keep tests in a temp folder."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
NOTEBOOKS_DIR = ROOT / "notebooks"


def load_module(path: Path, name: str):
    """Import the Python file at ``path`` under the module name ``name``."""
    if name in sys.modules:
        return sys.modules[name]
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def notebook_pages() -> list[Path]:
    """Return the folder of every notebook that has a ``meta.toml``."""
    return sorted(path.parent for path in NOTEBOOKS_DIR.glob("*/meta.toml"))


@pytest.fixture(autouse=True)
def _work_in_tmp_path(tmp_path, monkeypatch):
    """Run every test inside its own temp folder, so none writes into the checkout."""
    monkeypatch.chdir(tmp_path)


@pytest.fixture(scope="session")
def website_build():
    """The site builder, ``website/build.py``, as a module."""
    return load_module(ROOT / "website" / "build.py", "website_build")
