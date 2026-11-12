"""Rules every published notebook must meet."""

from __future__ import annotations

import hashlib
import json
import re
import subprocess
import sys
import tomllib

import pytest
from conftest import ROOT, notebook_pages

KB = 1024
MAX_PUBLIC_FILE = 5 * KB * KB
MAX_PUBLIC_TOTAL = 10 * KB * KB

# Files in public/ that prepare.py does not write, so its manifest does not list them.
UNLISTED_PUBLIC_FILES = {"manifest.json", "reading_guide.md", "validation.json"}

# The inline-metadata block defined by PEP 723.
SCRIPT_BLOCK = re.compile(
    r"(?m)^# /// script$\s(?P<content>(^#(| .*)$\s)+)^# ///$",
)

PAGES = notebook_pages()
PAGE_IDS = [page.name for page in PAGES]


def read_meta(page):
    with (page / "meta.toml").open("rb") as handle:
        return tomllib.load(handle)


def inline_metadata(notebook_path):
    """Return a notebook's PEP 723 script block, parsed."""
    match = SCRIPT_BLOCK.search(notebook_path.read_text(encoding="utf-8"))
    assert match, f"{notebook_path} has no PEP 723 script block"
    content = "".join(
        line[2:] if line.startswith("# ") else line[1:]
        for line in match.group("content").splitlines(keepends=True)
    )
    return tomllib.loads(content)


def inline_dependencies(notebook_path):
    """Return the package names in a notebook's PEP 723 ``dependencies`` list."""
    dependencies = inline_metadata(notebook_path)["dependencies"]
    return {
        re.split(r"[<>=!~ \[;]", item, maxsplit=1)[0].lower() for item in dependencies
    }


def workspace_dependencies():
    """Return the default pixi dependencies: what a notebook may import, and at which version."""
    with (ROOT / "pyproject.toml").open("rb") as handle:
        return tomllib.load(handle)["tool"]["pixi"]["dependencies"]


def browser_packages():
    """Return the names of the packages a notebook may import."""
    return {name.lower() for name in workspace_dependencies()}


@pytest.mark.parametrize("page", PAGES, ids=PAGE_IDS)
def test_marimo_check_passes(page):
    notebook = page / read_meta(page)["notebook"]
    result = subprocess.run(
        [sys.executable, "-m", "marimo", "check", "--strict", str(notebook)],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stdout + result.stderr


@pytest.mark.parametrize("page", PAGES, ids=PAGE_IDS)
def test_notebook_runs_as_a_script(page):
    notebook = page / read_meta(page)["notebook"]
    result = subprocess.run(
        [sys.executable, str(notebook)], capture_output=True, text=True, cwd=page
    )
    assert result.returncode == 0, result.stdout + result.stderr


@pytest.mark.parametrize("page", PAGES, ids=PAGE_IDS)
def test_inline_dependencies_are_browser_packages(page):
    notebook = page / read_meta(page)["notebook"]
    extra = inline_dependencies(notebook) - browser_packages()
    assert not extra, (
        f"{notebook.name} declares {sorted(extra)}, which the default pixi "
        "environment does not hold. Only packages that exist in Pyodide belong there."
    )


@pytest.mark.parametrize("page", PAGES, ids=PAGE_IDS)
def test_inline_pixi_pins_match_the_workspace(page):
    """`marimo edit --sandbox pixi` must build the same packages as `pixi run`."""
    notebook = page / read_meta(page)["notebook"]
    inline = inline_metadata(notebook)["tool"]["pixi"]["dependencies"]
    workspace = workspace_dependencies()
    differing = {
        name: spec for name, spec in inline.items() if workspace.get(name) != spec
    }
    assert not differing, (
        f"{notebook.name} pins differ from pyproject.toml: {differing}"
    )
    missing = inline_dependencies(notebook) - {name.lower() for name in inline}
    assert not missing, f"{notebook.name} has no pixi pin for {sorted(missing)}"


@pytest.mark.parametrize("page", PAGES, ids=PAGE_IDS)
def test_public_files_stay_small(page):
    files = [path for path in (page / "public").rglob("*") if path.is_file()]
    too_big = [path.name for path in files if path.stat().st_size > MAX_PUBLIC_FILE]
    assert not too_big, f"over 5 MB: {too_big}"
    assert sum(path.stat().st_size for path in files) <= MAX_PUBLIC_TOTAL


@pytest.mark.parametrize("page", PAGES, ids=PAGE_IDS)
def test_public_files_match_their_manifest(page):
    manifest_path = page / "public" / "manifest.json"
    if not manifest_path.is_file():
        pytest.skip("no manifest.json")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    listed = {entry["file"]: entry["sha256"] for entry in manifest["files"]}
    present = {
        path.name
        for path in (page / "public").iterdir()
        if path.is_file() and path.name not in UNLISTED_PUBLIC_FILES
    }
    assert present == set(listed)
    for name, expected in listed.items():
        digest = hashlib.sha256((page / "public" / name).read_bytes()).hexdigest()
        assert digest == expected, f"{name} differs from its manifest entry"


def test_scslide_landmarks_in_the_file_are_the_ones_the_notebook_draws():
    """The default landmark cells are stored for reference; the notebook must draw the same."""
    import numpy as np
    from conftest import NOTEBOOKS_DIR, load_module

    scslide = load_module(NOTEBOOKS_DIR / "scslide" / "scslide.py", "notebook_scslide")
    cells = np.load(NOTEBOOKS_DIR / "scslide" / "public" / "cells.npz")

    drawn = scslide.sample_landmarks(
        cells["leverage"], int(cells["landmark"].sum()), seed=0
    )

    assert np.array_equal(drawn, np.flatnonzero(cells["landmark"]))


def test_scslide_surroundings_in_the_file_are_the_ones_the_notebook_measures():
    """Every cell's surroundings are stored by prepare.py; the notebook function must agree."""
    import numpy as np
    from conftest import NOTEBOOKS_DIR, load_module

    scslide = load_module(NOTEBOOKS_DIR / "scslide" / "scslide.py", "notebook_scslide")
    cells = np.load(NOTEBOOKS_DIR / "scslide" / "public" / "cells.npz")
    manifest = json.loads(
        (NOTEBOOKS_DIR / "scslide" / "public" / "manifest.json").read_text(
            encoding="utf-8"
        )
    )

    offsets, widths = scslide.local_scales(
        [cells["pca"].astype(np.float32), cells["pls"].astype(np.float32)],
        rank=manifest["parameters"]["bandwidth_rank"],
    )

    assert np.allclose(offsets, cells["scale_offset"], atol=1e-5)
    assert np.allclose(widths, cells["scale_width"], atol=1e-5)
