"""Build the site: export every notebook to WebAssembly and render the landing page.

Each folder under ``notebooks/`` holding a ``meta.toml`` is one page. The
notebook it names is exported with ``marimo export html-wasm`` to
``<out>/<folder name>/``, and the landing page lists the pages whose
``listed`` flag is not false.

Run it through pixi: ``pixi run build``.
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
import tomllib
from pathlib import Path

import jinja2

ROOT = Path(__file__).resolve().parent.parent
NOTEBOOKS_DIR = ROOT / "notebooks"
TEMPLATES_DIR = Path(__file__).resolve().parent / "templates"
STATIC_DIR = Path(__file__).resolve().parent / "static"

REQUIRED_KEYS = ("title", "summary", "notebook")


def load_pages(notebooks_dir: Path = NOTEBOOKS_DIR) -> list[dict]:
    """Read every ``meta.toml`` under ``notebooks_dir``.

    Parameters
    ----------
    notebooks_dir
        Folder whose sub-folders each hold one notebook and its ``meta.toml``.

    Returns
    -------
    list of dict
        One entry per page, sorted by ``order`` then title. Each holds the
        keys of its ``meta.toml`` plus ``slug`` (the folder name) and
        ``source`` (the path of the notebook file).

    Raises
    ------
    ValueError
        If a ``meta.toml`` lacks a required key or names a missing notebook.
    """
    pages = []
    for meta_path in sorted(notebooks_dir.glob("*/meta.toml")):
        with meta_path.open("rb") as handle:
            meta = tomllib.load(handle)
        missing = [key for key in REQUIRED_KEYS if key not in meta]
        if missing:
            raise ValueError(f"{meta_path}: missing {', '.join(missing)}")
        source = meta_path.parent / meta["notebook"]
        if not source.is_file():
            raise ValueError(f"{meta_path}: notebook {meta['notebook']!r} not found")
        meta.setdefault("listed", True)
        meta.setdefault("execute", True)
        meta.setdefault("edit_copy", False)
        meta.setdefault("order", 100)
        meta.setdefault("tags", [])
        meta["slug"] = meta_path.parent.name
        meta["source"] = source
        pages.append(meta)
    return sorted(pages, key=lambda page: (page["order"], page["title"]))


def export_page(page: dict, out_dir: Path, execute: bool = True) -> Path:
    """Export one notebook as a read-only WebAssembly page.

    Parameters
    ----------
    page
        An entry from :func:`load_pages`.
    out_dir
        The site folder; the page lands in ``out_dir / page["slug"]``.
    execute
        Run the notebook first and embed its outputs, so a reader sees the
        page before the Python runtime has loaded in the browser.

    Returns
    -------
    pathlib.Path
        The folder the page was written to.

    Raises
    ------
    subprocess.CalledProcessError
        If the export fails.

    Notes
    -----
    A page whose ``edit_copy`` flag is set also gets ``edit.html`` in the same
    folder: the notebook with its code visible and editable. It shares the
    read-only page's assets, so it adds almost nothing to the site's size.
    """
    target = out_dir / page["slug"]
    export = [
        sys.executable,
        "-m",
        "marimo",
        "export",
        "html-wasm",
        str(page["source"]),
    ]
    command = [*export, "-o", str(target), "--mode", "run", "--no-show-code", "--force"]
    if execute and page["execute"]:
        command.append("--execute")
    subprocess.run(command, check=True, cwd=ROOT)
    if page["edit_copy"]:
        command = [
            *export,
            "-o",
            str(target / "edit.html"),
            "--mode",
            "edit",
            "--force",
        ]
        subprocess.run(command, check=True, cwd=ROOT)
    # marimo drops its own assistant prompt into every export; it is not part
    # of the page.
    (target / "CLAUDE.md").unlink(missing_ok=True)
    return target


def render_index(pages: list[dict], out_dir: Path) -> Path:
    """Write the landing page and its stylesheet into ``out_dir``.

    Parameters
    ----------
    pages
        Entries from :func:`load_pages`; those with ``listed`` false are left out.
    out_dir
        The site folder.

    Returns
    -------
    pathlib.Path
        The path of the written ``index.html``.
    """
    environment = jinja2.Environment(
        loader=jinja2.FileSystemLoader(TEMPLATES_DIR),
        autoescape=jinja2.select_autoescape(["html", "j2"]),
        undefined=jinja2.StrictUndefined,
    )
    template = environment.get_template("index.html.j2")
    listed = [page for page in pages if page["listed"]]
    out_dir.mkdir(parents=True, exist_ok=True)
    index_path = out_dir / "index.html"
    index_path.write_text(template.render(pages=listed), encoding="utf-8")
    shutil.copytree(STATIC_DIR, out_dir / "static", dirs_exist_ok=True)
    # GitHub Pages would otherwise run Jekyll, which drops paths starting with "_".
    (out_dir / ".nojekyll").touch()
    return index_path


def main(argv: list[str] | None = None) -> int:
    """Build the site. Returns the process exit code."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--out",
        type=Path,
        default=ROOT / "_site",
        help="folder to write the site to (default: _site)",
    )
    parser.add_argument(
        "--only",
        metavar="SLUG",
        action="append",
        help="export only this notebook folder; repeat for several",
    )
    parser.add_argument(
        "--no-export",
        action="store_true",
        help="render the landing page only",
    )
    parser.add_argument(
        "--no-execute",
        action="store_true",
        help="skip running the notebooks, so pages carry no pre-computed outputs",
    )
    args = parser.parse_args(argv)

    pages = load_pages()
    known = {page["slug"] for page in pages}
    unknown = set(args.only or []) - known
    if unknown:
        parser.error(f"unknown notebook folder: {', '.join(sorted(unknown))}")

    if not args.no_export:
        for page in pages:
            if args.only and page["slug"] not in args.only:
                continue
            print(f"exporting {page['slug']}", flush=True)
            export_page(page, args.out, execute=not args.no_execute)
    index_path = render_index(pages, args.out)
    print(f"wrote {index_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
