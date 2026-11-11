"""The site builder reads ``meta.toml`` files and renders the landing page."""

from __future__ import annotations

import pytest


def write_page(folder, meta, notebook="nb.py"):
    folder.mkdir(parents=True)
    (folder / "meta.toml").write_text(meta, encoding="utf-8")
    if notebook:
        (folder / notebook).write_text("", encoding="utf-8")


def test_load_pages_fills_defaults_and_sorts(website_build, tmp_path):
    write_page(
        tmp_path / "zebra", 'title = "Zebra"\nsummary = "z"\nnotebook = "nb.py"\n'
    )
    write_page(
        tmp_path / "first",
        'title = "Later title"\nsummary = "f"\nnotebook = "nb.py"\norder = 1\n',
    )

    pages = website_build.load_pages(tmp_path)

    assert [page["slug"] for page in pages] == ["first", "zebra"]
    assert pages[1]["listed"] is True
    assert pages[1]["execute"] is True
    assert pages[1]["edit_copy"] is False
    assert pages[1]["source"] == tmp_path / "zebra" / "nb.py"


@pytest.mark.parametrize(
    ("meta", "notebook", "message"),
    [
        ('title = "T"\nnotebook = "nb.py"\n', "nb.py", "missing summary"),
        ('title = "T"\nsummary = "s"\nnotebook = "gone.py"\n', "nb.py", "not found"),
    ],
)
def test_load_pages_rejects_a_broken_meta(
    website_build, tmp_path, meta, notebook, message
):
    write_page(tmp_path / "page", meta, notebook)

    with pytest.raises(ValueError, match=message):
        website_build.load_pages(tmp_path)


def test_render_index_lists_only_listed_pages(website_build, tmp_path):
    pages = [
        {
            "slug": "shown",
            "title": "Shown page",
            "summary": "s",
            "listed": True,
            "tags": [],
            "edit_copy": True,
        },
        {
            "slug": "hidden",
            "title": "Hidden page",
            "summary": "s",
            "listed": False,
            "tags": [],
            "edit_copy": False,
        },
    ]

    index_path = website_build.render_index(pages, tmp_path / "out")

    html = index_path.read_text(encoding="utf-8")
    assert 'href="shown/"' in html
    assert 'href="shown/edit.html"' in html
    assert "Hidden page" not in html
    assert (tmp_path / "out" / ".nojekyll").is_file()
    assert (tmp_path / "out" / "static" / "site.css").is_file()
