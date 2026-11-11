# NNF-Common

Notebooks and tutorials for the Novo Nordisk Foundation Center for Genomic
Mechanisms of Disease (NNFC). Each one explains an analysis in plain words and
lets the reader try it. The notebooks are
[marimo](https://marimo.io) notebooks, published as pages that run in the
reader's browser: there is nothing to install and no server.

## Notebooks

| Notebook | What it covers |
| --- | --- |
| [`notebooks/scslide`](notebooks/scslide/scslide.py) | scSLIDE, a method for comparing samples instead of cells (Jiang et al., Cell, [doi:10.1016/j.cell.2026.10.018](https://doi.org/10.1016/j.cell.2026.10.018)). Each step explained with an example to play with, then replicated in Python on a CRISPR screen and checked against Seurat 6. |
| [`notebooks/hello`](notebooks/hello/hello.py) | The smallest page that exercises the publishing pipeline. Copy it to start a notebook. It is built with the site but not listed on it. |

## Read a notebook

The site is built from this repository and served by GitHub Pages at
<https://nnfc-gmd.github.io/NNF-Common/>. Pages has to be switched on for the
repository first (see [Publishing](#publishing)). The site is public even
though the repository is private.

To read a notebook without the site, build it locally:

```bash
pixi run build      # export every notebook to _site/
pixi run serve      # then open http://localhost:8000
```

## Run a notebook on your computer

Install [pixi](https://pixi.sh), then:

```bash
pixi run edit notebooks/scslide/scslide.py
```

This opens the notebook in marimo's editor, with its code. The first run
downloads the environment.

A single notebook file also runs on its own, without this repository's
environment, because it lists what it needs at its top. `pixi exec` fetches a
marimo recent enough to read that list (0.25 or later):

```bash
pixi exec marimo edit --sandbox pixi notebooks/scslide/scslide.py
```

## Add a notebook

A notebook is a folder under `notebooks/`:

```text
notebooks/<name>/
    <name>.py       the marimo notebook
    meta.toml       its card on the landing page
    public/         small data files the page loads
    prepare.py      optional: builds public/ from raw data
```

`meta.toml` keys:

| Key | Meaning |
| --- | --- |
| `title`, `summary` | Shown on the landing page. Required. |
| `notebook` | File name of the notebook. Required. |
| `audience`, `minutes`, `tags` | Shown on the card. |
| `[paper]` `citation`, `doi` | The paper the notebook is about. |
| `order` | Position on the landing page; lower comes first. Default 100. |
| `listed` | `false` builds the page but keeps it off the landing page. |
| `execute` | `false` skips running the notebook at build time, so the page shows nothing until Python has started in the browser. Default `true`. |
| `edit_copy` | `true` also publishes `edit.html`, the page with its code visible and editable. |
| `browser_text` | A text the page shows only once Python has run in the browser. The site check waits for it. |

Three rules, each enforced by a test:

- **A notebook imports only what runs in a browser.** That is the list under
  `[tool.pixi.dependencies]` in `pyproject.toml`. Anything heavier goes in a
  `prepare.py`, which runs on a computer and writes small files to `public/`.
- **`public/` is published to the world.** Put only data there that is already
  public. A file may be at most 5 MB, a notebook's `public/` folder 10 MB.
- **The top of the notebook lists its packages twice**: once for the browser
  (`dependencies`) and once for marimo's pixi sandbox
  (`[tool.pixi.dependencies]`), with the versions of `pyproject.toml`.

## Check the site before it is published

```bash
pixi run build
pixi run -e e2e check-site
```

`check-site` opens every page in headless Google Chrome and fails on a page
error, a console error, a failed request, a failed notebook cell, or a page
that never finishes computing. It serves the site under a sub-path, as GitHub
Pages does, so a link that only works at the root of a server fails too.

## Publishing

`.github/workflows/pages.yml` builds the site, checks it and publishes it when
a push to `main` changes a notebook, the site builder or the environment. It
can also be started by hand from the Actions tab. One setting has to be made
once, by a repository administrator: Settings, Pages, Build and deployment,
Source: **GitHub Actions**.

## Development

```bash
pixi run -e dev test             # the tests
pixi run -e dev lint             # every pre-commit hook on every tracked file
pixi run -e dev install-hooks    # once per clone
```

The scSLIDE notebook's data is rebuilt, from a 336 MB download, with:

```bash
pixi run -e prep prepare-scslide
```

## License

MIT. See `LICENSE`.

The scSLIDE code was written from the paper's Methods. No Seurat source code
was copied or translated. Figures are drawn from public data; none is taken
from a paper.
