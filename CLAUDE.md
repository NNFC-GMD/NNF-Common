# NNF-Common

Notebooks and tutorials for the Novo Nordisk Foundation Center for Genomic
Mechanisms of Disease (NNFC), written as marimo notebooks and published as
WebAssembly pages on GitHub Pages. Readers include people who do not code.

`README.md` is the user-facing reference: how to read, run, add and publish a
notebook, and what each `meta.toml` key means. This file is the
contributor-facing complement. Don't duplicate the README here; point at it.

## Layout

```text
notebooks/<name>/        one folder per notebook
    <name>.py            the marimo notebook; browser-safe imports only
    meta.toml            its card on the landing page
    prepare.py           optional; native; raw data -> public/
    public/              small derived files, published with the page
website/
    build.py             exports every notebook, renders the landing page
    check_page.py        opens one page in headless Chrome
    check_site.py        serves _site/ and checks every page
    templates/, static/  the landing page
tests/                   notebook rules, method functions, site builder
validation/scslide/      Seurat 6 reference run and the comparison with Python
.github/workflows/       ci.yml (checks), pages.yml (build, check, publish)
data/                    raw downloads; ignored by git
_site/                   the built site; ignored by git
.r-lib/                  isolated R library for the reference run; ignored by git
```

## Commands

```bash
pixi run edit notebooks/scslide/scslide.py    # marimo editor
pixi exec marimo edit --sandbox pixi notebooks/scslide/scslide.py   # one file, own environment
pixi run build                                # export everything to _site/
pixi run build --only scslide --no-execute    # one notebook, faster
pixi run serve                                # http://localhost:8000
pixi run -e e2e check-site                    # every page in headless Chrome
pixi run -e dev test
pixi run -e dev lint
pixi run -e dev marimo check --fix notebooks/scslide/scslide.py
pixi run -e prep prepare-scslide              # rebuild notebooks/scslide/public/
Rscript validation/scslide/seurat_reference.R --seed 123
pixi run -e dev python validation/scslide/compare.py
pixi lock                                     # after any dependency change
```

pixi is the one environment manager. `default` holds what a notebook may
import; `dev` adds the test and lint tools; `prep` adds the native packages
that only `prepare.py` scripts use; `e2e` adds playwright.

## Conventions

- **The browser decides what a notebook may import.** The published pages run
  on Pyodide. The default pixi environment is held at the package versions of
  the Pyodide build that the pinned marimo ships (marimo 0.25.1, Pyodide
  v314.0.0), so a notebook that runs here runs the same way there. Moving the
  marimo pin means moving those pins with it: read the versions from
  `https://cdn.jsdelivr.net/pyodide/<version>/full/pyodide-lock.json`.
  scanpy, anndata, numba, umap-learn and POT do not exist in Pyodide.
- **Heavy work goes in `prepare.py`.** It runs natively, downloads the raw
  data into `data/`, and writes a few megabytes to `public/` with a
  `manifest.json` (provenance, parameters, checksums). Its output must not
  depend on when it runs: `write_npz` in the scSLIDE one fixes the time stamps
  that `numpy.savez` would add.
- **`public/` is published to the world.** Only data that is already public
  goes there, never unpublished center data. At most 5 MB per file and 10 MB
  per notebook. Everywhere else the limit is 512 KB.
- **Method code lives in the notebook, as top-level functions**
  (`@app.function`, with imports in `with app.setup:`). Readers can see it,
  the tests import it from the notebook file, and `prepare.py` imports it too,
  so nothing is implemented twice. There is no shared package yet.
- **Every number in the prose is computed or attributed.** Write it as an
  f-string from the data, or say whose result it is. Do not type a result
  into the text. This is the docs-accuracy rule: every concrete detail must
  be confirmable from source, and if it cannot be verified it is left out.
- **Write for someone who has never coded.** Short sentences and plain words.
  Explain a term where it first appears. Each step of a method gets the same
  four parts: the idea, a toy to play with, the real data, and what would
  happen if we did something else.
- **A replication says what it could not replicate.** Report the comparison
  with the paper and with the reference implementation, including where they
  differ and what was not tested.
- **Nothing is taken from a paper's figures or text.** Figures are drawn from
  data, explanations are in our own words, and the DOI is linked.
- **scSLIDE here is an independent implementation.** It was written from the
  paper's Methods and the public documentation, without Seurat's source. That
  independence is what gives the comparison with Seurat its meaning.
  `validation/` runs Seurat as a black box and reads only its outputs.
- **Charts: Okabe-Ito for categories, cividis for continuous scales.** At most
  eight categories get a colour; with more, highlight one group against grey.
  The economist-design skill is not used in this repository.
- **Indentation: 4 spaces, no tabs**; YAML at 2. ruff is configured in
  `pyproject.toml` and pinned through pre-commit.
- **Commit messages: UPPERCASE verb prefix** -- `ADD`, `FIX`, `UPDATE`,
  `REWRITE`, `REMOVE`. `main` changes through pull requests. Commit and push
  only when asked.

## Gotchas

- The site generator is `website/`, not `site/`: `site` is a module of the
  Python standard library.
- `marimo export` drops a file named `CLAUDE.md`, its own assistant prompt,
  into every exported folder. `build.py` deletes it.
- A page exported with `--execute` first shows outputs computed at build time,
  then computes everything again in the browser. A text that is in the
  preview therefore proves nothing about the browser. `browser_text` in
  `meta.toml` must be something only the browser run produces.
- `--execute` runs the notebook in a uv environment that marimo creates, so
  the build needs uv and network access.
- A cell that fails in the browser leaves no trace in a read-only page:
  marimo hides it and everything that depends on it. The error only reaches
  the browser console, as a log message. `check_page.py` reads the console
  for that reason. To see a failure yourself, open the browser's console.
- Do not use `mo.ui.table` for a small table. In the preview it asks a
  notebook that is not running yet for its rows and logs an error. Use
  `markdown_table`, as the scSLIDE notebook does.
- Text that marimo interpolates into an indented `mo.md(f"""...""")` must be
  a single line. A multi-line value breaks the common indentation and the
  whole cell renders as a code block.
- `mo.notebook_location()` is a folder on a computer and a URL in the
  browser. Read files through `read_public`, as both notebooks do.
- The browser runs Python 3.14 and the pixi environment 3.13. Keep to syntax
  that both accept.
- `marimo check --strict` fails on a Markdown cell whose text is indented
  more than one level. `marimo check --fix` rewrites it.
- An exported page is about 26 MB, nearly all of it marimo's own assets.
  `edit.html` shares them, so the editable copy costs almost nothing.
- `edit.html` opens marimo's editor and computes nothing until the reader
  presses the run button. It also logs errors about a language server that a
  static page cannot have. `check_site.py` does not check it.
- scSLIDE specifics:
    - The paper states no settings for the Perturb-seq dataset. Ours are the
      ones it gives for all its analyses, with the target gene as the PLS
      label. `manifest.json` records them.
    - Landmark cells do not vote by default. The paper's Methods count every
      cell, but Seurat's count table sums to five votes per non-landmark cell.
    - Correlation magnitudes depend strongly on the number of landmarks and
      of votes per cell. The pattern does not.
    - Python matches Seurat on the whole guide-by-guide table, with one
      exception: the RPS14/RPS15/RPS18 block, which swings with the balance
      between the two views. Seurat's value is what Python gives at a view
      weight near 0.45. The cause is not known; `compare.py` records the sweep.
    - `diffusion_map` links each sample to its nearest neighbours only. If
      the samples fall into separate groups, its components say which group a
      sample is in and nothing else.
    - The PLS coordinates are baked into `cells.npz`. Changing a PLS choice
      means running `prepare.py` again, then `compare.py`.
- The Seurat reference run uses an isolated R library in `.r-lib/`. Never
  install into the machine's own R library.
