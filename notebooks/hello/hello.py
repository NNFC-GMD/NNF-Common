# /// script
# requires-python = ">=3.12"
# dependencies = [
#     "altair",
#     "marimo",
#     "numpy",
#     "pandas",
# ]
#
# [tool.pixi.workspace]
# channels = ["conda-forge"]
# platforms = ["linux-64", "osx-64", "osx-arm64"]
#
# [tool.pixi.dependencies]
# python = "3.13.*"
# marimo = "==0.25.1"
# numpy = "2.4.*"
# pandas = "3.0.*"
# altair = "6.0.*"
# ///
"""The smallest page that exercises the whole publishing pipeline.

It reads a binary file from ``public/``, draws a chart, reacts to a slider and
reports which Python runs it. It is built with the site but not listed on the
landing page.
"""

import marimo

__generated_with = "0.25.1"
app = marimo.App(width="medium", app_title="Hello, WebAssembly")

with app.setup:
    import io
    import platform
    import sys
    import urllib.request

    import altair as alt
    import marimo as mo
    import numpy as np
    import pandas as pd

    # Okabe-Ito, the colour-blind-safe categorical palette used across the site.
    OKABE_ITO = [
        "#E69F00",
        "#56B4E9",
        "#009E73",
        "#F0E442",
        "#0072B2",
        "#D55E00",
        "#CC79A7",
        "#000000",
    ]


@app.function
def read_public(name):
    """Return the bytes of a file in this notebook's ``public/`` folder.

    On a laptop the folder is a directory next to the notebook. In the
    browser it is a URL next to the page, so the file is downloaded.
    """
    location = mo.notebook_location() / "public" / name
    if str(location).startswith(("http://", "https://")):
        with urllib.request.urlopen(str(location)) as response:
            return response.read()
    return location.read_bytes()


@app.function
def group_means(points, groups):
    """Return the mean position of each group as a data frame."""
    frame = pd.DataFrame(points, columns=["x", "y"]).assign(group=groups)
    return frame.groupby("group", as_index=False)[["x", "y"]].mean()


@app.cell(hide_code=True)
def _():
    mo.md("""
    # Hello, WebAssembly

    This page is a marimo notebook running in your browser. Nothing is
    installed on your computer and no server computes the results.
    """)
    return


@app.cell
def _():
    demo = np.load(io.BytesIO(read_public("demo.npz")))
    points = demo["points"].astype(np.float32)
    groups = demo["groups"]
    mo.md(
        f"Read **{len(points)} points** in **{len(set(groups))} groups** from `public/demo.npz`."
    )
    return groups, points


@app.cell
def _():
    spread = mo.ui.slider(
        0.2, 3.0, value=1.0, step=0.1, label="Spread", show_value=True
    )
    spread
    return (spread,)


@app.cell
def _(groups, points, spread):
    centres = group_means(points, groups).set_index("group")
    offsets = points - centres.loc[groups, ["x", "y"]].to_numpy()
    moved = centres.loc[groups, ["x", "y"]].to_numpy() + spread.value * offsets
    frame = pd.DataFrame(moved, columns=["x", "y"]).assign(group=groups)
    chart = (
        alt.Chart(frame)
        .mark_circle(size=40, opacity=0.7)
        .encode(
            x=alt.X("x:Q", title="x"),
            y=alt.Y("y:Q", title="y"),
            color=alt.Color(
                "group:N",
                title="Group",
                scale=alt.Scale(range=OKABE_ITO[: len(set(groups))]),
            ),
            tooltip=["group", "x", "y"],
        )
        .properties(width=480, height=320)
    )
    mo.ui.altair_chart(chart)
    return


@app.cell(hide_code=True)
def _():
    mo.md(
        f"""
        | | |
        |---|---|
        | Runs on | `{sys.platform}` |
        | Python | {platform.python_version()} |
        | marimo | {mo.__version__} |
        | numpy | {np.__version__} |
        | pandas | {pd.__version__} |
        | altair | {alt.__version__} |
        """
    )
    return


if __name__ == "__main__":
    app.run()
