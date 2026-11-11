"""The hello notebook's functions can be imported and used outside marimo."""

from __future__ import annotations

import numpy as np
from conftest import NOTEBOOKS_DIR, load_module


def test_group_means_returns_one_row_per_group():
    hello = load_module(NOTEBOOKS_DIR / "hello" / "hello.py", "notebook_hello")
    points = np.array([[0.0, 0.0], [2.0, 4.0], [10.0, 10.0]])
    groups = np.array(["a", "a", "b"])

    means = hello.group_means(points, groups).set_index("group")

    assert list(means.index) == ["a", "b"]
    assert means.loc["a", ["x", "y"]].tolist() == [1.0, 2.0]
    assert means.loc["b", ["x", "y"]].tolist() == [10.0, 10.0]
