"""Compare the Python scSLIDE results with the Seurat 6 reference run.

The comparison is made on the 138 x 138 matrix of correlations between guide
RNAs, the end product of the Perturb-seq analysis, and cell by cell on the
five landmarks that each cell votes for.

Four questions, in order:

1. How much do two runs of the *same* implementation differ, when only the
   random draw of landmark cells changes? That is the noise floor.
2. How close is Python to Seurat, judged against that floor?
3. Does a difference sit in the cell coordinates or in the later steps?
   Seurat's own cell coordinates and landmark cells are fed through the
   Python sample-level functions.
4. Which rule gives every cell the landmarks Seurat gives it? The paper's
   Methods and Seurat's code describe a cell's surroundings differently (see
   ``local_scales`` in the notebook). Both rules are run on Seurat's own
   coordinates and landmark cells. So is the code's rule with the numbers
   Seurat itself computed, which differ from the exact ones because Seurat
   searches for neighbours approximately.

The first comparison treated Seurat as a black box and left one difference
unexplained. Seurat's source (MIT licence) was then read to explain it;
question 4 is the result.

Run it after ``prepare.py`` and at least two reference runs:

    pixi run -e dev python validation/scslide/compare.py
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
NOTEBOOK_DIR = ROOT / "notebooks" / "scslide"

# Two distances closer than this are taken to be the same distance.
SAME_DISTANCE = 1e-5


def load_notebook_functions():
    """Import the notebook file, to use the method functions defined in it."""
    spec = importlib.util.spec_from_file_location(
        "scslide_notebook", NOTEBOOK_DIR / "scslide.py"
    )
    module = importlib.util.module_from_spec(spec)
    sys.modules["scslide_notebook"] = module
    spec.loader.exec_module(module)
    return module


def sample_level(
    methods,
    views,
    landmark_rows,
    guide,
    n_guides,
    scales=None,
    count_landmark_cells=False,
):
    """Run the sample-level steps of the notebook.

    Parameters
    ----------
    methods
        The notebook module.
    views : list of arrays of shape (n_cells, n_dims)
        The two cell embeddings.
    landmark_rows : int array
        Rows of the landmark cells. A landmark is numbered by its place here.
    guide : int array of shape (n_cells,)
        Guide of every cell.
    n_guides : int
    scales : pair of arrays of shape (2, n_cells), optional
        Every cell's surroundings, as ``local_scales`` returns them: the rule
        of Seurat's code. When None the surroundings are read off the
        landmarks: the rule of the paper's Methods.
    count_landmark_cells : bool
        Whether landmark cells vote. In Seurat they do not; the paper's
        Methods count every cell.

    Returns
    -------
    dict
        What ``run_sample_level`` returns: ``landmarks``, ``counts``,
        ``residuals`` and ``correlation``.
    """
    pool, affinities, repeated = methods.landmark_affinities(
        views, [view[landmark_rows] for view in views], scales=scales
    )
    return methods.run_sample_level(
        pool,
        affinities,
        repeated,
        guide,
        n_guides,
        landmark_rows,
        landmarks_vote=count_landmark_cells,
    )


def describe(methods, correlation, guides):
    """Summarise one correlation matrix with the numbers the notebook reports."""
    gene = guides["gene"].to_numpy()
    targeting = ~guides["control"].to_numpy()
    controls = np.flatnonzero(~targeting)
    activity = guides["relative_activity_day10"].to_numpy()
    predicted = methods.nearest_label_votes(1.0 - correlation, gene, k=3)
    near_controls = correlation[:, controls].mean(axis=1)
    active = targeting & (activity >= 0.5)

    def block(first, second):
        rows = np.flatnonzero(active & (gene == first))
        columns = np.flatnonzero(active & (gene == second))
        values = correlation[np.ix_(rows, columns)]
        if first == second:
            values = values[np.triu_indices(len(rows), 1)]
        return float(values.mean())

    ribosomal = [("RPS14", "RPS15"), ("RPS14", "RPS18"), ("RPS15", "RPS18")]
    return {
        "3-NN accuracy": float((predicted[targeting] == gene[targeting]).mean()),
        "Spearman, near controls vs activity": float(
            spearmanr(near_controls[targeting], activity[targeting]).statistic
        ),
        "control-control r": float(
            correlation[np.ix_(controls, controls)][
                np.triu_indices(len(controls), 1)
            ].mean()
        ),
        "same-gene r (active guides)": float(
            np.mean([block(name, name) for name in sorted(set(gene[targeting]))])
        ),
        "RPS14/15/18 cross r": float(np.mean([block(*pair) for pair in ribosomal])),
        "EIF2S1-HSPA9 r": block("EIF2S1", "HSPA9"),
    }


def upper_triangle_r(first, second):
    """Pearson correlation between the upper triangles of two square matrices."""
    upper = np.triu_indices(len(first), 1)
    return float(np.corrcoef(first[upper], second[upper])[0, 1])


def subspace_agreement(first, second):
    """Mean canonical correlation between two embeddings of the same cells.

    One means the two span the same directions, however rotated or scaled.
    """
    first_basis = np.linalg.svd(first - first.mean(axis=0), full_matrices=False)[0]
    second_basis = np.linalg.svd(second - second.mean(axis=0), full_matrices=False)[0]
    return float(np.linalg.svd(first_basis.T @ second_basis, compute_uv=False).mean())


def same_landmarks(ours, theirs):
    """Fraction of rows that hold the same landmarks, in any order."""
    return float((np.sort(ours, axis=1) == np.sort(theirs, axis=1)).all(axis=1).mean())


def read_reference(reference_dir, seed, kept_barcodes, guide_order):
    """Read what one reference run wrote, with its cells in the Python order.

    Parameters
    ----------
    reference_dir : pathlib.Path
        Folder holding the outputs of ``seurat_reference.R``.
    seed : int
    kept_barcodes : list of str
        Cell barcodes in the order of ``cells.npz``.
    guide_order : sequence of str
        Seurat's guide identifiers, in the order of the Python guide table.

    Returns
    -------
    dict
        ``cells`` (table), ``views`` (two arrays), ``landmark_rows`` (row of
        Seurat's landmark number 0, 1, ...), ``neighbours`` (n_cells, 20;
        Seurat's landmark numbers, best first), ``counts`` (landmark by
        guide), ``scales`` (offsets and widths, each (2, n_cells)) and
        ``correlation`` (guide by guide).
    """
    tag = f"seed{seed}"
    cells = (
        pd.read_csv(reference_dir / f"cells_{tag}.csv.gz")
        .set_index("cell_barcode")
        .loc[kept_barcodes]
    )
    views = [
        cells[[name for name in cells.columns if name.startswith(prefix)]].to_numpy(
            np.float32
        )
        for prefix in ("PC_", "pPLS_")
    ]
    row_of = pd.Series(np.arange(len(cells)), index=cells.index)
    landmark_order = (reference_dir / f"landmark_order_{tag}.txt").read_text().split()
    landmark_rows = row_of.loc[landmark_order].to_numpy()

    neighbours = (
        pd.read_csv(reference_dir / f"neighbours_{tag}.csv.gz")
        .set_index("cell_barcode")
        .loc[kept_barcodes]
    )
    index_columns = [name for name in neighbours.columns if name.startswith("index")]
    scales = (
        pd.read_csv(reference_dir / f"scales_{tag}.csv.gz")
        .set_index("cell_barcode")
        .loc[kept_barcodes]
    )
    return {
        "cells": cells,
        "views": views,
        "landmark_rows": landmark_rows,
        # R counts from one.
        "neighbours": neighbours[index_columns].to_numpy(np.int64) - 1,
        "counts": pd.read_csv(
            reference_dir / f"landmark_counts_{tag}.csv.gz", index_col=0
        )[list(guide_order)].to_numpy(),
        "scales": (
            scales[["offset_pca", "offset_pls"]].to_numpy(np.float32).T,
            scales[["width_pca", "width_pls"]].to_numpy(np.float32).T,
        ),
        "correlation": pd.read_csv(
            reference_dir / f"correlation_{tag}.csv", index_col=0
        )
        .loc[guide_order, guide_order]
        .to_numpy(),
    }


def main(argv: list[str] | None = None) -> int:
    """Print the comparison and save its numbers. Returns the process exit code."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--reference-dir",
        type=Path,
        default=HERE / "out",
        help="folder holding the outputs of seurat_reference.R",
    )
    parser.add_argument(
        "--save",
        type=Path,
        default=NOTEBOOK_DIR / "public" / "validation.json",
        help="where to save the numbers the notebook shows",
    )
    args = parser.parse_args(argv)
    methods = load_notebook_functions()

    data = np.load(NOTEBOOK_DIR / "public" / "cells.npz")
    guides = pd.read_csv(NOTEBOOK_DIR / "public" / "sgrna.csv")
    views = [data["pca"].astype(np.float32), data["pls"].astype(np.float32)]
    guide = data["guide"].astype(int)
    leverage = data["leverage"]
    n_landmarks = int(data["landmark"].sum())
    n_guides = len(guides)
    # The notebook's default: every cell's surroundings among all cells, as in
    # Seurat's code, and no votes from landmark cells, again as in Seurat.
    scales = (data["scale_offset"], data["scale_width"])

    # Seurat's guides, reordered to the Python guide order.
    reference_guides = pd.read_csv(args.reference_dir / "guides.csv")
    guide_order = (
        reference_guides.set_index("guide_identity")
        .loc[guides["guide_identity"], "id"]
        .to_numpy()
    )
    seeds = sorted(
        int(path.stem.removeprefix("correlation_seed"))
        for path in args.reference_dir.glob("correlation_seed*.csv")
    )
    if len(seeds) < 2:
        parser.error(
            f"need two reference runs in {args.reference_dir}, found {len(seeds)}"
        )
    kept_barcodes = (
        (ROOT / "data" / "GSE132080" / "kept_barcodes.txt").read_text().split()
    )
    reference = {
        seed: read_reference(args.reference_dir, seed, kept_barcodes, guide_order)
        for seed in seeds
    }
    seurat = {seed: run["correlation"] for seed, run in reference.items()}
    summary = dict(
        line.split(" ", 1)
        for line in (args.reference_dir / f"summary_seed{seeds[0]}.txt")
        .read_text()
        .splitlines()
    )

    # --- 1 and 2: the end product, at different draws of landmarks ----------
    python = {}
    paper_text = {}
    for seed in (0, 1, 2):
        rows = methods.sample_landmarks(leverage, n_landmarks, seed=seed)
        python[seed] = sample_level(
            methods, views, rows, guide, n_guides, scales=scales
        )["correlation"]
        paper_text[seed] = sample_level(methods, views, rows, guide, n_guides)[
            "correlation"
        ]

    # Does it matter that landmarks are drawn by leverage score rather than at random?
    uniform_rows = np.sort(
        np.random.default_rng(0).choice(len(guide), n_landmarks, replace=False)
    )
    uniform = sample_level(
        methods, views, uniform_rows, guide, n_guides, scales=scales
    )["correlation"]
    uniform_vs_leverage = upper_triangle_r(uniform, python[0])

    def pairwise(runs):
        keys = sorted(runs)
        return [
            upper_triangle_r(runs[a], runs[b])
            for index, a in enumerate(keys)
            for b in keys[index + 1 :]
        ]

    def against_seurat(runs):
        return [
            upper_triangle_r(run, other)
            for other in seurat.values()
            for run in runs.values()
        ]

    def spread(values):
        return {
            "mean": round(float(np.mean(values)), 4),
            "min": round(min(values), 4),
            "max": round(max(values), 4),
        }

    python_floor = pairwise(python)
    seurat_floor = pairwise(seurat)
    across = against_seurat(python)
    across_text = against_seurat(paper_text)
    print("1. Noise floor: the same implementation, a different draw of landmarks")
    print("   Python vs Python:", ", ".join(f"{value:.3f}" for value in python_floor))
    print("   Seurat vs Seurat:", ", ".join(f"{value:.3f}" for value in seurat_floor))
    print("\n2. Python against Seurat (upper triangle of the 138 x 138 matrix)")
    print(
        f"   r = {np.mean(across):.3f} (range {min(across):.3f} to {max(across):.3f} "
        f"over {len(across)} pairs of runs)"
    )
    print(
        f"   with the rule of the paper's text: r = {np.mean(across_text):.3f} "
        f"(range {min(across_text):.3f} to {max(across_text):.3f})"
    )
    print(
        f"   Landmarks drawn at random instead of by leverage: r = {uniform_vs_leverage:.3f}"
    )

    table = {
        "Python": describe(methods, python[0], guides),
        "Python, following the paper's text": describe(methods, paper_text[0], guides),
    }
    table.update(
        {
            f"Seurat, seed {seed}": describe(methods, run, guides)
            for seed, run in seurat.items()
        }
    )

    # --- 3: coordinates or later steps? --------------------------------------
    print("\n3. Seurat's landmark cells in the Python pipeline")
    first = reference[seeds[0]]
    reference_guide = (
        reference_guides.set_index("id")
        .loc[first["cells"]["guide"], "guide_identity"]
        .map(guides.set_index("guide_identity")["guide"])
        .to_numpy()
    )
    if not np.array_equal(reference_guide, guide):
        raise ValueError("the reference cells do not carry the guides of cells.npz")
    agreement = {
        "pca": subspace_agreement(views[0], first["views"][0]),
        "pls": subspace_agreement(views[1], first["views"][1]),
    }
    print(
        "   Same directions? mean canonical correlation, Python vs Seurat: "
        f"PCA {agreement['pca']:.3f}, PLS {agreement['pls']:.3f}"
    )
    transplants = {"python_coordinates": [], "python_coordinates_every_cell_votes": []}
    for seed, run in reference.items():
        result = sample_level(
            methods, views, run["landmark_rows"], guide, n_guides, scales=scales
        )["correlation"]
        transplants["python_coordinates"].append(upper_triangle_r(result, seurat[seed]))
        if seed == seeds[0]:
            # The Python pipeline, changed in one thing only: which cells are landmarks.
            table["Python, with Seurat's landmark cells"] = describe(
                methods, result, guides
            )
        result = sample_level(
            methods,
            views,
            run["landmark_rows"],
            guide,
            n_guides,
            scales=scales,
            count_landmark_cells=True,
        )["correlation"]
        transplants["python_coordinates_every_cell_votes"].append(
            upper_triangle_r(result, seurat[seed])
        )
    transplants = {name: float(np.mean(values)) for name, values in transplants.items()}
    for name, value in transplants.items():
        print(f"   {name.replace('_', ' ')}: r = {value:.4f}")
    print("\n   What each matrix says:")
    print(pd.DataFrame(table).round(3).to_string())

    # --- 4: which rule gives Seurat's landmarks, cell by cell? --------------
    print("\n4. Cell by cell, on Seurat's own coordinates and landmark cells")
    rules = {
        "paper_text": "the paper's text: nearest landmarks",
        "code_exact": "Seurat's code: nearest cells, exact search",
        "code_seurat_values": "Seurat's code, with Seurat's own numbers",
    }
    by_rule = {name: [] for name in rules}
    table_follows = []
    search = {"width_larger": [], "width_same": [], "offset_same": [], "ratio": []}
    for seed, run in reference.items():
        voters = np.ones(len(guide), dtype=bool)
        voters[run["landmark_rows"]] = False
        theirs = run["neighbours"][voters][:, :5]
        # Is Seurat's count table what its own stored neighbours imply, once
        # the landmark cells are left out?
        implied = methods.landmark_counts(
            theirs, guide[voters], len(run["landmark_rows"]), n_guides
        )
        table_follows.append(bool(np.array_equal(implied, run["counts"])))

        exact = methods.local_scales(run["views"])
        for name, chosen in (
            ("paper_text", None),
            ("code_exact", exact),
            ("code_seurat_values", run["scales"]),
        ):
            result = sample_level(
                methods,
                run["views"],
                run["landmark_rows"],
                guide,
                n_guides,
                scales=chosen,
            )
            ours = result["landmarks"][voters]
            by_rule[name].append(
                {
                    "same_five": same_landmarks(ours, theirs),
                    "cells_differing": int(
                        (np.sort(ours, axis=1) != np.sort(theirs, axis=1))
                        .any(axis=1)
                        .sum()
                    ),
                    "votes_moved": int(
                        np.abs(result["counts"] - run["counts"]).sum() // 2
                    ),
                    "r_with_seurat": upper_triangle_r(
                        result["correlation"], seurat[seed]
                    ),
                    "ribosomal": describe(methods, result["correlation"], guides)[
                        "RPS14/15/18 cross r"
                    ],
                }
            )
        # What does the approximate search do to the widths and offsets?
        width_gap = run["scales"][1] - exact[1]
        offset_gap = run["scales"][0] - exact[0]
        search["width_larger"].append(float((width_gap > SAME_DISTANCE).mean()))
        search["width_same"].append(float((np.abs(width_gap) <= SAME_DISTANCE).mean()))
        search["offset_same"].append(
            float((np.abs(offset_gap) <= SAME_DISTANCE).mean())
        )
        search["ratio"].append(
            float(np.median((run["scales"][1] / exact[1])[width_gap > SAME_DISTANCE]))
        )
    rule_rows = []
    for name, label in rules.items():
        row = {
            "rule": name,
            "same_five": float(np.mean([run["same_five"] for run in by_rule[name]])),
            "cells_differing": [run["cells_differing"] for run in by_rule[name]],
            "votes_moved": float(
                np.mean([run["votes_moved"] for run in by_rule[name]])
            ),
            "r_with_seurat": float(
                np.mean([run["r_with_seurat"] for run in by_rule[name]])
            ),
            "ribosomal": float(np.mean([run["ribosomal"] for run in by_rule[name]])),
        }
        rule_rows.append(row)
        print(
            f"   {label:<44} same five landmarks {row['same_five']:.4f} | "
            f"votes placed differently {row['votes_moved']:>7.0f} | "
            f"r with Seurat {row['r_with_seurat']:.5f} | RPS14/15/18 {row['ribosomal']:.3f}"
        )
    seurat_ribosomal = float(
        np.mean(
            [
                describe(methods, run, guides)["RPS14/15/18 cross r"]
                for run in seurat.values()
            ]
        )
    )
    print(f"   {'Seurat itself':<44} RPS14/15/18 {seurat_ribosomal:.3f}")
    print(
        "   Seurat's count table follows from its stored neighbours, landmark "
        f"cells left out: {all(table_follows)}"
    )
    print(
        "   Seurat's approximate search: width larger than the exact one for "
        f"{np.mean(search['width_larger']):.3f} of cells and views (by a median factor "
        f"of {np.mean(search['ratio']):.3f}), the same for {np.mean(search['width_same']):.3f}; "
        f"offset the same for {np.mean(search['offset_same']):.3f}"
    )
    # How much narrower is the code's kernel than the text's, on the notebook's data?
    default_rows = methods.sample_landmarks(leverage, n_landmarks, seed=0)
    width_ratio = []
    for view, width in zip(views, scales[1], strict=True):
        rows = methods.l2_normalize(view)
        _, distances = methods.nearest_rows(rows, rows[default_rows], 20)
        width_ratio.append(width / (distances[:, 19] - distances[:, 0]))
    width_ratio = float(np.median(np.concatenate(width_ratio)))
    print(
        f"   The code's width is {width_ratio:.2f} of the text's (median over cells and views)"
    )

    saved = {
        "seurat_version": summary["seurat"],
        "seurat_seeds": seeds,
        "python_seeds": sorted(python),
        "landmarks": n_landmarks,
        "seurat_count_total": int(summary["landmark_count_total"]),
        "cells": len(guide),
        "python_vs_python": [round(value, 4) for value in python_floor],
        "seurat_vs_seurat": [round(value, 4) for value in seurat_floor],
        "python_vs_seurat": spread(across),
        "paper_text_vs_seurat": spread(across_text),
        "uniform_vs_leverage_landmarks": round(uniform_vs_leverage, 4),
        "subspace_agreement": {
            name: round(value, 4) for name, value in agreement.items()
        },
        "with_seurat_landmark_cells": {
            name: round(value, 4) for name, value in transplants.items()
        },
        "voting_cells": int(len(guide) - n_landmarks),
        "seurat_table_follows_from_its_neighbours": all(table_follows),
        "rules_on_seurat_coordinates": [
            {
                "rule": row["rule"],
                "same_five": round(row["same_five"], 4),
                "cells_differing": row["cells_differing"],
                "votes_moved": round(row["votes_moved"], 1),
                "r_with_seurat": round(row["r_with_seurat"], 4),
                "ribosomal": round(row["ribosomal"], 4),
            }
            for row in rule_rows
        ],
        "seurat_ribosomal": round(seurat_ribosomal, 4),
        "approximate_search": {
            "width_larger": round(float(np.mean(search["width_larger"])), 4),
            "width_same": round(float(np.mean(search["width_same"])), 4),
            "offset_same": round(float(np.mean(search["offset_same"])), 4),
            "median_factor_where_larger": round(float(np.mean(search["ratio"])), 4),
        },
        "code_width_over_text_width": round(width_ratio, 4),
        "what_each_matrix_says": {
            name: {key: round(value, 4) for key, value in row.items()}
            for name, row in table.items()
        },
    }
    args.save.write_text(json.dumps(saved, indent=4) + "\n", encoding="utf-8")
    print(f"\nsaved {args.save}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
