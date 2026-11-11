"""Compare the Python scSLIDE results with the Seurat 6 reference run.

Seurat is treated as a black box: only what ``seurat_reference.R`` writes is
read. The comparison is made on the 138 x 138 matrix of correlations between
guide RNAs, the end product of the Perturb-seq analysis.

Three questions, in order:

1. How much do two runs of the *same* implementation differ, when only the
   random draw of landmark cells changes? That is the noise floor.
2. How close is Python to Seurat, judged against that floor?
3. If they differ, where? Seurat's own cell coordinates and landmark cells are
   fed through the Python sample-level functions, which separates a difference
   in the coordinates from a difference in the later steps.

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


def load_notebook_functions():
    """Import the notebook file, to use the method functions defined in it."""
    spec = importlib.util.spec_from_file_location(
        "scslide_notebook", NOTEBOOK_DIR / "scslide.py"
    )
    module = importlib.util.module_from_spec(spec)
    sys.modules["scslide_notebook"] = module
    spec.loader.exec_module(module)
    return module


def guide_correlation(
    methods, views, landmark_rows, guide, n_guides, count_landmark_cells=True
):
    """Run the sample-level steps and return the guide-by-guide correlation matrix.

    Parameters
    ----------
    methods
        The notebook module.
    views : list of arrays of shape (n_cells, n_dims)
        The two cell embeddings.
    landmark_rows : int array
        Rows of the landmark cells.
    guide : int array of shape (n_cells,)
        Guide of every cell.
    n_guides : int
    count_landmark_cells : bool
        When False, landmark cells do not vote. The reference's count table
        sums to five times the number of non-landmark cells, so this is what
        Seurat does; the paper's Methods count every cell.
    """
    pool, affinities, repeated = methods.landmark_affinities(
        views, [view[landmark_rows] for view in views]
    )
    landmarks, _ = methods.top_landmarks(
        pool, affinities, repeated, weights=(0.5, 0.5), k=5
    )
    voters = np.ones(len(guide), dtype=bool)
    if not count_landmark_cells:
        voters[landmark_rows] = False
    counts = methods.landmark_counts(
        landmarks[voters], guide[voters], len(landmark_rows), n_guides
    )
    return methods.column_correlation(methods.chi_square_residuals(counts))


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

    # Seurat's guides, reordered to the Python guide order.
    reference_guides = pd.read_csv(args.reference_dir / "guides.csv")
    order = (
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
    seurat = {
        seed: pd.read_csv(
            args.reference_dir / f"correlation_seed{seed}.csv", index_col=0
        )
        .loc[order, order]
        .to_numpy()
        for seed in seeds
    }
    summary = dict(
        line.split(" ", 1)
        for line in (args.reference_dir / f"summary_seed{seeds[0]}.txt")
        .read_text()
        .splitlines()
    )

    # The notebook's default: landmark cells do not vote, as in Seurat's count table.
    python = {}
    for seed in (0, 1, 2):
        rows = methods.sample_landmarks(leverage, n_landmarks, seed=seed)
        python[seed] = guide_correlation(
            methods, views, rows, guide, len(guides), count_landmark_cells=False
        )

    # Does it matter that landmarks are drawn by leverage score rather than at random?
    uniform_rows = np.sort(
        np.random.default_rng(0).choice(len(guide), n_landmarks, replace=False)
    )
    uniform = guide_correlation(
        methods, views, uniform_rows, guide, len(guides), count_landmark_cells=False
    )
    uniform_vs_leverage = upper_triangle_r(uniform, python[0])

    def pairwise(runs):
        keys = sorted(runs)
        return [
            upper_triangle_r(runs[a], runs[b])
            for index, a in enumerate(keys)
            for b in keys[index + 1 :]
        ]

    python_floor = pairwise(python)
    seurat_floor = pairwise(seurat)
    across = [
        upper_triangle_r(run, reference)
        for reference in seurat.values()
        for run in python.values()
    ]
    print("1. Noise floor: the same implementation, a different draw of landmarks")
    print("   Python vs Python:", ", ".join(f"{value:.3f}" for value in python_floor))
    print("   Seurat vs Seurat:", ", ".join(f"{value:.3f}" for value in seurat_floor))
    print("\n2. Python against Seurat (upper triangle of the 138 x 138 matrix)")
    print(
        f"   r = {np.mean(across):.3f} (range {min(across):.3f} to {max(across):.3f} "
        f"over {len(across)} pairs of runs)"
    )

    table = {"Python": describe(methods, python[0], guides)}
    table.update(
        {
            f"Seurat, seed {seed}": describe(methods, run, guides)
            for seed, run in seurat.items()
        }
    )
    print("\n   What each matrix says:")
    print(pd.DataFrame(table).round(3).to_string())

    print(
        f"\n   Landmarks drawn at random instead of by leverage: r = {uniform_vs_leverage:.3f}"
    )

    print("\n3. Where the difference comes from")
    seed = seeds[0]
    kept_barcodes = (
        (ROOT / "data" / "GSE132080" / "kept_barcodes.txt").read_text().split()
    )
    cells = (
        pd.read_csv(args.reference_dir / f"cells_seed{seed}.csv.gz")
        .set_index("cell_barcode")
        .loc[kept_barcodes]
    )
    reference_guide = (
        reference_guides.set_index("id")
        .loc[cells["guide"], "guide_identity"]
        .map(guides.set_index("guide_identity")["guide"])
        .to_numpy()
    )
    if not np.array_equal(reference_guide, guide):
        raise ValueError("the reference cells do not carry the guides of cells.npz")
    reference_views = [
        cells[[name for name in cells.columns if name.startswith("PC_")]].to_numpy(
            np.float32
        ),
        cells[[name for name in cells.columns if name.startswith("pPLS_")]].to_numpy(
            np.float32
        ),
    ]
    reference_landmarks = np.flatnonzero(cells["landmark"].to_numpy())
    agreement = {
        "pca": subspace_agreement(views[0], reference_views[0]),
        "pls": subspace_agreement(views[1], reference_views[1]),
    }
    print(
        "   Same directions? mean canonical correlation, Python vs Seurat: "
        f"PCA {agreement['pca']:.3f}, PLS {agreement['pls']:.3f}"
    )
    experiments = {
        "seurat_coordinates_landmarks_do_not_vote": (reference_views, False),
        "seurat_coordinates_every_cell_votes": (reference_views, True),
        "python_coordinates_landmarks_do_not_vote": (views, False),
        "python_coordinates_every_cell_votes": (views, True),
    }
    transplants = {}
    for name, (chosen_views, count_landmark_cells) in experiments.items():
        result = guide_correlation(
            methods,
            chosen_views,
            reference_landmarks,
            guide,
            len(guides),
            count_landmark_cells,
        )
        transplants[name] = upper_triangle_r(result, seurat[seed])
        print(
            f"   Seurat's landmark cells, {name.replace('_', ' ')}: r = {transplants[name]:.3f}"
        )
        if name == "python_coordinates_landmarks_do_not_vote":
            # The Python pipeline, changed in one thing only: which cells are landmarks.
            table["Python, with Seurat's landmark cells"] = describe(
                methods, result, guides
            )
    print("\n   The same summary, with Seurat's landmark cells in the Python pipeline:")
    print(pd.DataFrame(table).round(3).to_string())

    # The one block that differs, RPS14/15/18, sits on a knife edge between the
    # two views. Sweep the view weight, with Seurat's landmark cells, to show it.
    pool, affinities, repeated = methods.landmark_affinities(
        views, [view[reference_landmarks] for view in views]
    )
    sweep = []
    print("\n4. The weight of the label-aware view, with Seurat's landmark cells")
    for weight in (0.40, 0.45, 0.50, 0.55, 0.60):
        result = methods.run_sample_level(
            pool,
            affinities,
            repeated,
            guide,
            len(guides),
            reference_landmarks,
            weights=(1.0 - weight, weight),
        )["correlation"]
        row = {
            "weight": weight,
            "r_with_seurat": round(upper_triangle_r(result, seurat[seed]), 4),
            "ribosomal": round(
                describe(methods, result, guides)["RPS14/15/18 cross r"], 4
            ),
        }
        sweep.append(row)
        print(
            f"   weight {weight:.2f}: r with Seurat {row['r_with_seurat']:.3f}, RPS14/15/18 {row['ribosomal']:.3f}"
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
        "python_vs_seurat": {
            "mean": round(float(np.mean(across)), 4),
            "min": round(min(across), 4),
            "max": round(max(across), 4),
        },
        "uniform_vs_leverage_landmarks": round(uniform_vs_leverage, 4),
        "subspace_agreement": {
            name: round(value, 4) for name, value in agreement.items()
        },
        "with_seurat_landmark_cells": {
            name: round(value, 4) for name, value in transplants.items()
        },
        "weight_sweep_with_seurat_landmark_cells": sweep,
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
