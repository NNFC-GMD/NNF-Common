"""Turn the public GEO files of Jost et al. 2020 into the small files the notebook ships.

The raw data (GEO series GSE132080, a 336 MB count matrix) is far too large
for a browser, so this script does the heavy steps once, natively, and writes
a few megabytes of derived data to ``public/``:

    cells.npz       per-cell coordinates, landmark weights, each cell's
                    surroundings and target-gene expression
    sgrna.csv       one row per guide RNA, with the activity scores of the original study
    manifest.json   provenance, parameters, summary numbers and checksums

Run it in the ``prep`` environment:

    pixi run -e prep python notebooks/scslide/prepare.py

The method functions (canonical PLS, leverage scores, nearest neighbours) are
imported from the notebook itself, so this script and the published page run
the same code.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import importlib.metadata
import importlib.util
import io
import json
import sys
import urllib.request
import warnings
import zipfile
from pathlib import Path

import anndata
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.io
import umap
from scipy import sparse
from sklearn.decomposition import PCA

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent

GEO_URL = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE132nnn/GSE132080/suppl/"
# File name -> SHA-256 of the file as downloaded on 2026-11-11.
GEO_FILES = {
    "GSE132080_10X_barcodes.tsv.gz": "5ac798d460d21a0ca190f3c11d952ea98fa15ec2f4b2a0ec6d1cde2e7202aa65",
    "GSE132080_10X_genes.tsv.gz": "cb459e6f14a8008e68f2ae535ac5411b6651d278355e9b517ff6a7e43587e603",
    "GSE132080_10X_matrix.mtx.gz": "7c908fbb76feaccac285209bf2578970873f8b8ec819e2a88ac6b83e3553239d",
    "GSE132080_cell_identities.csv.gz": "db9a2bab685004537fccd7aa2b5e19c6c1fc782eca9a4eb1d9275ffe2a19b0a0",
    "GSE132080_sgRNA_barcode_sequences_and_phenotypes.csv.gz": "9ddc95fe2b668fe2bad0acd74c9c65daac6d88b4b414c41a363f240fadb2b505",
}

# The numbers the paper reports for this dataset after quality control.
PAPER_CELLS = 19_587
PAPER_GUIDES = 138

PARAMETERS = {
    "min_genes": 500,
    "min_umis": 1000,
    "scale_factor": 10_000,
    "n_variable_genes": 2000,
    "label_gene_fraction": 0.01,
    "label_gene_min_detection": 0.1,
    "scale_clip": 10.0,
    "n_pca": 30,
    "n_pls": 20,
    "pls_response": "target gene",
    "n_landmarks": 5000,
    "umap_neighbours": 20,
    "neighbour_candidates": 200,
    "bandwidth_rank": 20,
    "seed": 0,
}


def load_notebook_functions():
    """Import the notebook file, to use the method functions defined in it."""
    spec = importlib.util.spec_from_file_location(
        "scslide_notebook", HERE / "scslide.py"
    )
    module = importlib.util.module_from_spec(spec)
    sys.modules["scslide_notebook"] = module
    spec.loader.exec_module(module)
    return module


def sha256_of(path: Path) -> str:
    """Return the SHA-256 of a file, as hexadecimal."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def download(data_dir: Path) -> list[dict]:
    """Fetch the GEO files that are missing and check every file's checksum.

    Returns
    -------
    list of dict
        ``file``, ``sha256`` and ``bytes`` for every source file.

    Raises
    ------
    ValueError
        If a file does not have the recorded checksum.
    """
    data_dir.mkdir(parents=True, exist_ok=True)
    records = []
    for name, expected in GEO_FILES.items():
        path = data_dir / name
        if not path.is_file():
            print(f"downloading {name}", flush=True)
            partial = path.with_suffix(path.suffix + ".part")
            urllib.request.urlretrieve(GEO_URL + name, partial)
            partial.rename(path)
        actual = sha256_of(path)
        if expected and actual != expected:
            raise ValueError(f"{name}: SHA-256 is {actual}, expected {expected}")
        records.append({"file": name, "sha256": actual, "bytes": path.stat().st_size})
    return records


def read_counts(data_dir: Path) -> sparse.csr_matrix:
    """Read the count matrix as cells by genes.

    The parsed matrix is cached next to the download, because reading the
    336 MB text file takes a while.
    """
    cache = data_dir / "counts_cells_by_genes.npz"
    if cache.is_file():
        return sparse.load_npz(cache)
    print("parsing the count matrix", flush=True)
    with gzip.open(data_dir / "GSE132080_10X_matrix.mtx.gz", "rb") as handle:
        genes_by_cells = scipy.io.mmread(handle)
    counts = sparse.csr_matrix(genes_by_cells.T, dtype=np.float32)
    sparse.save_npz(cache, counts)
    return counts


def read_cell_table(data_dir: Path) -> pd.DataFrame:
    """Return one row per barcode of the matrix, with its guide assignment.

    Barcodes without an entry in the identities file keep missing values.
    """
    barcodes = pd.read_csv(
        data_dir / "GSE132080_10X_barcodes.tsv.gz", header=None, names=["cell_barcode"]
    )
    identities = pd.read_csv(data_dir / "GSE132080_cell_identities.csv.gz")
    return barcodes.merge(
        identities, on="cell_barcode", how="left", validate="one_to_one"
    )


def split_guide_identity(identity: str) -> tuple[str, str, str]:
    """Split a guide identity into target gene, short label and sgRNA name.

    ``RPS15_RPS15_+_1438413.23-P1P2_01`` gives ``("RPS15", "RPS15_01",
    "RPS15_+_1438413.23-P1P2_01")``. A non-targeting control such as
    ``neg_ctrl_non-targeting_00028`` gives ``("non-targeting", "NT_00028",
    "non-targeting_00028")``.
    """
    suffix = identity.rsplit("_", 1)[1]
    if identity.startswith("neg_ctrl_"):
        return "non-targeting", f"NT_{suffix}", identity.removeprefix("neg_ctrl_")
    gene, name = identity.split("_", 1)
    return gene, f"{gene}_{suffix}", name


def build_guide_table(cells: pd.DataFrame, data_dir: Path) -> pd.DataFrame:
    """Return one row per guide, ordered by gene then label, controls last."""
    phenotypes = pd.read_csv(
        data_dir / "GSE132080_sgRNA_barcode_sequences_and_phenotypes.csv.gz"
    )
    guides = cells["guide_identity"].value_counts().rename("n_cells").reset_index()
    parts = [split_guide_identity(identity) for identity in guides["guide_identity"]]
    guides[["gene", "label", "sgRNA_name"]] = pd.DataFrame(parts, index=guides.index)
    guides["control"] = guides["gene"] == "non-targeting"
    guides = guides.merge(
        phenotypes.drop(columns="gene"),
        on="sgRNA_name",
        how="left",
        validate="one_to_one",
    )
    missing = guides.loc[~guides["control"], "relative_activity_day10"].isna().sum()
    if missing:
        raise ValueError(f"{missing} targeting guides have no activity score")
    guides = guides.sort_values(["control", "gene", "label"]).reset_index(drop=True)
    guides.insert(0, "guide", np.arange(len(guides)))
    return guides


def label_correlated_genes(
    normalised: sparse.csr_matrix,
    labels: np.ndarray,
    fraction: float,
    min_detection: float,
) -> np.ndarray:
    """Pick the genes most correlated with each label, one label against the rest.

    Parameters
    ----------
    normalised
        Log-normalised expression, cells by genes.
    labels
        One label per cell.
    fraction
        Share of the eligible genes kept for each label, ranked by the
        absolute Pearson correlation with that label's 0/1 indicator.
    min_detection
        A gene is eligible when at least this share of cells express it.

    Returns
    -------
    numpy.ndarray
        Sorted column indices of the union of the picked genes.
    """
    n_cells = normalised.shape[0]
    categories, codes = np.unique(labels, return_inverse=True)
    indicator = np.zeros((n_cells, len(categories)))
    indicator[np.arange(n_cells), codes] = 1.0

    gene_mean = np.asarray(normalised.mean(axis=0)).ravel()
    gene_sq_mean = np.asarray(normalised.multiply(normalised).mean(axis=0)).ravel()
    gene_sd = np.sqrt(np.maximum(gene_sq_mean - gene_mean**2, 0.0))
    label_mean = indicator.mean(axis=0)
    label_sd = indicator.std(axis=0)
    covariance = np.asarray(normalised.T @ indicator) / n_cells - np.outer(
        gene_mean, label_mean
    )
    with np.errstate(divide="ignore", invalid="ignore"):
        correlation = covariance / np.outer(gene_sd, label_sd)

    detection = np.asarray((normalised > 0).mean(axis=0)).ravel()
    eligible = np.flatnonzero((detection >= min_detection) & (gene_sd > 0))
    n_top = int(np.ceil(fraction * len(eligible)))
    strength = np.abs(correlation[eligible])
    top = np.argsort(-strength, axis=0, kind="stable")[:n_top]
    return np.unique(eligible[top.ravel()])


def standardise(
    normalised: sparse.csr_matrix, columns: np.ndarray, clip: float
) -> np.ndarray:
    """Return the chosen genes with mean zero and standard deviation one, capped at ``clip``."""
    dense = normalised[:, columns].toarray().astype(np.float64)
    dense -= dense.mean(axis=0)
    spread = dense.std(axis=0, ddof=1)
    dense /= np.where(spread == 0, 1.0, spread)
    np.minimum(dense, clip, out=dense)
    return dense


def joint_cell_neighbours(
    methods, views: list[np.ndarray], k: int, candidates: int, bandwidth_rank: int
) -> tuple[np.ndarray, np.ndarray]:
    """Find every cell's ``k`` nearest other cells under the equal-weight joint affinity.

    This is the cell-to-cell version of what the notebook does between cells
    and landmarks: the same kernel in each view, averaged with equal weights.
    It is only used to draw the cell-level UMAP.

    Returns
    -------
    indices : int array of shape (n_cells, k + 1)
        The cell itself first, then its neighbours, most similar first.
    distances : float array of shape (n_cells, k + 1)
        Zero for the cell itself, then one minus the joint affinity.
    """
    normalised = [methods.l2_normalize(view) for view in views]
    n_cells = len(normalised[0])
    nearest = [methods.nearest_rows(view, view, candidates + 1) for view in normalised]
    # Column 0 of each search is the cell itself, at distance zero.
    pool = np.concatenate([indices[:, 1:] for indices, _ in nearest], axis=1)
    pool.sort(axis=1)
    repeated = np.zeros(pool.shape, dtype=bool)
    repeated[:, 1:] = pool[:, 1:] == pool[:, :-1]
    joint = np.zeros(pool.shape, dtype=np.float32)
    for view, (_, distances) in zip(normalised, nearest, strict=True):
        first = distances[:, 1]
        width = np.maximum(distances[:, bandwidth_rank] - first, 1e-12)
        for start in range(0, n_cells, 512):
            stop = start + 512
            offsets = view[pool[start:stop]] - view[start:stop, None, :]
            distance = np.sqrt(np.einsum("ijk,ijk->ij", offsets, offsets))
            beyond_first = np.maximum(distance - first[start:stop, None], 0.0)
            joint[start:stop] += 0.5 * np.exp(-beyond_first / width[start:stop, None])
    joint[repeated] = -np.inf
    best = np.argsort(-joint, axis=1, kind="stable")[:, :k]
    indices = np.column_stack(
        [np.arange(n_cells), np.take_along_axis(pool, best, axis=1)]
    )
    distances = np.column_stack(
        [np.zeros(n_cells), 1.0 - np.take_along_axis(joint, best, axis=1)]
    )
    return indices.astype(np.int32), distances.astype(np.float32)


def write_npz(path: Path, arrays: dict[str, np.ndarray]) -> None:
    """Write a compressed ``.npz`` whose bytes depend only on the arrays.

    ``numpy.savez_compressed`` stamps each entry with the current time, so
    two runs would give two different files. Here the time stamp is fixed.
    """
    with zipfile.ZipFile(path, "w", zipfile.ZIP_DEFLATED, compresslevel=9) as archive:
        for name, array in arrays.items():
            buffer = io.BytesIO()
            np.lib.format.write_array(buffer, np.asanyarray(array), allow_pickle=False)
            entry = zipfile.ZipInfo(f"{name}.npy", date_time=(1980, 1, 1, 0, 0, 0))
            entry.compress_type = zipfile.ZIP_DEFLATED
            entry.external_attr = 0o644 << 16
            archive.writestr(entry, buffer.getvalue(), compresslevel=9)


def main(argv: list[str] | None = None) -> int:
    """Run the preparation. Returns the process exit code."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=ROOT / "data" / "GSE132080",
        help="folder for the raw GEO files (default: data/GSE132080)",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=HERE / "public",
        help="folder for the derived files (default: public/ next to this script)",
    )
    args = parser.parse_args(argv)
    methods = load_notebook_functions()
    seed = PARAMETERS["seed"]

    sources = download(args.data_dir)
    counts = read_counts(args.data_dir)
    cells = read_cell_table(args.data_dir)
    if counts.shape[0] != len(cells):
        raise ValueError(f"{counts.shape[0]} matrix rows but {len(cells)} barcodes")
    genes = pd.read_csv(
        args.data_dir / "GSE132080_10X_genes.tsv.gz",
        sep="\t",
        header=None,
        names=["id", "symbol"],
    )

    # Quality control, in the two moves the paper describes.
    cells["n_umi"] = np.asarray(counts.sum(axis=1)).ravel()
    cells["n_genes"] = np.asarray((counts > 0).sum(axis=1)).ravel()
    assigned = (
        (cells["number_of_cells"] == 1)
        & cells["good_coverage"].astype("boolean").fillna(False)
        & cells["guide_identity"].notna()
        & (cells["guide_identity"] != "*")
    )
    n_assigned = int(assigned.sum())
    n_assigned_guides = cells.loc[assigned, "guide_identity"].nunique()
    if (n_assigned, n_assigned_guides) != (PAPER_CELLS, PAPER_GUIDES):
        raise ValueError(
            f"{n_assigned} cells and {n_assigned_guides} guides pass the assignment "
            f"filter; the paper reports {PAPER_CELLS} and {PAPER_GUIDES}"
        )
    deep_enough = (cells["n_genes"] >= PARAMETERS["min_genes"]) & (
        cells["n_umi"] >= PARAMETERS["min_umis"]
    )
    kept = (assigned & deep_enough).to_numpy()
    removed_by_thresholds = int((assigned & ~deep_enough).sum())
    print(
        f"{n_assigned} cells with one assigned guide; the gene and UMI thresholds "
        f"remove {removed_by_thresholds} of them",
        flush=True,
    )
    cells = cells.loc[kept].reset_index(drop=True)
    counts = counts[kept]
    # The barcodes are not shipped; validation/scslide/compare.py reads them
    # here to line its cells up with these.
    (args.data_dir / "kept_barcodes.txt").write_text(
        "\n".join(cells["cell_barcode"]) + "\n", encoding="utf-8"
    )

    guides = build_guide_table(cells, args.data_dir)
    guide_code = cells["guide_identity"].map(
        guides.set_index("guide_identity")["guide"]
    )
    cell_gene = guide_code.map(guides["gene"]).to_numpy()

    # Normalise, and choose the genes each embedding is built from.
    totals = np.asarray(counts.sum(axis=1)).ravel()
    normalised = sparse.diags(PARAMETERS["scale_factor"] / totals) @ counts
    normalised = normalised.tocsr()
    normalised.data = np.log1p(normalised.data)

    adata = anndata.AnnData(counts.copy())
    sc.pp.highly_variable_genes(
        adata, flavor="seurat_v3", n_top_genes=PARAMETERS["n_variable_genes"]
    )
    variable = np.flatnonzero(adata.var["highly_variable"].to_numpy())
    correlated = label_correlated_genes(
        normalised,
        cell_gene,
        PARAMETERS["label_gene_fraction"],
        PARAMETERS["label_gene_min_detection"],
    )
    pls_features = np.union1d(variable, correlated)
    print(
        f"{len(variable)} variable genes, {len(correlated)} label-correlated genes, "
        f"{len(pls_features)} in their union",
        flush=True,
    )

    # The "type" view: PCA of the variable genes.
    scaled_variable = standardise(normalised, variable, PARAMETERS["scale_clip"])
    pca = PCA(n_components=PARAMETERS["n_pca"], svd_solver="full").fit_transform(
        scaled_variable
    )
    del scaled_variable

    # The "state" view: canonical PLS against the target gene.
    genes_of_guides = np.array(sorted(set(cell_gene)))
    responses = (cell_gene[:, None] == genes_of_guides[None, :]).astype(np.float64)
    scaled_pls = standardise(normalised, pls_features, PARAMETERS["scale_clip"])
    print("fitting canonical PLS", flush=True)
    pls = methods.canonical_pls(scaled_pls, responses, PARAMETERS["n_pls"])["scores"]
    del scaled_pls

    # Landmark weights: leverage scores of the variable-gene matrix.
    variable_matrix = normalised[:, variable]
    leverage = methods.sketched_leverage_scores(variable_matrix, seed=seed)
    exact = methods.exact_leverage_scores(variable_matrix.toarray())
    leverage_agreement = float(np.corrcoef(leverage, exact)[0, 1])
    # Stored in single precision, and the landmarks are drawn from the stored
    # values, so the notebook draws exactly the same cells from the file.
    leverage = leverage.astype(np.float32)
    landmark = np.zeros(len(cells), dtype=bool)
    landmark[
        methods.sample_landmarks(leverage, PARAMETERS["n_landmarks"], seed=seed)
    ] = True
    print(
        f"sketched and exact leverage scores correlate {leverage_agreement:.3f}",
        flush=True,
    )

    # Every cell's surroundings among all cells, which Seurat's code uses to
    # judge how near a landmark is. Measured on the coordinates as they are
    # stored, so the notebook gets the same numbers from the file.
    print("measuring every cell's surroundings", flush=True)
    stored_views = [view.astype(np.float16).astype(np.float32) for view in (pca, pls)]
    scale_offset, scale_width = methods.local_scales(
        stored_views, rank=PARAMETERS["bandwidth_rank"]
    )

    # A two-dimensional picture of the cells, from the joint view.
    print("computing the UMAP", flush=True)
    neighbours = joint_cell_neighbours(
        methods,
        [pca, pls],
        k=PARAMETERS["umap_neighbours"] - 1,
        candidates=PARAMETERS["neighbour_candidates"],
        bandwidth_rank=PARAMETERS["bandwidth_rank"],
    )
    joint_view = np.hstack([methods.l2_normalize(pca), methods.l2_normalize(pls)])
    reducer = umap.UMAP(
        n_neighbors=PARAMETERS["umap_neighbours"],
        precomputed_knn=(*neighbours, None),
        random_state=seed,
        n_jobs=1,
        # Keep the layout pipeline of umap-learn 0.5, whatever the default becomes.
        compatibility_layout=True,
    )
    with warnings.catch_warnings():
        # No search index is passed, which only matters for embedding new cells.
        warnings.filterwarnings("ignore", message=".*knn_search_index.*")
        embedding = reducer.fit_transform(joint_view)

    # Expression of every target gene in every cell, for the knock-down check.
    target_genes = [gene for gene in genes_of_guides if gene != "non-targeting"]
    symbol_to_column = pd.Series(np.arange(len(genes)), index=genes["symbol"])
    missing = [gene for gene in target_genes if gene not in symbol_to_column.index]
    if missing:
        raise ValueError(f"target genes absent from the matrix: {missing}")
    target_columns = [
        int(np.atleast_1d(symbol_to_column[gene])[0]) for gene in target_genes
    ]
    target_expression = normalised[:, target_columns].toarray()

    args.out.mkdir(parents=True, exist_ok=True)
    write_npz(
        args.out / "cells.npz",
        {
            "pca": pca.astype(np.float16),
            "pls": pls.astype(np.float16),
            "umap": embedding.astype(np.float16),
            "leverage": leverage,
            "landmark": landmark,
            "scale_offset": scale_offset.astype(np.float32),
            "scale_width": scale_width.astype(np.float32),
            "guide": guide_code.to_numpy().astype(np.int16),
            "target_expression": target_expression.astype(np.float16),
            "target_genes": np.array(target_genes),
            "n_umi": cells["n_umi"].to_numpy().astype(np.int32),
            "n_genes": cells["n_genes"].to_numpy().astype(np.int16),
            "gemgroup": cells["gemgroup"].to_numpy().astype(np.int8),
        },
    )
    guides[
        [
            "guide",
            "label",
            "gene",
            "control",
            "n_cells",
            "sequence",
            "gamma_day5",
            "gamma_day10",
            "relative_activity_day5",
            "relative_activity_day10",
            "guide_identity",
        ]
    ].to_csv(args.out / "sgrna.csv", index=False, lineterminator="\n")

    manifest = {
        "source": {
            "accession": "GSE132080",
            "citation": "Jost et al., Nature Biotechnology 38, 355-364 (2020)",
            "doi": "10.1038/s41587-019-0387-5",
            "url": GEO_URL,
            "files": sources,
        },
        "parameters": PARAMETERS,
        "summary": {
            "cells_with_one_assigned_guide": n_assigned,
            "cells_removed_by_gene_and_umi_thresholds": removed_by_thresholds,
            "cells": len(cells),
            "guides": len(guides),
            "target_genes": len(target_genes),
            "control_guides": int(guides["control"].sum()),
            "variable_genes": len(variable),
            "label_correlated_genes": len(correlated),
            "pls_features": len(pls_features),
            "landmarks": int(landmark.sum()),
            "leverage_sketch_vs_exact_correlation": round(leverage_agreement, 4),
        },
        "versions": {
            name: importlib.metadata.version(name)
            for name in (
                "numpy",
                "scipy",
                "pandas",
                "scikit-learn",
                "scanpy",
                "anndata",
                "umap-learn",
            )
        },
        "files": [
            {
                "file": name,
                "sha256": sha256_of(args.out / name),
                "bytes": (args.out / name).stat().st_size,
            }
            for name in ("cells.npz", "sgrna.csv")
        ],
    }
    (args.out / "manifest.json").write_text(
        json.dumps(manifest, indent=4) + "\n", encoding="utf-8"
    )
    for entry in manifest["files"]:
        print(f"wrote {entry['file']}: {entry['bytes'] / 1e6:.2f} MB")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
