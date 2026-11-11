# /// script
# requires-python = ">=3.12"
# dependencies = [
#     "altair",
#     "marimo",
#     "numpy",
#     "pandas",
#     "scikit-learn",
#     "scipy",
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
# scipy = "1.17.*"
# pandas = "3.0.*"
# altair = "6.0.*"
# scikit-learn = "1.8.*"
# ///
"""scSLIDE explained and replicated in Python.

The method is from Jiang et al., "Reconstructing developmental and disease
progression with sample-level embeddings" (Cell, doi 10.1016/j.cell.2026.10.018).
Every function below is written from that paper's STAR Methods. Nothing is
copied or translated from the Seurat source code.
"""

import marimo

__generated_with = "0.25.1"
app = marimo.App(width="medium", app_title="scSLIDE, step by step")

with app.setup:
    import io
    import json
    import platform
    import sys
    import urllib.request

    import altair as alt
    import marimo as mo
    import numpy as np
    import pandas as pd
    from scipy import sparse
    from scipy.cluster import hierarchy
    from scipy.stats import rankdata
    from sklearn.cluster import KMeans
    from sklearn.manifold import TSNE

    # The charts here hold more rows than Altair's default limit of 5000.
    alt.data_transformers.disable_max_rows()

    # Okabe-Ito, a colour-blind-safe palette for categories.
    OKABE_ITO = {
        "orange": "#E69F00",
        "sky": "#56B4E9",
        "green": "#009E73",
        "yellow": "#F0E442",
        "blue": "#0072B2",
        "vermillion": "#D55E00",
        "purple": "#CC79A7",
        "black": "#000000",
    }
    GREY = "#B8B8B8"


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
def l2_normalize(matrix):
    """Scale every row to length one. A row of zeros stays a row of zeros."""
    matrix = np.asarray(matrix, dtype=np.float32)
    norms = np.linalg.norm(matrix, axis=1, keepdims=True)
    return matrix / np.where(norms == 0, 1.0, norms)


@app.function
def nearest_rows(queries, references, k, chunk_size=1024):
    """Find the ``k`` nearest reference rows of every query row.

    Parameters
    ----------
    queries : array of shape (n_queries, n_dims)
    references : array of shape (n_references, n_dims)
    k : int
        Number of neighbours to keep; capped at ``n_references``.
    chunk_size : int
        Query rows handled at a time, which bounds the memory used.

    Returns
    -------
    indices : int array of shape (n_queries, k)
        Reference rows, nearest first.
    distances : float array of shape (n_queries, k)
        The matching Euclidean distances.
    """
    queries = np.asarray(queries, dtype=np.float32)
    references = np.asarray(references, dtype=np.float32)
    k = min(k, len(references))
    reference_sq = np.einsum("ij,ij->i", references, references)
    indices = np.empty((len(queries), k), dtype=np.int32)
    distances = np.empty((len(queries), k), dtype=np.float32)
    for start in range(0, len(queries), chunk_size):
        block = queries[start : start + chunk_size]
        squared = (
            np.einsum("ij,ij->i", block, block)[:, None]
            + reference_sq[None, :]
            - 2.0 * (block @ references.T)
        )
        np.maximum(squared, 0.0, out=squared)
        if k < squared.shape[1]:
            nearest = np.argpartition(squared, k - 1, axis=1)[:, :k]
        else:
            nearest = np.tile(np.arange(k), (len(block), 1))
        nearest_sq = np.take_along_axis(squared, nearest, axis=1)
        order = np.argsort(nearest_sq, axis=1, kind="stable")
        indices[start : start + chunk_size] = np.take_along_axis(nearest, order, axis=1)
        distances[start : start + chunk_size] = np.sqrt(
            np.take_along_axis(nearest_sq, order, axis=1)
        )
    return indices, distances


@app.function
def landmark_affinities(
    cell_views, landmark_views, candidates=200, bandwidth_rank=20, chunk_size=512
):
    """Score how close every cell is to its candidate landmarks, in each view.

    Each view is one embedding of the same cells, for example PCA and PLS.
    Rows are L2-normalised first. In every view a cell's candidates are its
    ``candidates`` nearest landmarks; the candidate pool of a cell is the
    union over views. The affinity of a cell to a landmark in one view is

        exp(-max(d - d_first, 0) / (d_bandwidth - d_first))

    where ``d`` is their Euclidean distance, ``d_first`` the cell's distance
    to its nearest landmark and ``d_bandwidth`` its distance to its
    ``bandwidth_rank``-th nearest landmark, all in that view.

    Parameters
    ----------
    cell_views : sequence of arrays of shape (n_cells, n_dims_view)
    landmark_views : sequence of arrays of shape (n_landmarks, n_dims_view)
        The landmark rows of the same embeddings, in the same order.
    candidates : int
        Nearest landmarks kept per cell and view.
    bandwidth_rank : int
        Rank of the neighbour whose distance sets the width of the kernel.
    chunk_size : int
        Cells handled at a time, which bounds the memory used.

    Returns
    -------
    pool : int array of shape (n_cells, n_views * candidates)
        Candidate landmarks of every cell, sorted by index. A landmark that
        is a candidate in several views appears several times.
    affinities : float array of shape (n_views, n_cells, n_views * candidates)
        Affinity of each cell to each pooled landmark, per view.
    repeated : bool array, same shape as ``pool``
        True where an entry repeats the landmark to its left.
    """
    cells = [l2_normalize(view) for view in cell_views]
    landmarks = [l2_normalize(view) for view in landmark_views]
    n_cells = len(cells[0])
    candidates = min(candidates, len(landmarks[0]))
    bandwidth_rank = min(bandwidth_rank, candidates)

    nearest = [
        nearest_rows(cell, landmark, candidates, chunk_size=2 * chunk_size)
        for cell, landmark in zip(cells, landmarks, strict=True)
    ]
    pool = np.concatenate([indices for indices, _ in nearest], axis=1)
    pool.sort(axis=1)
    repeated = np.zeros(pool.shape, dtype=bool)
    repeated[:, 1:] = pool[:, 1:] == pool[:, :-1]

    affinities = np.empty((len(cells), *pool.shape), dtype=np.float32)
    for view, (cell, landmark, (_, distances)) in enumerate(
        zip(cells, landmarks, nearest, strict=True)
    ):
        first = distances[:, 0]
        width = np.maximum(distances[:, bandwidth_rank - 1] - first, 1e-12)
        for start in range(0, n_cells, chunk_size):
            stop = start + chunk_size
            offsets = landmark[pool[start:stop]] - cell[start:stop, None, :]
            distance = np.sqrt(np.einsum("ijk,ijk->ij", offsets, offsets))
            beyond_first = np.maximum(distance - first[start:stop, None], 0.0)
            affinities[view, start:stop] = np.exp(
                -beyond_first / width[start:stop, None]
            )
    return pool, affinities, repeated


@app.function
def top_landmarks(pool, affinities, repeated, weights=(0.5, 0.5), k=5):
    """Pick every cell's ``k`` landmarks with the highest weighted affinity.

    Parameters
    ----------
    pool, affinities, repeated
        The output of :func:`landmark_affinities`.
    weights : sequence of float
        One weight per view. scSLIDE fixes them at 0.5 and 0.5.
    k : int
        Landmarks kept per cell.

    Returns
    -------
    landmarks : int array of shape (n_cells, k)
        Landmark indices, highest joint affinity first.
    joint : float array of shape (n_cells, k)
        The matching joint affinities.
    """
    weights = np.asarray(weights, dtype=np.float32)
    joint = np.tensordot(weights, affinities, axes=1)
    joint[repeated] = -np.inf
    best = np.argpartition(-joint, k - 1, axis=1)[:, :k]
    best_joint = np.take_along_axis(joint, best, axis=1)
    order = np.argsort(-best_joint, axis=1, kind="stable")
    landmarks = np.take_along_axis(
        np.take_along_axis(pool, best, axis=1), order, axis=1
    )
    return landmarks, np.take_along_axis(best_joint, order, axis=1)


@app.function
def landmark_counts(landmarks, sample_codes, n_landmarks, n_samples):
    """Count, for every sample, how many of its cells sit next to each landmark.

    Parameters
    ----------
    landmarks : int array of shape (n_cells, k)
        Each cell's landmarks, from :func:`top_landmarks`.
    sample_codes : int array of shape (n_cells,)
        The sample of every cell, as an integer from 0 to ``n_samples - 1``.
    n_landmarks, n_samples : int

    Returns
    -------
    int array of shape (n_landmarks, n_samples)
        Every cell adds one to each of its ``k`` landmarks, in its sample's column.
    """
    counts = np.zeros((n_landmarks, n_samples), dtype=np.int64)
    columns = np.repeat(np.asarray(sample_codes), landmarks.shape[1])
    np.add.at(counts, (landmarks.ravel(), columns), 1)
    return counts


@app.function
def chi_square_residuals(counts):
    """Turn a landmark-by-sample count table into chi-square residuals.

    The expected count of a cell of the table is ``row total * column total /
    grand total``: what we would see if every sample spread its cells over
    the landmarks in the same proportions. The residual is ``(observed -
    expected) / sqrt(expected)``. It is positive where a sample has more
    cells near a landmark than expected, negative where it has fewer.

    Parameters
    ----------
    counts : array of shape (n_landmarks, n_samples)

    Returns
    -------
    float array of the same shape
        Zero wherever the expected count is zero.
    """
    counts = np.asarray(counts, dtype=np.float64)
    total = counts.sum()
    expected = np.outer(counts.sum(axis=1), counts.sum(axis=0)) / total
    residuals = np.zeros_like(counts)
    filled = expected > 0
    residuals[filled] = (counts[filled] - expected[filled]) / np.sqrt(expected[filled])
    return residuals


@app.function
def column_correlation(matrix):
    """Pearson correlation between the columns of ``matrix``.

    A constant column correlates zero with everything, itself included.
    """
    matrix = np.asarray(matrix, dtype=np.float64)
    centred = matrix - matrix.mean(axis=0, keepdims=True)
    norms = np.linalg.norm(centred, axis=0)
    unit = centred / np.where(norms == 0, 1.0, norms)
    return unit.T @ unit


@app.function
def cosine_distance(matrix):
    """One minus the cosine similarity between the columns of ``matrix``."""
    matrix = np.asarray(matrix, dtype=np.float64)
    norms = np.linalg.norm(matrix, axis=0)
    unit = matrix / np.where(norms == 0, 1.0, norms)
    return np.clip(1.0 - unit.T @ unit, 0.0, 2.0)


@app.function
def diffusion_map(distances, n_components=3, neighbours=7):
    """Embed samples with a diffusion map.

    A random walk is set up on the samples: from each sample it steps to a
    nearby one, with a probability that falls off with their distance. The
    leading eigenvectors of that walk are coordinates along which the walk
    mixes slowest. They follow the chain of neighbouring samples, even when
    that chain bends.

    The kernel is ``exp(-d_ij^2 / (s_i * s_j))`` with ``s_i`` the distance of
    sample ``i`` to its ``neighbours``-th nearest sample. It is kept only
    between a sample and its ``neighbours`` nearest samples, so that samples
    far apart are not linked at all. It is then divided by the local density
    on both sides, so unevenly sampled regions do not bend the result, and
    normalised into transition probabilities.

    Parameters
    ----------
    distances : array of shape (n_samples, n_samples)
        Symmetric pairwise distances with a zero diagonal.
    n_components : int
        Diffusion components to return.
    neighbours : int
        Rank of the neighbour that sets each sample's kernel width, and the
        number of neighbours each sample stays linked to. Seven is the usual
        choice for such self-tuning kernels (Zelnik-Manor and Perona, 2004).

    Returns
    -------
    components : float array of shape (n_samples, n_components)
        Each component is its eigenvector scaled by the square root of its
        eigenvalue, and flipped so that its largest entry is positive.
    eigenvalues : float array of shape (n_components,)

    Notes
    -----
    The samples must form one connected group. If they fall into separate
    groups with no links between them, the leading components only say which
    group a sample belongs to.
    """
    distances = np.asarray(distances, dtype=np.float64)
    neighbours = int(np.clip(neighbours, 1, len(distances) - 1))
    reach = np.sort(distances, axis=1)[:, neighbours]
    width = np.maximum(reach, 1e-12)
    kernel = np.exp(-(distances**2) / np.outer(width, width))
    near = distances <= reach[:, None]
    kernel = kernel * (near | near.T)
    density = kernel.sum(axis=1)
    kernel = kernel / np.outer(density, density)
    degree = kernel.sum(axis=1)
    symmetric = kernel / np.sqrt(np.outer(degree, degree))
    eigenvalues, eigenvectors = np.linalg.eigh(symmetric)
    order = np.argsort(eigenvalues)[::-1][1 : n_components + 1]
    eigenvalues = eigenvalues[order]
    vectors = eigenvectors[:, order] / np.sqrt(degree)[:, None]
    vectors = vectors / np.linalg.norm(vectors, axis=0, keepdims=True)
    largest = np.abs(vectors).argmax(axis=0)
    vectors = vectors * np.sign(vectors[largest, np.arange(vectors.shape[1])])
    return vectors * np.sqrt(np.clip(eigenvalues, 0.0, None)), eigenvalues


@app.function
def local_linear_smooth(positions, values, at, span=0.4):
    """Smooth ``values`` against ``positions`` with local straight-line fits.

    Each fit uses the ``span`` fraction of points nearest to the place being
    estimated, weighted so that nearer points count more.

    Parameters
    ----------
    positions : array of shape (n_points,)
    values : array of shape (n_points, n_dims)
    at : array of shape (n_out,)
        Positions at which to estimate.
    span : float
        Fraction of the points used by each fit.

    Returns
    -------
    float array of shape (n_out, n_dims)
    """
    positions = np.asarray(positions, dtype=np.float64)
    values = np.asarray(values, dtype=np.float64)
    at = np.asarray(at, dtype=np.float64)
    n_local = int(np.clip(np.ceil(span * len(positions)), 3, len(positions)))
    gaps = np.abs(at[:, None] - positions[None, :])
    # A hair wider than the farthest of the local points, so that points tied
    # at that distance keep a little weight and the fit is always defined.
    reach = np.maximum(np.sort(gaps, axis=1)[:, n_local - 1], 1e-12) * 1.001
    weights = np.clip(1.0 - (gaps / reach[:, None]) ** 3, 0.0, None) ** 3
    total = weights.sum(axis=1)
    mean_position = (weights @ positions) / total
    mean_value = (weights @ values) / total[:, None]
    offset = positions[None, :] - mean_position[:, None]
    spread = np.maximum((weights * offset**2).sum(axis=1), 1e-12)
    slope = ((weights * offset) @ values) / spread[:, None]
    return mean_value + slope * (at - mean_position)[:, None]


@app.function
def project_onto_curve(points, curve, stretch=2.0):
    """Find the closest place on a polyline for every point.

    Parameters
    ----------
    points : array of shape (n_points, n_dims)
    curve : array of shape (n_vertices, n_dims)
        Consecutive vertices of the polyline.
    stretch : float
        How far a point may project beyond either end, in lengths of the end
        segment. Without it, every point lying past an end would be given
        the same position.

    Returns
    -------
    float array of shape (n_points,)
        The distance along the polyline, from its first vertex, of the place
        closest to each point. It is negative for a point that projects
        before the first vertex.
    """
    starts, ends = curve[:-1], curve[1:]
    segment = ends - starts
    length_sq = np.maximum(np.einsum("ij,ij->i", segment, segment), 1e-24)
    along = np.einsum("pij,ij->pi", points[:, None, :] - starts[None, :, :], segment)
    lowest = np.zeros(len(segment))
    highest = np.ones(len(segment))
    lowest[0] -= stretch
    highest[-1] += stretch
    fraction = np.clip(along / length_sq[None, :], lowest[None, :], highest[None, :])
    closest = starts[None, :, :] + fraction[:, :, None] * segment[None, :, :]
    gap_sq = ((points[:, None, :] - closest) ** 2).sum(axis=2)
    best = gap_sq.argmin(axis=1)
    lengths = np.sqrt(length_sq)
    before = np.concatenate([[0.0], np.cumsum(lengths)])[best]
    return before + fraction[np.arange(len(points)), best] * lengths[best]


@app.function
def principal_curve(points, start=None, iterations=10, span=0.4, n_vertices=100):
    """Fit a smooth curve through the middle of a cloud of points.

    This is the principal-curve idea of Hastie and Stuetzle: start from a
    first guess of every point's position, then alternate between smoothing
    each coordinate against the position and re-projecting every point onto
    the new curve.

    Parameters
    ----------
    points : array of shape (n_points, n_dims)
    start : array of shape (n_points,), optional
        First guess of the positions. The default is the projection onto the
        first principal component. For diffusion components, pass the first
        component: it already orders the samples.
    iterations : int
    span : float
        Smoothness, passed to :func:`local_linear_smooth`.
    n_vertices : int
        Vertices of the returned polyline.

    Returns
    -------
    position : float array of shape (n_points,)
        Position of every point along the curve, from 0 to 1.
    curve : float array of shape (n_vertices, n_dims)
    """
    points = np.asarray(points, dtype=np.float64)
    if start is None:
        centred = points - points.mean(axis=0)
        _, _, directions = np.linalg.svd(centred, full_matrices=False)
        position = centred @ directions[0]
    else:
        position = np.asarray(start, dtype=np.float64)
    curve = None
    for _ in range(iterations):
        grid = np.linspace(position.min(), position.max(), n_vertices)
        curve = local_linear_smooth(position, points, grid, span=span)
        position = project_onto_curve(points, curve)
    length = max(position.max() - position.min(), 1e-12)
    return (position - position.min()) / length, curve


@app.function
def canonical_pls(features, responses, n_components):
    """Canonical partial least squares (Indahl, Liland and Naes, 2009).

    PLS looks for directions in the features that co-vary with the responses.
    The canonical variant builds each direction in two moves. It first forms
    one candidate score per response column, ``Z = X X' Y``. It then takes
    the combination of those candidates that has the highest canonical
    correlation with the responses. The found score is removed from the
    features before the next direction is built.

    Parameters
    ----------
    features : array of shape (n_samples, n_features)
        Usually standardised columns. Centred here in any case.
    responses : array of shape (n_samples, n_responses)
        For a categorical label, one 0/1 column per category.
    n_components : int

    Returns
    -------
    dict
        ``scores`` (n_samples, n_components): the new coordinates.
        ``weights`` (n_features, n_components): the unit-length directions.
        ``loadings`` (n_features, n_components): the features regressed on each score.
        ``projection`` (n_features, n_components): multiply centred features
        by it to get the scores of new samples.
        ``feature_means`` (n_features,): the column means removed.
    """
    remaining = np.array(features, dtype=np.float64)
    feature_means = remaining.mean(axis=0)
    remaining -= feature_means
    responses = np.asarray(responses, dtype=np.float64)
    responses = responses - responses.mean(axis=0)
    response_basis, response_values, _ = np.linalg.svd(responses, full_matrices=False)
    response_basis = response_basis[:, response_values > response_values[0] * 1e-10]

    n_samples, n_features = remaining.shape
    scores = np.zeros((n_samples, n_components))
    weights = np.zeros((n_features, n_components))
    loadings = np.zeros((n_features, n_components))
    for component in range(n_components):
        covariance = remaining.T @ responses
        candidates = remaining @ covariance
        basis, values, right = np.linalg.svd(candidates, full_matrices=False)
        kept = values > values[0] * 1e-10
        leading, _, _ = np.linalg.svd(
            basis[:, kept].T @ response_basis, full_matrices=False
        )
        weight = covariance @ ((right[kept].T / values[kept]) @ leading[:, 0])
        weight /= np.linalg.norm(weight)
        score = remaining @ weight
        loading = remaining.T @ score / (score @ score)
        for start in range(0, n_features, 256):
            remaining[:, start : start + 256] -= np.outer(
                score, loading[start : start + 256]
            )
        scores[:, component] = score
        weights[:, component] = weight
        loadings[:, component] = loading
    projection = weights @ np.linalg.inv(loadings.T @ weights)
    return {
        "scores": scores,
        "weights": weights,
        "loadings": loadings,
        "projection": projection,
        "feature_means": feature_means,
    }


@app.function
def exact_leverage_scores(matrix):
    """Leverage score of every row: how much it alone shapes the column space.

    It is the squared length of the row in an orthonormal basis of the
    columns. The scores lie between 0 and 1 and sum to the rank of the matrix.
    """
    basis, values, _ = np.linalg.svd(
        np.asarray(matrix, dtype=np.float64), full_matrices=False
    )
    basis = basis[:, values > values[0] * 1e-10]
    return np.einsum("ij,ij->i", basis, basis)


@app.function
def sketched_leverage_scores(matrix, sketch_rows=None, projection_size=64, seed=0):
    """Approximate the leverage score of every row without a full decomposition.

    Three moves, as in the paper's Methods. A CountSketch folds the rows into
    far fewer rows by adding each one, with a random sign, to a random bucket.
    A QR decomposition of that small matrix gives a triangular factor ``R``
    for which the rows of ``A R^-1`` are close to orthonormal coordinates.
    A random projection to ``projection_size`` columns keeps their lengths
    nearly intact, and the squared lengths are the scores.

    Parameters
    ----------
    matrix : dense array or scipy sparse matrix of shape (n_rows, n_cols)
    sketch_rows : int, optional
        Buckets of the CountSketch. The default is four times the number of
        columns, and never more than the number of rows.
    projection_size : int
        Columns of the random projection.
    seed : int

    Returns
    -------
    float array of shape (n_rows,)
    """
    rng = np.random.default_rng(seed)
    n_rows, n_cols = matrix.shape
    if sketch_rows is None:
        sketch_rows = 4 * n_cols
    sketch_rows = int(min(sketch_rows, n_rows))
    buckets = rng.integers(0, sketch_rows, size=n_rows)
    signs = rng.choice(np.array([-1.0, 1.0]), size=n_rows)
    count_sketch = sparse.csr_matrix(
        (signs, (buckets, np.arange(n_rows))), shape=(sketch_rows, n_rows)
    )
    folded = count_sketch @ matrix
    if sparse.issparse(folded):
        folded = folded.toarray()
    triangular = np.linalg.qr(np.asarray(folded, dtype=np.float64), mode="r")
    gaussian = rng.normal(size=(triangular.shape[0], projection_size)) / np.sqrt(
        projection_size
    )
    # Least squares instead of a triangular solve, so a rank-deficient R is handled.
    reduced = np.linalg.lstsq(triangular, gaussian, rcond=None)[0]
    coordinates = np.asarray(matrix @ reduced)
    return np.einsum("ij,ij->i", coordinates, coordinates)


@app.function
def sample_landmarks(scores, n_landmarks, seed=0):
    """Draw landmark rows without replacement, with probability proportional to ``scores``.

    Returns the sorted row indices.
    """
    rng = np.random.default_rng(seed)
    scores = np.clip(np.asarray(scores, dtype=np.float64), 0.0, None)
    n_landmarks = min(n_landmarks, int((scores > 0).sum()))
    chosen = rng.choice(
        len(scores), size=n_landmarks, replace=False, p=scores / scores.sum()
    )
    return np.sort(chosen)


@app.function
def nearest_label_votes(distances, labels, k=3):
    """Predict every sample's label from its ``k`` nearest other samples.

    The prediction is the most common label among the neighbours; a tie goes
    to the label of the nearer neighbour.

    Parameters
    ----------
    distances : array of shape (n_samples, n_samples)
    labels : sequence of length n_samples
    k : int

    Returns
    -------
    object array of shape (n_samples,)
    """
    distances = np.array(distances, dtype=np.float64)
    np.fill_diagonal(distances, np.inf)
    labels = np.asarray(labels, dtype=object)
    predictions = np.empty(len(labels), dtype=object)
    for row, order in enumerate(np.argsort(distances, axis=1, kind="stable")[:, :k]):
        neighbour_labels = labels[order]
        tally = {}
        for rank, label in enumerate(neighbour_labels):
            votes, best_rank = tally.get(label, (0, rank))
            tally[label] = (votes + 1, best_rank)
        predictions[row] = max(
            tally, key=lambda label: (tally[label][0], -tally[label][1])
        )
    return predictions


@app.function
def spearman(first, second):
    """Spearman rank correlation between two sequences."""
    return float(np.corrcoef(rankdata(first), rankdata(second))[0, 1])


@app.function
def run_sample_level(
    pool,
    affinities,
    repeated,
    sample_codes,
    n_samples,
    landmark_rows,
    weights=(0.5, 0.5),
    k=5,
    landmarks_vote=False,
    keep=None,
):
    """Go from cell-to-landmark affinities to the sample-by-sample correlation.

    Parameters
    ----------
    pool, affinities, repeated
        The output of :func:`landmark_affinities`.
    sample_codes : int array of shape (n_cells,)
    n_samples : int
    landmark_rows : int array
        Rows of the cells that serve as landmarks.
    weights : sequence of float
        One weight per view.
    k : int
        Landmarks every cell votes for.
    landmarks_vote : bool
        Whether the landmark cells themselves are counted. A landmark's
        nearest landmark is itself, so it would always vote for itself.
    keep : bool array of shape (n_cells,), optional
        Cells to count; all of them by default.

    Returns
    -------
    dict
        ``landmarks`` (n_cells, k), ``counts`` and ``residuals``
        (n_landmarks, n_samples), and ``correlation`` (n_samples, n_samples).
    """
    landmarks, _ = top_landmarks(pool, affinities, repeated, weights=weights, k=k)
    voters = np.ones(len(sample_codes), dtype=bool) if keep is None else keep.copy()
    if not landmarks_vote:
        voters[landmark_rows] = False
    counts = landmark_counts(
        landmarks[voters], sample_codes[voters], len(landmark_rows), n_samples
    )
    residuals = chi_square_residuals(counts)
    return {
        "landmarks": landmarks,
        "counts": counts,
        "residuals": residuals,
        "correlation": column_correlation(residuals),
    }


@app.function
def score_distances(distances, gene, control, activity, k=3):
    """Two checks on a guide-by-guide distance matrix.

    Parameters
    ----------
    distances : array of shape (n_guides, n_guides)
    gene : array of shape (n_guides,)
        Target gene of every guide.
    control : bool array of shape (n_guides,)
        True for non-targeting guides.
    activity : array of shape (n_guides,)
        Independent activity score of every targeting guide.
    k : int

    Returns
    -------
    dict
        ``same_gene``: share of targeting guides whose ``k`` nearest guides
        vote for the right gene. ``activity``: Spearman correlation, over
        targeting guides, between closeness to the controls and activity.
        A good method gives a strongly negative value: weak guides sit near
        the controls.
    """
    targeting = ~control
    predicted = nearest_label_votes(distances, gene, k=k)
    near_controls = -distances[:, control].mean(axis=1)
    return {
        "same_gene": float((predicted[targeting] == gene[targeting]).mean()),
        "activity": spearman(near_controls[targeting], activity[targeting]),
    }


@app.function
def scaled_cosine_distance(profiles):
    """Cosine distance between rows, after standardising every column across rows."""
    profiles = np.asarray(profiles, dtype=np.float64)
    spread = profiles.std(axis=0)
    scaled = (profiles - profiles.mean(axis=0)) / np.where(spread == 0, 1.0, spread)
    return cosine_distance(scaled.T)


@app.function
def wasserstein_on_a_line(first, second):
    """Earth mover's distance between two histograms over evenly spaced positions.

    On a line the cheapest way to turn one pile into the other is known in
    closed form: the area between the two cumulative curves.
    """
    first = np.asarray(first, dtype=np.float64)
    second = np.asarray(second, dtype=np.float64)
    return float(
        np.abs(np.cumsum(first / first.sum()) - np.cumsum(second / second.sum())).sum()
    )


@app.function
def sinkhorn_distances(histograms, cost, regularisation=0.05, iterations=300):
    """Optimal-transport cost between every pair of histograms.

    The cost of moving one unit of mass from bin ``i`` to bin ``j`` is
    ``cost[i, j]``. The Sinkhorn iteration finds, for every pair of
    histograms, a transport plan that is nearly the cheapest; the small
    blur set by ``regularisation`` is what makes it fast.

    Parameters
    ----------
    histograms : array of shape (n_samples, n_bins)
        Non-negative; every row is normalised to sum to one.
    cost : array of shape (n_bins, n_bins)
    regularisation : float
        Blur, as a fraction of the largest cost.
    iterations : int

    Returns
    -------
    float array of shape (n_samples, n_samples)
    """
    histograms = np.asarray(histograms, dtype=np.float64) + 1e-9
    histograms = histograms / histograms.sum(axis=1, keepdims=True)
    cost = np.asarray(cost, dtype=np.float64)
    kernel = np.exp(-cost / (regularisation * cost.max()))
    source = np.broadcast_to(
        histograms[:, None, :], (len(histograms), *histograms.shape)
    )
    target = np.broadcast_to(histograms[None, :, :], source.shape)
    left = np.ones(source.shape)
    for _ in range(iterations):
        right = target / (left @ kernel)
        left = source / (right @ kernel.T)
    distances = np.einsum("abi,ij,abj->ab", left, kernel * cost, right)
    distances = 0.5 * (distances + distances.T)
    np.fill_diagonal(distances, 0.0)
    return distances


@app.function
def cluster_then_order(correlation, gene, control):
    """Order guides for a heat map: genes by similarity, guides within a gene by similarity.

    Returns
    -------
    order : int array
        Guide rows in display order; the controls come last.
    gene_order : list of str
        The genes in display order.
    """
    genes = sorted(set(gene[~control]))
    between = np.array(
        [
            [
                correlation[np.ix_(gene == first, gene == second)].mean()
                for second in genes
            ]
            for first in genes
        ]
    )
    np.fill_diagonal(between, 1.0)
    condensed = (1.0 - between)[np.triu_indices(len(genes), 1)]
    gene_order = [
        genes[i] for i in hierarchy.leaves_list(hierarchy.linkage(condensed, "average"))
    ]
    gene_order.append(gene[control][0])
    order = []
    for name in gene_order:
        rows = np.flatnonzero(gene == name)
        if len(rows) > 2:
            within = (1.0 - correlation[np.ix_(rows, rows)])[
                np.triu_indices(len(rows), 1)
            ]
            rows = rows[hierarchy.leaves_list(hierarchy.linkage(within, "average"))]
        order.extend(rows.tolist())
    return np.array(order), gene_order


@app.function
def markdown_table(frame):
    """Show a small data frame as a plain table.

    A plain table is part of the page itself. marimo's interactive table
    asks the running notebook for its rows, which fails in the preview a
    reader sees before Python has started.
    """
    header = "| " + " | ".join(str(name) for name in frame.columns) + " |"
    rule = "|" + "|".join(" --- " for _ in frame.columns) + "|"
    rows = [
        "| " + " | ".join(str(value) for value in row) + " |"
        for row in frame.itertuples(index=False)
    ]
    return mo.md("\n".join([header, rule, *rows]))


@app.function
def highlight_layers(frame, x, y, focus_field, focus_title, tooltip, size=28):
    """Scatter plot with grey background points and coloured focus points.

    ``frame`` needs a ``layer`` column holding "background", "control" or
    "focus". Focus points are coloured by ``focus_field``.
    """
    base = alt.Chart(frame).encode(
        x=alt.X(f"{x}:Q", axis=None), y=alt.Y(f"{y}:Q", axis=None), tooltip=tooltip
    )
    background = base.transform_filter(alt.datum.layer == "background").mark_circle(
        size=size * 0.4, color=GREY, opacity=0.5
    )
    control = base.transform_filter(alt.datum.layer == "control").mark_circle(
        size=size * 0.5, color=OKABE_ITO["black"], opacity=0.55
    )
    focus = (
        base.transform_filter(alt.datum.layer == "focus")
        .mark_circle(size=size, opacity=0.9)
        .encode(
            color=alt.Color(
                f"{focus_field}:N",
                title=focus_title,
                scale=alt.Scale(range=list(OKABE_ITO.values())[:7]),
            )
        )
    )
    return background + control + focus


@app.cell(hide_code=True)
def _():
    mo.md("""
    # scSLIDE, step by step

    A single-cell experiment measures thousands of cells. Yet many questions
    are about **samples**: patients, embryos, donor lines, or the guide RNAs
    of a CRISPR screen. Which samples resemble each other? Can they be put in
    order, from early to late, or from mild to severe?

    **scSLIDE** is a method built for these questions. It was published by
    Jiang and colleagues
    ([Cell, in press](https://doi.org/10.1016/j.cell.2026.10.018)) and is
    part of version 6 of the Seurat software. This page does three things.

    1. It explains each step with a small example you can play with.
    2. It repeats one analysis of the paper in Python, on a CRISPR screen,
       and checks the result against the paper and against Seurat.
    3. It says what the method means for donor lines grown in a cell village.

    You do not need to write code. Move a slider and the page recomputes
    in your browser. The settings in the panel on the left apply to every
    chart of the real data.
    """)
    return


@app.cell(hide_code=True)
def _():
    _steps = [
        (
            "1. Two views",
            "Describe every cell twice: once by what kind of cell it is, once by what the experiment did to it.",
        ),
        (
            "2. Join them",
            "Combine the two descriptions so that both kinds of difference count.",
        ),
        (
            "3. Landmarks",
            "Pick a few thousand reference cells. Count how many cells of each sample sit next to each one.",
        ),
        ("4. Normalise", "Ask where a sample has more cells than expected, or fewer."),
        (
            "5. Compare",
            "Samples with similar patterns are similar. Then cluster them or put them in order.",
        ),
    ]
    mo.vstack(
        [
            mo.md("## The method in one row"),
            mo.hstack(
                [
                    mo.md(f"**{title}**\n\n{text}").callout(kind="neutral")
                    for title, text in _steps
                ],
                widths="equal",
                align="stretch",
                wrap=True,
            ),
            mo.md(
                """
                Optimal transport, which you may have seen mentioned next to
                scSLIDE, is **not** one of these steps. It is what a competing
                method uses for step 5. We come to it there.
                """
            ),
        ]
    )
    return


@app.cell(hide_code=True)
def _():
    _archive = np.load(io.BytesIO(read_public("cells.npz")))
    cells = {name: _archive[name] for name in _archive.files}
    guides = pd.read_csv(io.BytesIO(read_public("sgrna.csv")))
    manifest = json.loads(read_public("manifest.json"))
    validation = json.loads(read_public("validation.json"))
    reading_guide = read_public("reading_guide.md").decode("utf-8")
    return cells, guides, manifest, reading_guide, validation


@app.cell(hide_code=True)
def _(cells, guides):
    views = [cells["pca"].astype(np.float32), cells["pls"].astype(np.float32)]
    guide_of_cell = cells["guide"].astype(np.int64)
    gene_of_guide = guides["gene"].to_numpy(dtype=object)
    label_of_guide = guides["label"].to_numpy(dtype=object)
    is_control = guides["control"].to_numpy(dtype=bool)
    activity = guides["relative_activity_day10"].to_numpy(dtype=float)
    target_genes = [str(name) for name in cells["target_genes"]]
    return (
        activity,
        gene_of_guide,
        guide_of_cell,
        is_control,
        label_of_guide,
        target_genes,
        views,
    )


@app.cell(hide_code=True)
def _(target_genes):
    focus_gene = mo.ui.dropdown(target_genes, value="RPS15", label="Gene to highlight")
    view_weight = mo.ui.slider(
        0.0,
        1.0,
        value=0.5,
        step=0.1,
        show_value=True,
        label="Weight of the label-aware view",
    )
    landmark_count = mo.ui.dropdown(
        {"1,000": 1000, "2,000": 2000, "5,000 (paper)": 5000},
        value="5,000 (paper)",
        label="Landmark cells",
    )
    votes_per_cell = mo.ui.dropdown(
        {"5 (paper)": 5, "20": 20}, value="5 (paper)", label="Landmarks per cell"
    )
    landmarks_vote = mo.ui.checkbox(False, label="Landmark cells vote too")
    return (
        focus_gene,
        landmark_count,
        landmarks_vote,
        view_weight,
        votes_per_cell,
    )


@app.cell(hide_code=True)
def _(focus_gene, landmark_count, landmarks_vote, view_weight, votes_per_cell):
    mo.sidebar(
        [
            mo.md("### Settings for the real data"),
            focus_gene,
            view_weight,
            landmark_count,
            votes_per_cell,
            landmarks_vote,
            mo.md(
                """
                The paper's settings are a weight of 0.5, 5,000 landmarks and
                5 landmarks per cell. Each setting is explained at its step.
                """
            ),
        ]
    )
    return


@app.cell(hide_code=True)
def _(guides, manifest):
    _summary = manifest["summary"]
    _strip = (
        alt.Chart(guides[~guides["control"]])
        .mark_circle(size=70, opacity=0.8, color=OKABE_ITO["blue"])
        .encode(
            x=alt.X(
                "relative_activity_day10:Q",
                title="Activity of the guide, relative to the best guide for its gene",
            ),
            y=alt.Y("gene:N", title=None, sort="ascending"),
            tooltip=[
                alt.Tooltip("label:N", title="Guide"),
                alt.Tooltip(
                    "relative_activity_day10:Q", title="Relative activity", format=".2f"
                ),
                alt.Tooltip("n_cells:Q", title="Cells"),
            ],
        )
        .properties(width=520, height=420)
    )
    mo.vstack(
        [
            mo.md(
                f"""
                ## The experiment we replicate

                The paper applies scSLIDE to several datasets. We use the smallest:
                a CRISPR screen read out in single cells, from
                [Jost and colleagues (2020)](https://doi.org/10.1038/s41587-019-0387-5).
                It uses CRISPR interference, which turns a gene down instead of
                cutting it. Each of {_summary["target_genes"]} genes is targeted
                by five or six guide RNAs. One matches the gene perfectly. The
                others carry deliberate mismatches, so they turn the gene down
                less, or not at all. Ten more guides target nothing and serve
                as controls.

                In scSLIDE's language, **each guide is a sample**: the group of
                cells that received it. That gives {_summary["guides"]} samples
                and {_summary["cells"]:,} cells, the same numbers as in the paper.
                """
            ),
            mo.hstack(
                [
                    mo.stat(f"{_summary['cells']:,}", label="cells"),
                    mo.stat(str(_summary["guides"]), label="guides, our samples"),
                    mo.stat(str(_summary["target_genes"]), label="target genes"),
                    mo.stat(str(_summary["control_guides"]), label="control guides"),
                ],
                justify="start",
                gap=2,
            ),
            mo.md(
                """
                The original study also measured how well each guide works, in a
                separate experiment based on cell growth. That number, the
                *relative activity*, is 1 for the perfectly matched guide and
                near 0 for a guide that does nothing. scSLIDE never sees it,
                so we can use it later as an independent check.
                """
            ),
            _strip,
        ]
    )
    return


@app.cell(hide_code=True)
def _():
    toy_shift = mo.ui.slider(
        0.0,
        3.0,
        value=0.0,
        step=0.25,
        show_value=True,
        label="Change in state: move sample B's type-2 cells",
    )
    toy_share = mo.ui.slider(
        10,
        90,
        value=50,
        step=5,
        show_value=True,
        label="Change in abundance: % of type-2 cells in sample B",
    )
    return toy_share, toy_shift


@app.cell(hide_code=True)
def _(toy_share, toy_shift):
    def _draw(rng, share_two, shift_two, n_cells=300):
        n_two = int(round(n_cells * share_two))
        one = rng.normal([0.0, 0.0], 0.6, size=(n_cells - n_two, 2))
        two = rng.normal([4.0, shift_two], 0.6, size=(n_two, 2))
        return np.vstack([one, two])

    _rng = np.random.default_rng(0)
    _a = _draw(_rng, 0.5, 0.0)
    _b = _draw(_rng, toy_share.value / 100, toy_shift.value)
    _pooled = np.vstack([_a, _b])
    _sample = np.repeat([0, 1], [len(_a), len(_b)])

    _average_gap = float(np.linalg.norm(_a.mean(axis=0) - _b.mean(axis=0)))
    _share_gap = abs(toy_share.value - 50)
    _landmark_rows = np.sort(
        np.random.default_rng(1).choice(len(_pooled), 30, replace=False)
    )
    _nearest, _ = nearest_rows(_pooled, _pooled[_landmark_rows], 5)
    _counts = landmark_counts(_nearest, _sample, 30, 2)
    _profile = _counts / _counts.sum(axis=0)
    _moved = float(0.5 * np.abs(_profile[:, 0] - _profile[:, 1]).sum())

    _cells = pd.DataFrame(_pooled, columns=["x", "y"]).assign(
        sample=np.where(_sample == 0, "Sample A", "Sample B")
    )
    _landmarks = pd.DataFrame(_pooled[_landmark_rows], columns=["x", "y"]).assign(
        difference=_profile[:, 1] - _profile[:, 0]
    )
    _picture = (
        alt.Chart(_cells)
        .mark_circle(size=22, opacity=0.5)
        .encode(
            x=alt.X("x:Q", title="Gene program 1", scale=alt.Scale(domain=[-2.5, 6.5])),
            y=alt.Y("y:Q", title="Gene program 2", scale=alt.Scale(domain=[-2.5, 5.5])),
            color=alt.Color(
                "sample:N",
                title=None,
                scale=alt.Scale(range=[OKABE_ITO["blue"], OKABE_ITO["orange"]]),
            ),
        )
        + alt.Chart(_landmarks)
        .mark_point(
            shape="diamond", size=110, filled=True, stroke="black", strokeWidth=1
        )
        .encode(
            x="x:Q",
            y="y:Q",
            fill=alt.Fill(
                "difference:Q",
                title=["Landmark:", "more B (yellow)", "or more A (blue)"],
                scale=alt.Scale(scheme="cividis", domainMid=0),
            ),
        )
    ).properties(width=430, height=330)

    mo.vstack(
        [
            mo.md(
                """
                ## Why not simply average, or count cell types?

                Two simple ways to describe a sample come to mind. Both lose
                something. Here are two samples, each with two cell types.
                Change sample B in two different ways and watch three summaries.
                """
            ),
            mo.hstack([toy_shift, toy_share], justify="start", gap=2),
            mo.hstack(
                [
                    _picture,
                    mo.vstack(
                        [
                            mo.stat(
                                f"{_average_gap:.2f}",
                                label="Average cell: distance between the two averages",
                            ),
                            mo.stat(
                                f"{_share_gap} points",
                                label="Cell-type counts: difference in % of type 2",
                            ),
                            mo.stat(
                                f"{_moved:.0%}",
                                label="Landmarks: share of cells sitting somewhere else",
                            ),
                        ]
                    ),
                ],
                justify="start",
                gap=2,
            ),
            mo.md(
                """
                **Try this.** Move only the first slider. Cells of type 2 change
                their state, but there are as many of them as before, so the
                cell-type counts stay blind. Put it back and move only the second
                slider. Now the average cell moves too, although no cell changed
                at all: only the mix did. From the average alone you cannot tell
                the two situations apart.

                The diamonds are **landmarks**: reference cells spread over the
                picture. Counting each sample's cells around every landmark sees
                both kinds of change, and shows where they happen. That is the
                heart of scSLIDE. With both sliders at rest the last number is
                still not zero: two samples drawn from the same population
                differ a little by chance.
                """
            ),
        ]
    )
    return


@app.cell(hide_code=True)
def _():
    effect_size = mo.ui.slider(
        0.0,
        3.0,
        value=1.5,
        step=0.25,
        show_value=True,
        label="Size of the treatment effect",
    )
    return (effect_size,)


@app.cell(hide_code=True)
def _(effect_size):
    _rng = np.random.default_rng(3)
    _n = 600
    _cell_type = _rng.integers(0, 2, _n)
    _treated = _rng.integers(0, 2, _n)
    _features = _rng.normal(size=(_n, 12))
    # Six genes mark the cell type; three other genes respond to the treatment.
    _features[:, :6] += 2.5 * _cell_type[:, None]
    _features[:, 8:11] += effect_size.value * _treated[:, None]
    _features = (_features - _features.mean(axis=0)) / _features.std(axis=0)

    _pca = np.linalg.svd(_features, full_matrices=False)
    _pca = _pca[0][:, :2] * _pca[1][:2]
    _pls = canonical_pls(_features, np.eye(2)[_treated], 2)["scores"]

    def _frame(coordinates, name):
        return pd.DataFrame(coordinates, columns=["first", "second"]).assign(
            view=name,
            treatment=np.where(_treated == 1, "Treated", "Untreated"),
            kind=np.where(_cell_type == 1, "Cell type 2", "Cell type 1"),
        )

    def _explained(coordinates, labels):
        return float(np.corrcoef(coordinates[:, 0], labels)[0, 1] ** 2)

    _both = pd.concat(
        [
            _frame(_pca, "Unsupervised view (PCA)"),
            _frame(_pls, "Label-aware view (PLS)"),
        ]
    )
    _picture = (
        alt.Chart(_both)
        .mark_point(size=28, opacity=0.6, filled=True)
        .encode(
            x=alt.X("first:Q", title="First axis"),
            y=alt.Y("second:Q", title="Second axis"),
            color=alt.Color(
                "treatment:N",
                title="Treatment",
                scale=alt.Scale(range=[OKABE_ITO["vermillion"], OKABE_ITO["sky"]]),
            ),
            shape=alt.Shape("kind:N", title="Cell type"),
        )
        .properties(width=290, height=260)
        .facet(
            column=alt.Column("view:N", title=None, sort=["Unsupervised view (PCA)"])
        )
        .resolve_scale(x="independent", y="independent")
    )
    _table = pd.DataFrame(
        {
            "First axis of": ["Unsupervised view (PCA)", "Label-aware view (PLS)"],
            "explained by cell type": [
                f"{_explained(_pca, _cell_type):.0%}",
                f"{_explained(_pls, _cell_type):.0%}",
            ],
            "explained by treatment": [
                f"{_explained(_pca, _treated):.0%}",
                f"{_explained(_pls, _treated):.0%}",
            ],
        }
    )
    mo.vstack(
        [
            mo.md(
                """
                ## Step 1. Two ways to look at a cell

                Counting cells around landmarks needs a map on which cells that
                are alike sit close together. Which map? scSLIDE uses two.

                The **unsupervised view** is what every single-cell analysis
                starts with: principal component analysis (PCA), or a reference
                atlas. It finds the directions in which cells differ most. Those
                are usually cell types, so it is good at keeping cell types apart.
                A subtle effect of disease or treatment, carried by few genes,
                can be invisible in it.

                The **label-aware view** uses partial least squares (PLS). PLS is
                told a label for every cell, such as *treated* or *untreated*,
                and finds the directions that follow that label. It picks up
                subtle effects, but it sees only what it was told to look for
                and it blurs the cell types.

                In this toy, six genes mark the cell type and three other
                genes respond to a treatment.
                """
            ),
            effect_size,
            _picture,
            markdown_table(_table),
            mo.md(
                """
                The unsupervised view separates the two cell types (the two
                shapes) and mixes the colours. The label-aware view pulls the
                colours apart along its first axis. Slide the effect to zero:
                the label-aware view then has nothing real to find, and the
                little it shows is chance.

                Neither view is enough alone, which is why scSLIDE keeps both.
                """
            ),
        ]
    )
    return


@app.cell(hide_code=True)
def _(
    cells,
    focus_gene,
    gene_of_guide,
    guide_of_cell,
    is_control,
    label_of_guide,
):
    _rng = np.random.default_rng(0)
    _gene_of_cell = gene_of_guide[guide_of_cell]
    _in_focus = _gene_of_cell == focus_gene.value
    _is_control_cell = is_control[guide_of_cell]
    _others = np.flatnonzero(~_in_focus & ~_is_control_cell)
    _controls = np.flatnonzero(_is_control_cell)
    _rows = np.concatenate(
        [
            _rng.choice(_others, 4000, replace=False),
            _rng.choice(_controls, 800, replace=False),
            np.flatnonzero(_in_focus),
        ]
    )
    _layer = np.where(
        _in_focus[_rows],
        "focus",
        np.where(_is_control_cell[_rows], "control", "background"),
    )
    _frame = pd.DataFrame(
        {
            "umap_1": cells["umap"][_rows, 0].astype(float),
            "umap_2": cells["umap"][_rows, 1].astype(float),
            "layer": _layer,
            "guide": label_of_guide[guide_of_cell[_rows]],
            "gene": _gene_of_cell[_rows],
        }
    )
    mo.vstack(
        [
            mo.md(
                f"""
                ### On the real data

                The screen has no cell types: all cells are of one kind. Its
                unsupervised view is a PCA of the 2,000 most variable genes. Its
                label-aware view is a PLS that was told the target gene of every
                cell. The map below is drawn from both views together. Each dot
                is a cell. Coloured dots carry a guide against
                **{focus_gene.value}**, black dots a control guide, grey dots
                anything else.
                """
            ),
            highlight_layers(
                _frame,
                "umap_1",
                "umap_2",
                "guide",
                f"Guide against {focus_gene.value}",
                [alt.Tooltip("guide:N", title="Guide")],
            ).properties(width=560, height=440),
            mo.md(
                """
                Look at where each colour sits. Guides that work move their cells
                to a region of their own. Guides that do not work leave their
                cells in the middle, among the controls. Pick another gene in
                the panel on the left to see the same thing again.

                One warning belongs here. The label-aware view was told the
                target gene. That cells with the same target sit together is
                therefore partly built in. Later we rely on checks that the
                method could not have been told.
                """
            ),
        ]
    )
    return


@app.cell(hide_code=True)
def _(cells, landmark_count, views):
    landmark_rows = sample_landmarks(cells["leverage"], landmark_count.value, seed=0)
    affinity = landmark_affinities(views, [view[landmark_rows] for view in views])
    return affinity, landmark_rows


@app.cell(hide_code=True)
def _(
    affinity,
    guide_of_cell,
    guides,
    landmark_rows,
    landmarks_vote,
    votes_per_cell,
):
    def sample_level(weight, keep=None):
        """Run the sample-level steps with the given weight on the label-aware view."""
        return run_sample_level(
            *affinity,
            guide_of_cell,
            len(guides),
            landmark_rows,
            weights=(1.0 - weight, weight),
            k=votes_per_cell.value,
            landmarks_vote=landmarks_vote.value,
            keep=keep,
        )

    pca_only = sample_level(0.0)
    pls_only = sample_level(1.0)
    return pca_only, pls_only, sample_level


@app.cell(hide_code=True)
def _(sample_level, view_weight):
    current = sample_level(view_weight.value)
    return (current,)


@app.cell(hide_code=True)
def _(current, pca_only, pls_only, view_weight, votes_per_cell):
    def _shared(first, second):
        both = np.concatenate([first, second], axis=1)
        both.sort(axis=1)
        return float((both[:, 1:] == both[:, :-1]).sum(axis=1).mean())

    _k = votes_per_cell.value
    _with_pca = _shared(current["landmarks"], pca_only["landmarks"])
    _with_pls = _shared(current["landmarks"], pls_only["landmarks"])
    _between = _shared(pca_only["landmarks"], pls_only["landmarks"])
    mo.vstack(
        [
            mo.md(
                """
                ## Step 2. Joining the two views

                Now each cell has two sets of coordinates. scSLIDE joins them
                with a recipe its authors invented for another purpose: cells
                measured in two ways at once, for example RNA and protein. It is
                called weighted nearest neighbours (WNN).

                The idea is to ask each view the same question: *how close is
                this cell to that landmark?* Each view answers with a score
                between 0 and 1. The score is 1 for the cell's nearest landmark
                and fades with distance, judged against the cell's own
                surroundings in that view. The two scores are then averaged.
                scSLIDE fixes the weights at one half each.
                """
            ),
            mo.accordion(
                {
                    "The formula, for the curious": mo.md(
                        r"""
                        In one view, let $d$ be the distance from a cell to a
                        landmark, $d_1$ the distance to its nearest landmark and
                        $d_{20}$ the distance to its 20th nearest. The score is

                        $$\theta = \exp\!\left(-\,\frac{\max(d - d_1,\,0)}{d_{20} - d_1}\right)$$

                        and the joint score is
                        $0.5\,\theta_{\text{unsupervised}} + 0.5\,\theta_{\text{label-aware}}$.
                        Distances are measured after scaling every cell's
                        coordinates to length one. Each cell keeps the landmarks
                        with the highest joint score.
                        """
                    ),
                    "Why not just glue the two sets of coordinates together?": mo.md(
                        """
                        Because distances in the two views are not on the same
                        scale. One view may spread cells widely, the other pack
                        them tightly, and the wider one would decide everything.
                        The score above is relative: it compares each distance
                        with the cell's own neighbourhood in that view, so both
                        views get an equal say.
                        """
                    ),
                }
            ),
            mo.md(
                f"""
                ### On the real data

                The slider *Weight of the label-aware view* in the left panel
                sets the balance: 0 uses the unsupervised view only, 1 the
                label-aware view only. It is now at **{view_weight.value}**.

                Each cell keeps {_k} landmarks. Left alone, the two views
                disagree a good deal: on average they share only
                **{_between:.1f} of {_k}** landmarks for a cell. With the current
                weight, a cell keeps on average **{_with_pca:.1f}** of the
                landmarks the unsupervised view would choose and
                **{_with_pls:.1f}** of those the label-aware view would choose.
                """
            ),
        ]
    )
    return


@app.cell(hide_code=True)
def _():
    rare_draw = mo.ui.slider(
        1, 20, value=1, step=1, show_value=True, label="Draw the landmarks again"
    )
    return (rare_draw,)


@app.cell(hide_code=True)
def _(rare_draw):
    _rng = np.random.default_rng(11)
    _n_common, _n_rare, _n_landmarks = 1455, 45, 30
    _common = np.abs(_rng.normal(size=(_n_common, 5))) * np.array(
        [1.0, 1.0, 1.0, 0.03, 0.03]
    )
    _rare = np.abs(_rng.normal(size=(_n_rare, 5))) * np.array(
        [0.03, 0.03, 0.03, 1.0, 1.0]
    )
    _matrix = np.vstack([_common, _rare])
    _is_rare = np.repeat([False, True], [_n_common, _n_rare])
    _leverage = exact_leverage_scores(_matrix)
    _draw = np.random.default_rng(rare_draw.value)
    _uniform = np.sort(_draw.choice(len(_matrix), _n_landmarks, replace=False))
    _by_leverage = sample_landmarks(_leverage, _n_landmarks, seed=rare_draw.value)

    def _frame(rows, name):
        frame = pd.DataFrame(
            {
                "common": _matrix[:, 0],
                "rare": _matrix[:, 3],
                "landmark": np.isin(np.arange(len(_matrix)), rows),
                "way": name,
            }
        )
        return frame

    _both = pd.concat(
        [
            _frame(_uniform, "Landmarks drawn at random"),
            _frame(_by_leverage, "Landmarks drawn by leverage score"),
        ]
    )
    _base = alt.Chart(_both).encode(
        x=alt.X("common:Q", title="Activity of a gene program of the common cells"),
        y=alt.Y("rare:Q", title="Activity of a program of the rare cells"),
    )
    _picture = (
        (
            _base.transform_filter(~alt.datum.landmark).mark_circle(
                size=14, color=GREY, opacity=0.5
            )
            + _base.transform_filter(alt.datum.landmark).mark_point(
                shape="diamond",
                size=90,
                filled=True,
                color=OKABE_ITO["vermillion"],
                stroke="black",
                strokeWidth=0.8,
            )
        )
        .properties(width=290, height=240)
        .facet(
            column=alt.Column("way:N", title=None, sort=["Landmarks drawn at random"])
        )
    )
    mo.vstack(
        [
            mo.md(
                """
                ## Step 3. Landmarks, and counting the cells around them

                A **landmark** is an ordinary cell chosen as a reference point.
                scSLIDE uses 5,000 of them. Together they should cover every
                kind of cell in the data, including the rare kinds.

                Picking cells at random does not do that: a population that
                makes up 3% of the cells gets 3% of the landmarks, which may be
                none. scSLIDE picks cells by their **leverage score** instead.
                The score is high for a cell that is unlike most others, so
                unusual cells are picked more often.

                In this toy, 45 of 1,500 cells belong to a rare population.
                Thirty landmarks are drawn, in each of the two ways.
                """
            ),
            rare_draw,
            _picture,
            mo.md(
                f"""
                Drawn at random, **{int(_is_rare[_uniform].sum())}** of the 30
                landmarks fall in the rare population. Drawn by leverage score,
                **{int(_is_rare[_by_leverage].sum())}** do. Draw again a few
                times: the first number stays at zero, one or two, the second
                near ten.
                """
            ),
        ]
    )
    return


@app.cell(hide_code=True)
def _(focus_gene, gene_of_guide, label_of_guide):
    _labels = [
        str(label) for label in label_of_guide[gene_of_guide == focus_gene.value]
    ]
    focus_guide = mo.ui.dropdown(_labels, value=_labels[0], label="Guide")
    return (focus_guide,)


@app.cell(hide_code=True)
def _(
    cells,
    current,
    focus_guide,
    guides,
    label_of_guide,
    landmark_count,
    landmark_rows,
    validation,
    votes_per_cell,
):
    _column = int(np.flatnonzero(label_of_guide == focus_guide.value)[0])
    _frame = pd.DataFrame(
        {
            "umap_1": cells["umap"][landmark_rows, 0].astype(float),
            "umap_2": cells["umap"][landmark_rows, 1].astype(float),
            "votes": current["counts"][:, _column],
        }
    )
    _base = alt.Chart(_frame).encode(
        x=alt.X("umap_1:Q", axis=None), y=alt.Y("umap_2:Q", axis=None)
    )
    _picture = (
        _base.transform_filter(alt.datum.votes == 0).mark_circle(
            size=8, color=GREY, opacity=0.35
        )
        + _base.transform_filter(alt.datum.votes > 0)
        .mark_circle(opacity=0.85, stroke="black", strokeWidth=0.3)
        .encode(
            size=alt.Size("votes:Q", title="Votes", scale=alt.Scale(range=[25, 220])),
            color=alt.Color(
                "votes:Q", title="Votes", scale=alt.Scale(scheme="cividis")
            ),
            tooltip=[alt.Tooltip("votes:Q", title="Votes")],
        )
    ).properties(width=560, height=420)
    _guide = guides.iloc[_column]
    mo.vstack(
        [
            mo.md(
                f"""
                ### Counting, on the real data

                Now the counting. Every cell votes for the
                {votes_per_cell.value} landmarks it is closest to, by the joint
                score of step 2. (The landmark cells themselves are an
                exception, explained under the picture.) Adding the votes of all cells that carry one
                guide tells us where that guide's cells sit.

                Below, each dot is one of the {landmark_count.value:,} landmarks,
                drawn on the same map as before. Coloured dots received votes
                from cells with the chosen guide.
                """
            ),
            focus_guide,
            _picture,
            mo.md(
                f"""
                **{_guide["label"]}** was found in {int(_guide["n_cells"])} cells
                and has a relative activity of
                {_guide["relative_activity_day10"]:.2f}. Compare a guide with
                high activity and one with low activity: the first gathers its
                votes in one region, the second scatters them over the middle.
                """
            ),
            mo.accordion(
                {
                    "Do the landmark cells vote?": mo.md(
                        f"""
                        A landmark is itself a cell, and the landmark nearest
                        to it is itself. The paper's Methods count every cell.
                        Seurat's own table tells a different story: its total is
                        {validation["seurat_count_total"]:,}, which is five votes
                        for each of the {validation["cells"] - validation["landmarks"]:,}
                        cells that are not landmarks. So in Seurat the landmark
                        cells do not vote, and we do the same by default. The
                        tick box in the left panel lets them vote. It changes
                        little.
                        """
                    ),
                    "Does the leverage score matter here?": mo.md(
                        f"""
                        Little. This screen has no rare cell type to miss.
                        Drawing the same number of landmarks at random gives a
                        guide-by-guide result that correlates
                        {validation["uniform_vs_leverage_landmarks"]:.2f} with
                        the one from leverage scores. Two draws by leverage
                        score correlate {np.mean(validation["python_vs_python"]):.2f}
                        with each other. In tissue, where rare cell types exist,
                        the toy above is the better guide.
                        """
                    ),
                }
            ),
        ]
    )
    return


@app.cell(hide_code=True)
def _():
    size_ratio = mo.ui.slider(
        1,
        10,
        value=4,
        step=1,
        show_value=True,
        label="Sample B has this many times the cells of sample A",
    )
    return (size_ratio,)


@app.cell(hide_code=True)
def _(size_ratio):
    _names = ["1", "2", "3", "4"]
    _same = np.array([0.40, 0.30, 0.20, 0.10])
    _other = np.array([0.10, 0.20, 0.30, 0.40])
    _counts = np.column_stack(
        [200 * _same, 200 * size_ratio.value * _same, 200 * _other]
    ).round()
    _residuals = chi_square_residuals(_counts)
    _similarity = column_correlation(_residuals)
    _raw = pd.DataFrame(
        _counts.astype(int), index=_names, columns=["Sample A", "Sample B", "Sample C"]
    )
    _normalised = pd.DataFrame(
        _residuals.round(1), index=_names, columns=["Sample A", "Sample B", "Sample C"]
    )
    mo.vstack(
        [
            mo.md(
                """
                ## Step 4. Making samples comparable

                The raw counts cannot be compared yet. A sample with more cells
                has higher counts everywhere, and a landmark in a crowded region
                collects more votes from every sample.

                scSLIDE removes both effects with one classic formula. For every
                landmark and sample it works out the count it would **expect**
                if all samples spread their cells in the same way, then keeps
                the difference between observed and expected, scaled by the
                square root of the expected count. Statisticians call this a
                chi-square residual. A positive value means *more cells here
                than expected*, a negative one *fewer*.

                In this toy, samples A and B spread their cells in the same way,
                but B is larger. Sample C spreads its cells differently.
                """
            ),
            size_ratio,
            mo.hstack(
                [
                    mo.vstack(
                        [
                            mo.md("**Votes counted**"),
                            markdown_table(_raw.reset_index(names="Landmark")),
                        ]
                    ),
                    mo.vstack(
                        [
                            mo.md("**After normalising**"),
                            markdown_table(_normalised.reset_index(names="Landmark")),
                        ]
                    ),
                ],
                widths="equal",
                gap=2,
            ),
            mo.md(
                f"""
                In the raw counts, B looks nothing like A as soon as it is
                larger. After normalising, A and B rise and fall together
                (correlation **{_similarity[0, 1]:.2f}**) and C goes the other
                way (correlation with A: **{_similarity[0, 2]:.2f}**), whatever
                the size of B.
                """
            ),
        ]
    )
    return


@app.cell(hide_code=True)
def _():
    ot_shift = mo.ui.slider(
        0, 20, value=3, step=1, show_value=True, label="How far sample B has moved"
    )
    return (ot_shift,)


@app.cell(hide_code=True)
def _(ot_shift):
    _states = np.arange(40)

    def _bump(centre):
        profile = np.exp(-0.5 * ((_states - centre) / 2.0) ** 2)
        return profile / profile.sum()

    _a = _bump(9)
    _shifts = np.arange(0, 21)
    _curves = pd.DataFrame(
        {
            "shift": _shifts,
            "Correlation distance": [
                1.0 - np.corrcoef(_a, _bump(9 + shift))[0, 1] for shift in _shifts
            ],
            "Earth mover's distance": [
                wasserstein_on_a_line(_a, _bump(9 + shift)) for shift in _shifts
            ],
        }
    )
    for _name in ("Correlation distance", "Earth mover's distance"):
        _curves[_name] = _curves[_name] / _curves[_name].max()
    _long = _curves.melt("shift", var_name="measure", value_name="distance")
    _profiles = pd.DataFrame(
        {
            "state": np.tile(_states, 2),
            "share": np.concatenate([_a, _bump(9 + ot_shift.value)]),
            "sample": np.repeat(["Sample A", "Sample B"], len(_states)),
        }
    )
    _two_piles = (
        alt.Chart(_profiles)
        .mark_area(opacity=0.55, interpolate="monotone")
        .encode(
            x=alt.X("state:Q", title="Cell state, from early to late"),
            y=alt.Y("share:Q", title="Share of the sample's cells", stack=None),
            color=alt.Color(
                "sample:N",
                title=None,
                scale=alt.Scale(range=[OKABE_ITO["blue"], OKABE_ITO["orange"]]),
            ),
        )
        .properties(width=300, height=220)
    )
    _measures = (
        alt.Chart(_long)
        .mark_line(strokeWidth=2.5)
        .encode(
            x=alt.X("shift:Q", title="How far sample B has moved"),
            y=alt.Y("distance:Q", title="Distance, as a share of its largest value"),
            color=alt.Color(
                "measure:N",
                title=None,
                scale=alt.Scale(range=[OKABE_ITO["vermillion"], OKABE_ITO["green"]]),
            ),
            strokeDash=alt.StrokeDash("measure:N", legend=None),
        )
        + alt.Chart(pd.DataFrame({"shift": [ot_shift.value]}))
        .mark_rule(color="black", strokeDash=[4, 3])
        .encode(x="shift:Q")
    ).properties(width=300, height=220)
    mo.vstack(
        [
            mo.md(
                """
                ## Step 5. Comparing samples, and what optimal transport would do instead

                After step 4 each sample is a long list of numbers, one per
                landmark. scSLIDE calls two samples similar when their lists
                rise and fall together. For the screen it uses the plain
                correlation between the two lists. For the patient cohorts it
                uses a close relative, the cosine distance.

                Other methods compare samples differently. The best known
                alternative is **optimal transport**, used by a method called
                PILOT. Think of each sample as a pile of sand spread over the
                cell states. The earth mover's distance is the least work
                needed to reshape one pile into the other: how much sand moves,
                times how far.
                """
            ),
            ot_shift,
            mo.hstack([_two_piles, _measures], justify="start", gap=2),
            mo.md(
                """
                Move sample B away from sample A. While the two piles overlap,
                both measures grow. Once they no longer overlap, the correlation
                distance is stuck at its maximum: to it, *different* is all it
                can say. The earth mover's distance keeps growing, because it
                knows how far apart the states are.

                So optimal transport knows something correlation does not. Why
                does scSLIDE do without it? Three reasons, in decreasing order
                of how sure we are.

                - **It is not needed when samples overlap.** Look at the left
                  part of the chart. Neighbouring samples overlap, and there the
                  two measures agree. A chain of overlapping samples is enough
                  to order them all, as step 6 shows.
                - **The landmarks already carry some geometry.** Each cell votes
                  for several nearby landmarks, so nearby landmarks receive
                  similar votes.
                - **The evidence in the paper.** On two patient cohorts the
                  authors compared scSLIDE with 16 other approaches, PILOT
                  among them, and report that scSLIDE predicted the clinical
                  variables best. That is their result on their data, not a law.

                Below, in *What if we had done something else?*, we try optimal
                transport on the screen ourselves.
                """
            ),
        ]
    )
    return


@app.cell(hide_code=True)
def _():
    stage_spread = mo.ui.slider(
        0.0,
        0.15,
        value=0.08,
        step=0.01,
        show_value=True,
        label="How far apart lines of the same day are",
    )
    cells_per_donor = mo.ui.slider(
        50, 1000, value=300, step=50, show_value=True, label="Cells per donor line"
    )
    return cells_per_donor, stage_spread


@app.cell(hide_code=True)
def _(cells_per_donor, stage_spread):
    _rng = np.random.default_rng(5)
    _days = np.repeat(["Day 4", "Day 8", "Day 12"], 15)
    _nominal = np.repeat([0.2, 0.5, 0.8], 15)
    _progress = np.clip(
        _nominal + _rng.normal(scale=stage_spread.value, size=45), 0.0, 1.0
    )
    _states = np.arange(60)[:, None]
    _profile = np.exp(-0.5 * ((_states - (8 + 44 * _progress)[None, :]) / 7.0) ** 2)
    _profile = _profile / _profile.sum(axis=0)
    _counts = np.column_stack(
        [
            _rng.multinomial(cells_per_donor.value, _profile[:, donor])
            for donor in range(45)
        ]
    )
    _residuals = chi_square_residuals(_counts)
    _centred = _residuals - _residuals.mean(axis=1, keepdims=True)
    _components, _eigenvalues = diffusion_map(cosine_distance(_centred), n_components=2)
    # An eigenvalue of one means the lines fall into groups with no link between them.
    _connected = bool(_eigenvalues[0] < 1.0 - 1e-8)
    _position, _ = principal_curve(_components, start=_components[:, 0])
    # A trajectory has no built-in direction; the day labels give it one.
    if np.corrcoef(_position, _nominal)[0, 1] < 0:
        _position = 1.0 - _position
    _agreement = spearman(_position, _progress)
    _frame = pd.DataFrame(
        {
            "true": _progress,
            "inferred": _position,
            "day": _days,
            "donor": [f"Line {i + 1}" for i in range(45)],
        }
    )
    _scatter = (
        alt.Chart(_frame)
        .mark_circle(size=90, opacity=0.85)
        .encode(
            x=alt.X(
                "true:Q", title="True progress of the line (hidden from the method)"
            ),
            y=alt.Y("inferred:Q", title="Position found by the method"),
            color=alt.Color(
                "day:N",
                title="Day collected",
                sort=["Day 4", "Day 8", "Day 12"],
                scale=alt.Scale(
                    range=[OKABE_ITO["sky"], OKABE_ITO["orange"], OKABE_ITO["purple"]]
                ),
            ),
            tooltip=[
                "donor",
                "day",
                alt.Tooltip("true:Q", format=".2f"),
                alt.Tooltip("inferred:Q", format=".2f"),
            ],
        )
        .properties(width=330, height=300)
    )
    _order = np.argsort(_position)
    _heat = pd.DataFrame(
        {
            "rank": np.repeat(np.arange(45), 60),
            "state": np.tile(np.arange(60), 45),
            "share": (_counts[:, _order] / _counts[:, _order].sum(axis=0)).T.ravel(),
        }
    )
    _strip = (
        alt.Chart(_heat)
        .mark_rect()
        .encode(
            x=alt.X(
                "state:O",
                title="Cell state, from early to late",
                axis=alt.Axis(labels=False, ticks=False),
            ),
            y=alt.Y(
                "rank:O",
                title="Lines, in the order found",
                axis=alt.Axis(labels=False, ticks=False),
            ),
            color=alt.Color(
                "share:Q", title="Share of cells", scale=alt.Scale(scheme="cividis")
            ),
        )
        .properties(width=300, height=300)
    )
    _by_day = _frame.groupby("day")["inferred"].agg(
        lambda values: values.max() - values.min()
    )
    mo.vstack(
        [
            mo.md(
                """
                ## Step 6. Putting samples in order

                *This section is an illustration on made-up data. The paper does
                not order the guides of the screen, so there is nothing to
                replicate here.*

                Often the interesting thing about samples is not which groups
                they form but where they sit along a progression: of a disease,
                or of development. The paper does this for patients with
                COVID-19 and Alzheimer's disease and for zebrafish embryos.

                Imagine a cell village: 45 donor lines differentiating together,
                collected on three days. Lines collected on the same day are not
                all equally far along. Some run ahead and some lag. The method
                is given the cells and nothing else: not the true progress, and
                not even which day came first.
                """
            ),
            mo.hstack([stage_spread, cells_per_donor], justify="start", gap=2),
            mo.hstack([_scatter, _strip], justify="start", gap=2),
            mo.md(
                """
                **At these settings the lines fall into separate groups, with
                no line in between to connect them. The positions shown are
                not an order.** Increase the spread.
                """
            ).callout(kind="warn")
            if not _connected
            else mo.md(""),
            mo.md(
                f"""
                Each dot is a donor line. The found position follows the true
                progress closely (rank correlation **{_agreement:.2f}**), and
                lines of the same day are spread out along it: by
                {_by_day["Day 4"]:.0%}, {_by_day["Day 8"]:.0%} and
                {_by_day["Day 12"]:.0%} of the whole trajectory for the three
                days. The day label alone would have hidden that. The paper
                calls this a *pseudostage* and finds it among zebrafish embryos
                collected at the same hour.

                Give each line fewer cells and the order gets noisier. Set the
                spread to zero and the three days fall into three separate
                groups: with no line in between, there is nothing to connect
                them, and no order can be found. An ordering needs samples that
                overlap.
                """
            ),
        ]
    )
    return


@app.cell(hide_code=True)
def _():
    bend = mo.ui.slider(
        1.0,
        2.0,
        value=1.6,
        step=0.1,
        show_value=True,
        label="How far the path bends back on itself",
    )
    return (bend,)


@app.cell(hide_code=True)
def _(bend):
    _rng = np.random.default_rng(2)
    _angle = np.sort(_rng.uniform(0.0, bend.value * np.pi, 70))
    _points = np.column_stack([np.cos(_angle), np.sin(_angle)]) + _rng.normal(
        scale=0.04, size=(70, 2)
    )
    _centred = _points - _points.mean(axis=0)
    _first_pc = _centred @ np.linalg.svd(_centred, full_matrices=False)[2][0]
    _distances = np.linalg.norm(_points[:, None, :] - _points[None, :, :], axis=2)
    _first_dc = diffusion_map(_distances, n_components=1)[0][:, 0]
    _truth = np.arange(70)
    _frame = pd.DataFrame(
        {
            "x": _points[:, 0],
            "y": _points[:, 1],
            "true": _truth,
            "Straight-line axis (PCA)": rankdata(
                _first_pc * np.sign(np.corrcoef(_first_pc, _truth)[0, 1])
            ),
            "Diffusion map": rankdata(
                _first_dc * np.sign(np.corrcoef(_first_dc, _truth)[0, 1])
            ),
        }
    )
    _path = (
        alt.Chart(_frame)
        .mark_circle(size=60)
        .encode(
            x=alt.X("x:Q", axis=None),
            y=alt.Y("y:Q", axis=None),
            color=alt.Color(
                "true:Q", title="True order", scale=alt.Scale(scheme="cividis")
            ),
        )
        .properties(width=250, height=250, title="Samples, as two numbers each")
    )
    _recovered = (
        alt.Chart(
            _frame.melt(["x", "y", "true"], var_name="method", value_name="found")
        )
        .mark_circle(size=40, opacity=0.8)
        .encode(
            x=alt.X("true:Q", title="True order"),
            y=alt.Y("found:Q", title="Order found"),
            color=alt.Color(
                "method:N",
                title=None,
                scale=alt.Scale(range=[OKABE_ITO["green"], OKABE_ITO["vermillion"]]),
            ),
        )
        .properties(width=330, height=250)
    )
    mo.vstack(
        [
            mo.md(
                """
                ### Why a diffusion map, and not the usual straight-line axes?

                The usual way to reduce many numbers to a few is PCA. It draws
                straight axes through the data. scSLIDE uses a **diffusion
                map** instead. A diffusion map imagines a walker hopping from
                each sample to its closest neighbours, and asks which samples
                the walker reaches quickly from which. It follows the chain of
                neighbours, wherever the chain goes.

                The difference shows when a progression bends back on itself:
                late samples start to resemble early ones in the feature that
                varies most.
                """
            ),
            bend,
            mo.hstack([_path, _recovered], justify="start", gap=2),
            mo.md(
                f"""
                With the current bend, the straight-line axis agrees with the
                true order at a rank correlation of
                **{abs(spearman(_first_pc, _truth)):.2f}**, the diffusion map at
                **{abs(spearman(_first_dc, _truth)):.2f}**.

                At a bend of 1.0 the path is a half circle and both do well.
                As it bends further, the straight-line axis puts the two ends
                next to the middle, while the diffusion map keeps following the
                chain. At 2.0 the circle closes. The last samples are then
                neighbours of the first, and the diffusion map fails as well:
                when early and late samples look alike, no method can tell
                them apart. On the simple village example above, where nothing
                bends, a straight-line axis would have been enough.

                Once samples are laid out by the diffusion map, scSLIDE fits a
                smooth curve through them, called a principal curve, and reads
                each sample's position along it. That is the position plotted
                in the village example.
                """
            ),
        ]
    )
    return


@app.cell(hide_code=True)
def _():
    mo.md("""
    ## Step 7. From an order to genes and to abundance

    *This section only explains. Nothing here is computed on the screen.*

    An order of samples is a starting point. The paper builds three
    analyses on it, all of which treat the **sample**, not the cell, as the
    unit that is counted. That matters: a thousand cells from one donor are
    not a thousand independent observations.

    - **Which genes change along the order?** For one cell type at a time,
      the cells of each sample are added up into one profile per sample.
      A regression then asks, gene by gene, whether expression rises or
      falls with the sample's position. Because the position is a number
      and not a yes-or-no label, a gradual change is easier to detect. In
      the paper's Alzheimer's data, 79% of the genes found this way in one
      cohort were confirmed in a second cohort, against 18% for the usual
      case-control comparison.
    - **Which cell types become more or less abundant?** The paper uses a
      model called Hooke, which is built for counts of cells per type and
      allows for covariates such as age and sex.
    - **Which cell types are affected most overall?** The paper uses a
      method called TRADE, which summarises the changes of all genes in a
      cell type as one number.

    The same three questions apply to a cell village. They are the subject
    of a later notebook, on a dataset where the paper reports them.
    """)
    return


@app.cell(hide_code=True)
def _(activity, current, gene_of_guide, is_control):
    current_distance = 1.0 - current["correlation"]
    current_scores = score_distances(
        current_distance, gene_of_guide, is_control, activity
    )
    return (current_scores,)


@app.cell(hide_code=True)
def _(
    activity,
    current,
    focus_gene,
    gene_of_guide,
    guides,
    is_control,
    label_of_guide,
):
    _residuals = current["residuals"]
    _centred = (_residuals - _residuals.mean(axis=1, keepdims=True)).T
    _left, _values, _ = np.linalg.svd(_centred, full_matrices=False)
    _reduced = _left[:, :30] * _values[:30]
    _map = TSNE(
        n_components=2, perplexity=20, init="pca", random_state=0
    ).fit_transform(_reduced)
    _layer = np.where(
        gene_of_guide == focus_gene.value,
        "focus",
        np.where(is_control, "control", "background"),
    )
    _frame = pd.DataFrame(
        {
            "tsne_1": _map[:, 0],
            "tsne_2": _map[:, 1],
            "layer": _layer,
            "guide": label_of_guide,
            "gene": gene_of_guide,
            "activity": activity,
            "cells": guides["n_cells"].to_numpy(),
        }
    )
    _tooltip = [
        alt.Tooltip("guide:N", title="Guide"),
        alt.Tooltip("activity:Q", title="Relative activity", format=".2f"),
        alt.Tooltip("cells:Q", title="Cells"),
    ]
    mo.vstack(
        [
            mo.md(
                f"""
                ## The replication

                The paper shows the screen in its Supplementary Figure 14. We
                rebuild its panels here from the public data. The figure itself
                is in the paper's
                [supplement](https://doi.org/10.1016/j.cell.2026.10.018).

                ### Every guide as one dot

                Each dot below is a guide: a whole sample, summarised by its
                list of numbers from step 4. Guides with similar lists are
                placed close together. Black dots are the control guides,
                coloured dots the guides against **{focus_gene.value}**.
                """
            ),
            highlight_layers(
                _frame,
                "tsne_1",
                "tsne_2",
                "guide",
                f"Guide against {focus_gene.value}",
                _tooltip,
                size=110,
            ).properties(width=560, height=420),
            mo.md(
                """
                As in the paper's panel b, guides against the same gene sit
                together and the controls form a group of their own in the
                middle. Guides that do not work are found in the middle too.
                Point at a dot to see which guide it is. The exact layout of
                such a map changes from run to run; only who sits next to whom
                carries meaning.
                """
            ),
        ]
    )
    return


@app.cell(hide_code=True)
def _(current, gene_of_guide, is_control, label_of_guide):
    _order, _gene_order = cluster_then_order(
        current["correlation"], gene_of_guide, is_control
    )
    _labels = [str(label) for label in label_of_guide[_order]]
    _genes = gene_of_guide[_order]
    _n = len(_order)
    _matrix = current["correlation"][np.ix_(_order, _order)]
    _frame = pd.DataFrame(
        {
            "row": np.repeat(_labels, _n),
            "column": np.tile(_labels, _n),
            "r": _matrix.ravel(),
        }
    )
    # Name each gene once, at the middle one of its guides.
    _names = {}
    for _gene in _gene_order:
        _rows = np.flatnonzero(_genes == _gene)
        _names[_labels[_rows[len(_rows) // 2]]] = (
            "controls" if _gene == "non-targeting" else str(_gene)
        )
    _label_expression = f"({json.dumps(_names)})[datum.value] || ''"
    _axis = alt.Axis(
        labelExpr=_label_expression,
        ticks=False,
        labelFontSize=9,
        title=None,
        labelLimit=80,
    )
    _heat_map = (
        alt.Chart(_frame)
        .mark_rect()
        .encode(
            x=alt.X("column:O", sort=_labels, axis=_axis.copy()),
            y=alt.Y("row:O", sort=_labels, axis=_axis.copy()),
            color=alt.Color(
                "r:Q",
                title="Correlation",
                scale=alt.Scale(scheme="cividis", domain=[-0.2, 1.0], clamp=True),
            ),
            tooltip=[
                alt.Tooltip("row:N", title="Guide"),
                alt.Tooltip("column:N", title="and guide"),
                alt.Tooltip("r:Q", title="Correlation", format=".2f"),
            ],
        )
        .properties(width=600, height=600)
    )
    mo.vstack(
        [
            mo.md(
                """
                ### Every guide against every guide

                This is the paper's panel c. Each row and each column is a
                guide, and the colour says how alike two guides are: yellow for
                alike, dark blue for unrelated. Guides are grouped by gene, and
                genes that resemble each other are placed side by side. The
                paper draws the same picture in red and blue.
                """
            ),
            _heat_map,
        ]
    )
    return


@app.cell(hide_code=True)
def _(activity, current, gene_of_guide, is_control):
    _correlation = current["correlation"]
    _active = ~is_control & (np.nan_to_num(activity) >= 0.5)

    def _block(first, second):
        rows = np.flatnonzero(_active & (gene_of_guide == first))
        columns = np.flatnonzero(_active & (gene_of_guide == second))
        values = _correlation[np.ix_(rows, columns)]
        if first == second:
            values = values[np.triu_indices(len(rows), 1)]
        return float(values.mean()) if values.size else float("nan")

    _genes = sorted(set(gene_of_guide[~is_control]))
    _pairs = sorted(
        (
            (_block(first, second), first, second)
            for index, first in enumerate(_genes)
            for second in _genes[index + 1 :]
        ),
        reverse=True,
    )
    _controls = np.flatnonzero(is_control)
    _control_block = float(
        _correlation[np.ix_(_controls, _controls)][
            np.triu_indices(len(_controls), 1)
        ].mean()
    )
    _same_gene = float(np.nanmean([_block(gene, gene) for gene in _genes]))
    _ribosomal = float(
        np.mean(
            [
                _block("RPS14", "RPS15"),
                _block("RPS14", "RPS18"),
                _block("RPS15", "RPS18"),
            ]
        )
    )
    _table = pd.DataFrame(
        {
            "Pair of genes": [
                f"{first} and {second}" for _, first, second in _pairs[:6]
            ],
            "Average correlation between their guides": [
                f"{value:.2f}" for value, _, _ in _pairs[:6]
            ],
        }
    )
    mo.vstack(
        [
            mo.md(
                f"""
                **What to look for.** The bright squares on the diagonal are
                the genes: guides against the same gene are alike (average
                correlation **{_same_gene:.2f}** among guides with a relative
                activity of at least 0.5). A dark row cutting through a bright
                square is a guide that does not work. The large faint square is
                the controls, which are alike only weakly
                (**{_control_block:.2f}**): they all sit in the same broad
                middle, with nothing specific in common.

                Off the diagonal, a few pairs of *different* genes resemble
                each other. These are the strongest at the current settings:
                """
            ),
            markdown_table(_table),
            mo.md(
                f"""
                The paper's panel shows two such groups clearly: EIF2S1 with
                HSPA9, and the three ribosomal genes RPS14, RPS15 and RPS18.
                (The paper's text names RPL14, RPL15 and RPL18. Those genes
                are not in this screen; the labels in the paper's own figure
                read RPS.) Here the three ribosomal genes correlate
                **{_ribosomal:.2f}** on average.

                **Where we differ from the paper's picture.** With the settings
                the paper states for all its analyses, our squares off the
                diagonal are fainter than in its panel c, and so is the square
                of the controls. Seurat 6, run by us with the same settings,
                gives fainter values too, close to our Python (see *Does the
                Python version match Seurat?* below, which also reports the
                one number where the two differ). So this is not a fault of
                the translation. The paper does not say which settings it used
                for the screen. Two of them change the picture a great deal,
                and you can try both in the left panel.

                - **Fewer landmarks**, or **more landmarks per cell**, give
                  every landmark more votes. The numbers become less noisy and
                  all correlations rise.
                - **More weight on the unsupervised view** brings out the
                  ribosomal group. The label-aware view was told to keep
                  different genes apart, and it does, even where two genes
                  cause the same response. Move the weight to 0 to see it.
                """
            ),
        ]
    )
    return


@app.cell(hide_code=True)
def _(
    activity,
    cells,
    current,
    focus_gene,
    gene_of_guide,
    guide_of_cell,
    guides,
    is_control,
    label_of_guide,
    target_genes,
):
    _rows = np.flatnonzero(gene_of_guide == focus_gene.value)
    _rows = _rows[np.argsort(-activity[_rows])]
    _labels = [str(label) for label in label_of_guide[_rows]]
    _block = current["correlation"][np.ix_(_rows, _rows)]
    _square = pd.DataFrame(
        {
            "row": np.repeat(_labels, len(_rows)),
            "column": np.tile(_labels, len(_rows)),
            "r": _block.ravel(),
        }
    )
    _base = alt.Chart(_square).encode(
        x=alt.X("column:O", sort=_labels, title=None, axis=alt.Axis(labelAngle=-40)),
        y=alt.Y("row:O", sort=_labels, title=None),
    )
    _matrix = (
        _base.mark_rect().encode(
            color=alt.Color(
                "r:Q",
                title="Correlation",
                scale=alt.Scale(scheme="cividis", domain=[-0.2, 1.0], clamp=True),
            )
        )
        + _base.mark_text(fontSize=12).encode(
            text=alt.Text("r:Q", format=".2f"),
            color=alt.condition(
                alt.datum.r > 0.45, alt.value("black"), alt.value("white")
            ),
        )
    ).properties(width=250, height=250, title="How alike the guides are")

    _expression = cells["target_expression"][
        :, target_genes.index(focus_gene.value)
    ].astype(float)
    _cell_gene = gene_of_guide[guide_of_cell]
    _shown = (_cell_gene == focus_gene.value) | is_control[guide_of_cell]
    _group = np.where(
        is_control[guide_of_cell], "Controls", label_of_guide[guide_of_cell]
    )
    _levels = pd.DataFrame({"group": _group[_shown], "expression": _expression[_shown]})
    _boxes = (
        alt.Chart(_levels)
        .mark_boxplot(size=22, outliers=False, color=OKABE_ITO["sky"])
        .encode(
            x=alt.X(
                "group:N",
                sort=[*_labels, "Controls"],
                title=None,
                axis=alt.Axis(labelAngle=-40),
            ),
            y=alt.Y(
                "expression:Q",
                title=f"{focus_gene.value} RNA per cell (normalised, log scale)",
            ),
        )
        .properties(width=250, height=250, title=f"Is {focus_gene.value} turned down?")
    )
    _scores = guides.iloc[_rows].melt(
        id_vars="label",
        value_vars=["relative_activity_day5", "relative_activity_day10"],
        var_name="day",
        value_name="relative activity",
    )
    _scores["day"] = _scores["day"].map(
        {"relative_activity_day5": "Day 5", "relative_activity_day10": "Day 10"}
    )
    _bars = (
        alt.Chart(_scores)
        .mark_bar()
        .encode(
            x=alt.X("label:N", sort=_labels, title=None, axis=alt.Axis(labelAngle=-40)),
            xOffset=alt.XOffset("day:N", sort=["Day 5", "Day 10"]),
            y=alt.Y(
                "relative activity:Q", title="Relative activity (growth experiment)"
            ),
            color=alt.Color(
                "day:N",
                title=None,
                sort=["Day 5", "Day 10"],
                scale=alt.Scale(range=[OKABE_ITO["orange"], OKABE_ITO["blue"]]),
            ),
        )
        .properties(width=250, height=250, title="What the original study measured")
    )
    _controls = np.flatnonzero(is_control)
    _best = _rows[0]
    _summary = pd.DataFrame(
        {
            "Guide": _labels,
            "Relative activity": [f"{value:.2f}" for value in activity[_rows]],
            "Correlation with the best guide": [
                f"{current['correlation'][row, _best]:.2f}" for row in _rows
            ],
            "Average correlation with the controls": [
                f"{current['correlation'][row, _controls].mean():.2f}" for row in _rows
            ],
        }
    )
    mo.vstack(
        [
            mo.md(
                f"""
                ### One gene up close: {focus_gene.value}

                The paper's panels d to k zoom in on four genes: RPS15, MTOR,
                EIF2S1 and BCR. Choose any gene in the left panel. Guides are
                ordered from the most to the least active.
                """
            ),
            mo.hstack([_matrix, _boxes, _bars], justify="start", gap=1, wrap=True),
            markdown_table(_summary),
            mo.md(
                """
                Read the three pictures together. On the left, scSLIDE says
                which guides behave like the best one. In the middle, the cells
                themselves say whether the gene was turned down. On the right,
                an independent experiment says how active each guide is. The
                three agree: a guide that scSLIDE sets apart from the others
                leaves the gene on and scores low in the growth experiment.

                scSLIDE was told which gene each guide targets. It was told
                nothing about which guides work. That it finds them is the part
                of the result that could not have been built in.
                """
            ),
        ]
    )
    return


@app.cell(hide_code=True)
def _(activity, current, current_scores, is_control):
    _controls = np.flatnonzero(is_control)
    _frame = pd.DataFrame(
        {
            "activity": activity[~is_control],
            "near": current["correlation"][:, _controls].mean(axis=1)[~is_control],
        }
    )
    _picture = (
        alt.Chart(_frame)
        .mark_circle(size=55, opacity=0.75, color=OKABE_ITO["blue"])
        .encode(
            x=alt.X(
                "activity:Q", title="Relative activity (independent growth experiment)"
            ),
            y=alt.Y("near:Q", title="Average correlation with the control guides"),
        )
        .properties(width=420, height=300)
    )
    mo.vstack(
        [
            mo.md(
                f"""
                ### All guides at once

                Each dot is one of the 128 guides that target a gene. Guides
                that score low in the growth experiment resemble the controls,
                and guides that score high do not. The rank correlation is
                **{current_scores["activity"]:.2f}**.

                A second number sums up the squares on the diagonal. For each
                guide, take its three most similar guides and let them vote on
                its target gene. They get it right for
                **{current_scores["same_gene"]:.0%}** of the guides.
                """
            ),
            _picture,
        ]
    )
    return


@app.cell(hide_code=True)
def _(
    activity,
    current_scores,
    gene_of_guide,
    guide_of_cell,
    guides,
    is_control,
    pca_only,
    pls_only,
    sample_level,
    view_weight,
    views,
):
    def _score(distances, rows=None):
        if rows is None:
            return score_distances(distances, gene_of_guide, is_control, activity)
        return score_distances(
            distances[np.ix_(rows, rows)],
            gene_of_guide[rows],
            is_control[rows],
            activity[rows],
        )

    _n_guides = len(guides)
    _membership = np.zeros((len(guide_of_cell), _n_guides))
    _membership[np.arange(len(guide_of_cell)), guide_of_cell] = 1.0
    _sizes = _membership.sum(axis=0)

    def _average_cell(view):
        return scaled_cosine_distance((_membership.T @ view) / _sizes[:, None])

    # Twenty groups of cells, as a stand-in for cell types.
    _clustering = KMeans(n_clusters=20, n_init=1, random_state=0).fit(views[0])
    _per_group = np.zeros((_n_guides, 20))
    np.add.at(_per_group, (guide_of_cell, _clustering.labels_), 1.0)
    _shares = _per_group / _per_group.sum(axis=1, keepdims=True)
    _centres = _clustering.cluster_centers_
    _cost = np.linalg.norm(_centres[:, None, :] - _centres[None, :, :], axis=2)

    # Every guide cut down to the same number of cells.
    _rng = np.random.default_rng(0)
    _enough = np.flatnonzero(_sizes >= 50)
    _keep = np.zeros(len(guide_of_cell), dtype=bool)
    for _guide in _enough:
        _keep[
            _rng.choice(np.flatnonzero(guide_of_cell == _guide), 50, replace=False)
        ] = True
    _cut_down = sample_level(view_weight.value, keep=_keep)

    _results = {
        "scSLIDE, with the settings in the left panel": current_scores,
        "scSLIDE, unsupervised view only": _score(1.0 - pca_only["correlation"]),
        "scSLIDE, label-aware view only": _score(1.0 - pls_only["correlation"]),
        "scSLIDE, every guide cut down to 50 cells": _score(
            1.0 - _cut_down["correlation"], _enough
        ),
        "Average cell, unsupervised view": _score(_average_cell(views[0])),
        "Average cell, label-aware view": _score(_average_cell(views[1])),
        "Shares of 20 groups of cells": _score(
            scaled_cosine_distance(np.log1p(1e4 * _shares))
        ),
        "Optimal transport between those shares": _score(
            sinkhorn_distances(_shares, _cost)
        ),
    }
    alternative_scores = _results
    _table = pd.DataFrame(
        {
            "Way of comparing guides": list(_results),
            "Three nearest guides vote for the right gene": [
                f"{row['same_gene']:.0%}" for row in _results.values()
            ],
            "Rank correlation with the growth experiment": [
                f"{row['activity']:.2f}" for row in _results.values()
            ],
        }
    )
    mo.vstack(
        [
            mo.md(
                """
                ## What if we had done something else?

                It is easy to claim that a method is needed. Here we check, on
                this screen, with two yardsticks. The first asks whether guides
                against the same gene end up as neighbours. The second asks
                whether weak guides end up near the controls, judged by the
                growth experiment. For the second, a value further below zero
                is better.
                """
            ),
            markdown_table(_table),
            mo.md(
                """
                Read the table with three things in mind.

                - **The first yardstick is partly built in** wherever the
                  label-aware view is used, because that view was told the
                  target genes. Compare the rows that use it with the rows that
                  do not. The second yardstick is fair to all rows.
                - **Cutting every guide down to 50 cells** weakens the result
                  a little, as fewer cells mean more noise, but does not change
                  it. It is not a side effect of some guides having more cells
                  than others. (Four guides with fewer than 50 cells are left
                  out of that row.)
                - **Counting groups of cells does worst.** These cells are all
                  of one kind. Cutting them into 20 groups and keeping only the
                  shares throws most of the information away. Optimal transport
                  between the shares wins a little of it back, because it knows
                  which groups are close to each other.
                - **The average cell holds up well.** Among the ways that are
                  not told the target gene, the average cell does at least as
                  well as scSLIDE on this screen. That is not a surprise: these
                  are strong perturbations in a single kind of cell, where the
                  average says most of what there is to say. The paper makes
                  the same point, that specialised tools already handle screens
                  in cell lines. Its case for scSLIDE rests on patients and
                  embryos: many cell types, subtle effects, and changes in
                  abundance and in state at once. We have not tested those here.
                """
            ),
        ]
    )
    return (alternative_scores,)


@app.cell(hide_code=True)
def _(validation):
    _rows = validation["what_each_matrix_says"]
    _table = (
        pd.DataFrame(_rows).T.map(lambda value: f"{value:.2f}").reset_index(names="Run")
    )
    _table = _table.rename(
        columns={
            "3-NN accuracy": "Three nearest guides vote for the right gene",
            "Spearman, near controls vs activity": "Rank correlation with the growth experiment",
            "control-control r": "Controls with controls",
            "same-gene r (active guides)": "Same gene",
            "RPS14/15/18 cross r": "RPS14, RPS15, RPS18",
            "EIF2S1-HSPA9 r": "EIF2S1 with HSPA9",
        }
    )
    _swap = validation["with_seurat_landmark_cells"]
    _python = _rows["Python"]["RPS14/15/18 cross r"]
    _seurat = [
        row["RPS14/15/18 cross r"]
        for name, row in _rows.items()
        if name.startswith("Seurat")
    ]
    _sweep = pd.DataFrame(validation["weight_sweep_with_seurat_landmark_cells"])
    _sweep_table = pd.DataFrame(
        {
            "Weight of the label-aware view, in Python": [
                f"{value:.2f}" for value in _sweep["weight"]
            ],
            "RPS14, RPS15, RPS18": [f"{value:.2f}" for value in _sweep["ribosomal"]],
            "Whole table, correlation with Seurat": [
                f"{value:.3f}" for value in _sweep["r_with_seurat"]
            ],
        }
    )
    mo.vstack(
        [
            mo.md(
                f"""
                ## Does the Python version match Seurat?

                The code on this page was written from the Methods section of
                the paper, not copied from Seurat. To check it, we ran
                Seurat {validation["seurat_version"]} itself on the same cells
                with the same settings and compared the end product: the table
                of guide against guide, 138 by 138.

                The comparison needs a yardstick, because landmarks are drawn
                at random and no two runs agree exactly.

                | Compared | Correlation between the two tables |
                |---|---|
                | Seurat with Seurat, different draws of landmarks | {np.mean(validation["seurat_vs_seurat"]):.3f} |
                | Python with Python, different draws of landmarks | {np.mean(validation["python_vs_python"]):.3f} |
                | **Python with Seurat** | **{validation["python_vs_seurat"]["mean"]:.3f}** ({validation["python_vs_seurat"]["min"]:.3f} to {validation["python_vs_seurat"]["max"]:.3f}) |

                Python differs from Seurat about as much as Seurat differs
                from itself.

                Two further checks say where the small remaining difference
                comes from. The cell coordinates agree exactly: both the
                unsupervised and the label-aware view span the same directions
                in Python and in Seurat (agreement
                {validation["subspace_agreement"]["pca"]:.3f} and
                {validation["subspace_agreement"]["pls"]:.3f}, where 1 is
                identical). And when Python is given Seurat's own choice of
                landmark cells, its table correlates
                {_swap["python_coordinates_landmarks_do_not_vote"]:.3f} with
                Seurat's. So most of what is left is the random draw of
                landmarks.

                The numbers below are from these runs, at the paper's settings.
                """
            ),
            markdown_table(_table),
            mo.md(
                f"""
                **One number does not match.** The three ribosomal genes
                correlate {_python:.2f} with one another in Python and
                {min(_seurat):.2f} to {max(_seurat):.2f} in Seurat, and giving
                Python Seurat's landmark cells does not close the gap. This
                number turns out to be the most delicate one on the page. It
                is carried by the unsupervised view and erased by the
                label-aware view, so it swings with the balance between them:
                """
            ),
            markdown_table(_sweep_table),
            mo.md(
                """
                Seurat's value is what Python gives at a weight of about 0.45
                rather than 0.50. In other words, Seurat's way of joining the
                views lets the unsupervised view count a little more than a
                plain average of the two scores does. We did not find out why:
                the Methods describe a plain average, and we did not read
                Seurat's source code. For the screen as a whole the difference
                is small, as the last column shows.

                The scripts for this check are in the repository, under
                `validation/scslide/`.
                """
            ),
        ]
    )
    return


@app.cell(hide_code=True)
def _():
    mo.md("""
    ## What this means for village data

    *This section is our own reading. The paper does not analyse cell
    villages.*

    A cell village is built like this screen. Cells from many sources grow
    in one dish, each cell is assigned to its source afterwards, and the
    sources are what we want to compare. The translation is short.

    | In the screen | In a village |
    |---|---|
    | A guide RNA | A donor line, or a donor line on one day |
    | Cells that carry the guide | Cells assigned to that line |
    | A guide that does not work, sitting with the controls | A line that behaves like the average |
    | Guides against one gene, sitting together | Lines that differentiate alike |
    | The growth experiment, never shown to the method | Anything you know about the lines that the method is not told |

    Five things to decide before running it on a village.

    1. **What is a sample?** A line, or a line on a given day. The second
       lets a line move along the trajectory over time, as in step 6.
    2. **What does the label-aware view get told?** This choice matters
       most. Tell it the day, and it will look for what changes with time.
       Tell it the line, and lines will look distinct by construction, the
       same trap as the target gene here. The screen also showed the
       cost: told to keep genes apart, the label-aware view hid that three
       ribosomal genes act alike. The paper reports that for zebrafish
       development the method worked as well with no labels at all. For a
       strong process like differentiation, try that first.
    3. **Are there enough cells per sample?** The paper drops samples with
       fewer than 200 cells in its patient cohorts. Here the guides have
       about 120 cells each, and 50 were enough for these strong effects. A
       subtle difference between lines will need more.
    4. **Do samples overlap?** An order can only be found through samples
       that resemble their neighbours. Lines collected on days far apart,
       with nothing in between, will come out as separate groups.
    5. **What will you check it against?** Keep something back that the
       method never sees, as the growth experiment was kept back here: a
       marker-based staging, a known difference between lines, or a
       repeat of the village.

    A village has one advantage over a patient cohort. The paper warns
    against pooling samples that were collected or processed differently,
    because such differences change which cells are recovered and look
    like biology. In a village all lines share one dish and one protocol.
    """)
    return


@app.cell(hide_code=True)
def _(reading_guide):
    mo.md(reading_guide)
    return


@app.cell(hide_code=True)
def _(manifest):
    _source = manifest["source"]
    _parameters = manifest["parameters"]
    _summary = manifest["summary"]
    _files = pd.DataFrame(
        {
            "File": [f"`{entry['file']}`" for entry in _source["files"]],
            "Size": [f"{entry['bytes'] / 1e6:.1f} MB" for entry in _source["files"]],
            "SHA-256": [f"`{entry['sha256'][:16]}…`" for entry in _source["files"]],
        }
    )
    mo.vstack(
        [
            mo.md(
                f"""
                ## Where the data and the code come from

                **Data.** {_source["citation"]},
                [doi:{_source["doi"]}](https://doi.org/{_source["doi"]}), from
                the Gene Expression Omnibus, series
                [{_source["accession"]}](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={_source["accession"]}).
                """
            ),
            markdown_table(_files),
            mo.md(
                f"""
                The count matrix is far too large for a browser. A script,
                `prepare.py`, does the heavy steps once and stores a few
                megabytes of results, which this page loads. It keeps the
                {_summary["cells_with_one_assigned_guide"]:,} cells with exactly
                one assigned guide; the paper's thresholds of
                {_parameters["min_genes"]} genes and
                {_parameters["min_umis"]:,} molecules per cell then remove
                {_summary["cells_removed_by_gene_and_umi_thresholds"]} of them.
                It selects {_summary["variable_genes"]:,} variable genes, adds
                {_summary["pls_features"] - _summary["variable_genes"]:,} genes
                that follow the target gene, and computes
                {_parameters["n_pca"]} PCA and {_parameters["n_pls"]} PLS
                coordinates per cell.

                **Method.** Jiang et al., *Reconstructing developmental and
                disease progression with sample-level embeddings*, Cell,
                [doi:10.1016/j.cell.2026.10.018](https://doi.org/10.1016/j.cell.2026.10.018).
                The paper does not state its settings for this screen. Ours are
                the ones it gives for all its analyses, with the target gene as
                the label for the label-aware view.

                **Code.** Written for this page from the paper's Methods, under
                the MIT licence. The charts use the Okabe-Ito colours and the
                cividis scale, which stay readable with colour-blindness.
                """
            ),
        ]
    )
    return


@app.cell(hide_code=True)
def _(alternative_scores, current_scores):
    # Depends on the slowest cells, so it appears once the whole page has run.
    _finished = bool(alternative_scores) and bool(current_scores)
    _where = (
        "in your browser"
        if sys.platform == "emscripten"
        else "on a computer, ahead of time"
    )
    mo.md(
        f"""
        *Everything on this page was computed {_where}, with Python
        {platform.python_version()}.*
        """
        if _finished
        else ""
    )
    return


if __name__ == "__main__":
    app.run()
