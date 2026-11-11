"""The scSLIDE method functions, checked on small synthetic inputs."""

from __future__ import annotations

import numpy as np
import pytest
from conftest import NOTEBOOKS_DIR, load_module
from scipy import sparse


@pytest.fixture(scope="module")
def scslide():
    return load_module(NOTEBOOKS_DIR / "scslide" / "scslide.py", "notebook_scslide")


def test_l2_normalize_gives_unit_rows_and_keeps_zero_rows(scslide):
    normalised = scslide.l2_normalize(np.array([[3.0, 4.0], [0.0, 0.0]]))

    assert normalised[0].tolist() == pytest.approx([0.6, 0.8])
    assert normalised[1].tolist() == [0.0, 0.0]


def test_nearest_rows_matches_brute_force(scslide):
    rng = np.random.default_rng(1)
    queries = rng.normal(size=(57, 6))
    references = rng.normal(size=(40, 6))

    indices, distances = scslide.nearest_rows(queries, references, k=7, chunk_size=16)

    full = np.linalg.norm(queries[:, None, :] - references[None, :, :], axis=2)
    expected = np.argsort(full, axis=1)[:, :7]
    assert np.array_equal(indices, expected)
    assert distances == pytest.approx(
        np.take_along_axis(full, expected, axis=1), abs=1e-4
    )


def test_nearest_rows_caps_k_at_the_number_of_references(scslide):
    indices, _ = scslide.nearest_rows(np.eye(3), np.eye(3)[:2], k=10)

    assert indices.shape == (3, 2)


def test_landmark_affinity_follows_the_kernel_formula(scslide):
    # One cell at angle 0; landmarks at 10, 20 and 40 degrees on the unit circle.
    angles = np.deg2rad([10.0, 20.0, 40.0])
    cell = np.array([[1.0, 0.0]])
    landmarks = np.column_stack([np.cos(angles), np.sin(angles)])

    pool, affinities, repeated = scslide.landmark_affinities(
        [cell], [landmarks], candidates=3, bandwidth_rank=2
    )

    distance = 2.0 * np.sin(angles / 2.0)
    expected = np.exp(-(distance - distance[0]) / (distance[1] - distance[0]))
    assert pool.tolist() == [[0, 1, 2]]
    assert not repeated.any()
    assert affinities[0, 0] == pytest.approx(expected, rel=1e-4)


def test_top_landmarks_never_returns_a_landmark_twice(scslide):
    rng = np.random.default_rng(2)
    cells_a, cells_b = rng.normal(size=(30, 4)), rng.normal(size=(30, 3))
    chosen = np.arange(0, 30, 3)

    pool, affinities, repeated = scslide.landmark_affinities(
        [cells_a, cells_b],
        [cells_a[chosen], cells_b[chosen]],
        candidates=6,
        bandwidth_rank=3,
    )
    landmarks, joint = scslide.top_landmarks(pool, affinities, repeated, k=4)

    assert repeated.any(), "the two views should share some candidates in this set-up"
    assert all(len(set(row)) == 4 for row in landmarks.tolist())
    assert np.all(np.diff(joint, axis=1) <= 1e-6)
    # A landmark cell is at distance zero from itself in both views.
    assert landmarks[chosen, 0].tolist() == list(range(len(chosen)))
    assert joint[chosen, 0] == pytest.approx(1.0)


def test_a_zero_weight_view_has_no_say(scslide):
    rng = np.random.default_rng(3)
    cells_a, cells_b = rng.normal(size=(40, 5)), rng.normal(size=(40, 5))
    chosen = np.arange(0, 40, 2)
    pool, affinities, repeated = scslide.landmark_affinities(
        [cells_a, cells_b],
        [cells_a[chosen], cells_b[chosen]],
        candidates=20,
        bandwidth_rank=5,
    )

    landmarks, _ = scslide.top_landmarks(
        pool, affinities, repeated, weights=(1.0, 0.0), k=3
    )

    expected, _ = scslide.nearest_rows(
        scslide.l2_normalize(cells_a), scslide.l2_normalize(cells_a[chosen]), 3
    )
    assert np.array_equal(np.sort(landmarks, axis=1), np.sort(expected, axis=1))


def test_landmark_counts_adds_one_per_cell_and_landmark(scslide):
    landmarks = np.array([[0, 1], [1, 2], [1, 2]])
    sample_codes = np.array([0, 0, 1])

    counts = scslide.landmark_counts(
        landmarks, sample_codes, n_landmarks=3, n_samples=2
    )

    assert counts.tolist() == [[1, 0], [2, 1], [1, 1]]
    assert counts.sum(axis=0).tolist() == [4, 2]


def test_chi_square_residuals_are_zero_under_independence(scslide):
    counts = np.outer([10, 30, 60], [1, 2, 5])

    assert scslide.chi_square_residuals(counts) == pytest.approx(0.0)


def test_chi_square_residuals_match_the_formula_and_ignore_empty_rows(scslide):
    counts = np.array([[8.0, 2.0], [2.0, 8.0], [0.0, 0.0]])

    residuals = scslide.chi_square_residuals(counts)

    assert residuals[0].tolist() == pytest.approx([3 / np.sqrt(5), -3 / np.sqrt(5)])
    assert residuals[2].tolist() == [0.0, 0.0]
    assert np.isfinite(residuals).all()


def test_column_correlation_matches_numpy_and_handles_constant_columns(scslide):
    rng = np.random.default_rng(4)
    matrix = np.column_stack([rng.normal(size=(50, 3)), np.ones(50)])

    correlation = scslide.column_correlation(matrix)

    assert correlation[:3, :3] == pytest.approx(
        np.corrcoef(matrix[:, :3], rowvar=False)
    )
    assert correlation[3].tolist() == [0.0, 0.0, 0.0, 0.0]


def test_cosine_distance_is_zero_for_parallel_and_one_for_orthogonal_columns(scslide):
    matrix = np.array([[1.0, 2.0, 0.0], [0.0, 0.0, 3.0]])

    distance = scslide.cosine_distance(matrix)

    assert distance[0, 1] == pytest.approx(0.0)
    assert distance[0, 2] == pytest.approx(1.0)


def test_diffusion_map_orders_points_along_a_curve(scslide):
    # Points along three quarters of a circle: curved, so not a straight line in 2-D.
    angle = np.linspace(0.0, 1.5 * np.pi, 60)
    points = np.column_stack([np.cos(angle), np.sin(angle)])
    distances = np.linalg.norm(points[:, None, :] - points[None, :, :], axis=2)

    components, eigenvalues = scslide.diffusion_map(distances, n_components=2)

    # The first component is flat at the two ends, so compare ranks, not signs of steps.
    first = components[:, 0]
    rank_correlation = np.corrcoef(np.argsort(np.argsort(first)), np.arange(60))[0, 1]
    assert abs(rank_correlation) > 0.999
    assert eigenvalues[0] > eigenvalues[1] > 0
    assert components[np.abs(first).argmax(), 0] > 0


def test_pca_folds_a_curve_that_the_diffusion_map_unrolls(scslide):
    # Four fifths of a circle: the first principal component is a straight
    # line through it, so it puts the two ends next to the middle.
    angle = np.linspace(0.0, 1.6 * np.pi, 80)
    points = np.column_stack([np.cos(angle), np.sin(angle)])
    distances = np.linalg.norm(points[:, None, :] - points[None, :, :], axis=2)

    components, _ = scslide.diffusion_map(distances, n_components=1)

    def rank_correlation(values):
        return abs(np.corrcoef(np.argsort(np.argsort(values)), np.arange(80))[0, 1])

    centred = points - points.mean(axis=0)
    first_pc = centred @ np.linalg.svd(centred, full_matrices=False)[2][0]
    assert rank_correlation(components[:, 0]) > 0.99
    assert rank_correlation(first_pc) < 0.9


def test_diffusion_map_orders_samples_that_only_overlap_with_their_neighbours(scslide):
    # Forty samples along a line of states; each overlaps only the few next to it,
    # so most pairs are equally far apart.
    rng = np.random.default_rng(0)
    progress = np.sort(rng.uniform(0.0, 1.0, 40))
    states = np.arange(60)[:, None]
    profile = np.exp(-0.5 * ((states - (6 + 48 * progress)[None, :]) / 4.0) ** 2)
    counts = np.column_stack(
        [
            rng.multinomial(300, profile[:, donor] / profile[:, donor].sum())
            for donor in range(40)
        ]
    )
    residuals = scslide.chi_square_residuals(counts)
    centred = residuals - residuals.mean(axis=1, keepdims=True)

    components, _ = scslide.diffusion_map(
        scslide.cosine_distance(centred), n_components=1
    )

    assert abs(scslide.spearman(components[:, 0], progress)) > 0.98


def test_spearman_ignores_monotone_rescaling(scslide):
    values = np.array([0.1, 0.5, 0.2, 0.9, 0.7])

    assert scslide.spearman(values, np.exp(3 * values)) == pytest.approx(1.0)
    assert scslide.spearman(values, -values) == pytest.approx(-1.0)


def test_wasserstein_on_a_line_counts_how_far_mass_moves(scslide):
    pile = np.array([0.0, 1.0, 0.0, 0.0, 0.0])

    assert scslide.wasserstein_on_a_line(pile, np.roll(pile, 1)) == pytest.approx(1.0)
    assert scslide.wasserstein_on_a_line(pile, np.roll(pile, 3)) == pytest.approx(3.0)
    assert scslide.wasserstein_on_a_line(pile, 5 * pile) == pytest.approx(0.0)


def test_sinkhorn_distances_grow_with_how_far_mass_moves(scslide):
    positions = np.arange(6, dtype=float)
    cost = np.abs(positions[:, None] - positions[None, :])
    histograms = np.eye(6)[[0, 1, 4]]

    distances = scslide.sinkhorn_distances(histograms, cost, regularisation=0.02)

    assert distances == pytest.approx(distances.T)
    assert np.diag(distances) == pytest.approx(0.0)
    assert distances[0, 1] == pytest.approx(1.0, abs=0.05)
    assert distances[0, 2] == pytest.approx(4.0, abs=0.05)


def test_scaled_cosine_distance_removes_what_all_rows_share(scslide):
    profiles = np.array([[10.0, 1.0], [10.0, 3.0], [12.0, 1.0], [12.0, 3.0]])

    distance = scslide.scaled_cosine_distance(profiles)

    assert distance[0, 3] == pytest.approx(2.0)
    assert distance[0, 0] == pytest.approx(0.0)


def test_score_distances_rewards_a_matrix_that_groups_genes_and_tracks_activity(
    scslide,
):
    gene = np.array(["a", "a", "a", "b", "b", "b", "ctrl", "ctrl"], dtype=object)
    control = gene == "ctrl"
    activity = np.array([1.0, 0.9, 0.1, 0.95, 0.8, 0.2, np.nan, np.nan])
    # Guides sit on a line: gene a on the left, gene b on the right, controls in
    # the middle. A weaker guide sits closer to the middle.
    position = np.array([-10.0, -9.0, -1.0, 10.0, 8.0, 2.0, 0.0, 0.5])
    distances = np.abs(position[:, None] - position[None, :])

    scores = scslide.score_distances(distances, gene, control, activity, k=1)

    assert scores["activity"] == pytest.approx(-1.0)
    assert scores["same_gene"] == pytest.approx(4 / 6)


def test_cluster_then_order_keeps_genes_together_and_puts_controls_last(scslide):
    gene = np.array(["a", "b", "a", "ctrl", "b", "c", "c", "ctrl", "a"], dtype=object)
    control = gene == "ctrl"
    rng = np.random.default_rng(3)
    base = rng.normal(size=(9, 9))
    correlation = np.corrcoef(base)

    order, gene_order = scslide.cluster_then_order(correlation, gene, control)

    assert sorted(order.tolist()) == list(range(9))
    assert gene_order[-1] == "ctrl"
    ordered = gene[order].tolist()
    assert ordered[-2:] == ["ctrl", "ctrl"]
    # Every gene's guides are next to one another.
    assert all(
        ordered.index(name) + ordered.count(name) - 1
        == len(ordered) - 1 - ordered[::-1].index(name)
        for name in set(ordered)
    )


def test_run_sample_level_leaves_landmark_cells_out_unless_asked(scslide):
    rng = np.random.default_rng(4)
    view = rng.normal(size=(60, 4))
    landmark_rows = np.arange(0, 60, 5)
    samples = rng.integers(0, 3, size=60)
    pool, affinities, repeated = scslide.landmark_affinities(
        [view, view],
        [view[landmark_rows], view[landmark_rows]],
        candidates=8,
        bandwidth_rank=4,
    )

    without = scslide.run_sample_level(
        pool, affinities, repeated, samples, 3, landmark_rows, k=3
    )
    with_votes = scslide.run_sample_level(
        pool, affinities, repeated, samples, 3, landmark_rows, k=3, landmarks_vote=True
    )

    assert without["counts"].sum() == 3 * (60 - len(landmark_rows))
    assert with_votes["counts"].sum() == 3 * 60
    assert without["correlation"].shape == (3, 3)


def test_principal_curve_recovers_the_order_along_an_arc(scslide):
    rng = np.random.default_rng(5)
    angle = np.sort(rng.uniform(0.0, np.pi, size=80))
    points = np.column_stack([np.cos(angle), np.sin(angle)]) + rng.normal(
        scale=0.02, size=(80, 2)
    )

    position, curve = scslide.principal_curve(points)

    rank_correlation = np.corrcoef(np.argsort(np.argsort(position)), np.arange(80))[
        0, 1
    ]
    assert abs(rank_correlation) > 0.99
    assert position.min() == 0.0 and position.max() == pytest.approx(1.0)
    assert curve.shape == (100, 2)


def test_principal_curve_follows_an_arch_when_started_from_the_first_coordinate(
    scslide,
):
    # The shape two diffusion components make: the first orders the samples,
    # the second rises and falls, and both are spread about equally.
    rng = np.random.default_rng(6)
    order = np.sort(rng.uniform(-1.0, 1.0, 45))
    points = np.column_stack([order, 1.5 * (order**2 - 0.33)]) + rng.normal(
        scale=0.03, size=(45, 2)
    )

    position, _ = scslide.principal_curve(points, start=points[:, 0])

    assert abs(scslide.spearman(position, order)) > 0.98
    # Points past an end are told apart, not piled up on the last vertex.
    assert len(np.unique(position.round(6))) >= 40


def test_principal_curve_stays_finite_when_points_pile_up(scslide):
    points = np.repeat(np.array([[0.0, 0.0], [1.0, 1.0], [2.0, 0.0]]), 15, axis=0)

    position, curve = scslide.principal_curve(points, start=points[:, 0])

    assert np.isfinite(position).all() and np.isfinite(curve).all()


def test_project_onto_curve_extends_past_both_ends(scslide):
    curve = np.array([[0.0, 0.0], [1.0, 0.0], [2.0, 0.0]])
    points = np.array([[-0.5, 0.2], [0.5, -0.1], [2.5, 0.3]])

    assert scslide.project_onto_curve(points, curve).tolist() == pytest.approx(
        [-0.5, 0.5, 2.5]
    )
    assert scslide.project_onto_curve(
        points, curve, stretch=0.0
    ).tolist() == pytest.approx([0.0, 0.5, 2.0])


def test_canonical_pls_scores_are_uncorrelated_and_reproduced_by_the_projection(
    scslide,
):
    rng = np.random.default_rng(6)
    features = rng.normal(size=(120, 15))
    labels = rng.integers(0, 3, size=120)
    features[:, 0] += 2.0 * labels
    responses = np.eye(3)[labels]

    fit = scslide.canonical_pls(features, responses, n_components=4)

    scores = fit["scores"]
    gram = scores.T @ scores
    assert gram - np.diag(np.diag(gram)) == pytest.approx(0.0, abs=1e-8)
    assert (features - fit["feature_means"]) @ fit["projection"] == pytest.approx(
        scores
    )
    assert np.linalg.norm(fit["weights"], axis=0) == pytest.approx(1.0)


def test_canonical_pls_finds_a_label_direction_that_pca_misses(scslide):
    rng = np.random.default_rng(7)
    labels = np.repeat([0, 1], 150)
    features = rng.normal(size=(300, 20))
    # Five features share one strong factor that carries no label information.
    features[:, :5] = rng.normal(size=(300, 1)) + 0.3 * rng.normal(size=(300, 5))
    features[:, 10] += 1.5 * labels  # one feature, on its own, carries the label
    features = (features - features.mean(axis=0)) / features.std(axis=0)
    responses = np.eye(2)[labels]

    fit = scslide.canonical_pls(features, responses, n_components=1)

    first_pc = np.linalg.svd(features, full_matrices=False)[0][:, 0]
    assert abs(np.corrcoef(fit["scores"][:, 0], labels)[0, 1]) > 0.5
    assert abs(np.corrcoef(first_pc, labels)[0, 1]) < 0.2


def test_canonical_pls_with_one_response_points_along_the_covariance(scslide):
    rng = np.random.default_rng(8)
    features = rng.normal(size=(80, 6))
    response = features @ rng.normal(size=6) + rng.normal(scale=0.1, size=80)

    fit = scslide.canonical_pls(features, response[:, None], n_components=1)

    centred = features - features.mean(axis=0)
    covariance = centred.T @ (response - response.mean())
    assert abs(
        fit["weights"][:, 0] @ covariance / np.linalg.norm(covariance)
    ) == pytest.approx(1.0)


def test_exact_leverage_scores_sum_to_the_rank(scslide):
    rng = np.random.default_rng(9)
    matrix = rng.normal(size=(50, 4))
    matrix = np.column_stack([matrix, matrix[:, 0] + matrix[:, 1]])  # rank stays 4

    scores = scslide.exact_leverage_scores(matrix)

    assert scores.sum() == pytest.approx(4.0)
    assert scores.min() >= 0.0 and scores.max() <= 1.0 + 1e-12


def test_leverage_scores_single_out_a_rare_kind_of_row(scslide):
    rng = np.random.default_rng(10)
    common = rng.normal(size=(990, 8)) * np.array([1.0] * 4 + [0.01] * 4)
    rare = rng.normal(size=(10, 8)) * np.array([0.01] * 4 + [1.0] * 4)

    scores = scslide.exact_leverage_scores(np.vstack([common, rare]))

    assert scores[990:].mean() > 20 * scores[:990].mean()


@pytest.mark.parametrize("as_sparse", [False, True])
def test_sketched_leverage_scores_track_the_exact_ones(scslide, as_sparse):
    rng = np.random.default_rng(11)
    matrix = rng.normal(size=(3000, 12)) * rng.gamma(2.0, size=(3000, 1))
    exact = scslide.exact_leverage_scores(matrix)

    given = sparse.csr_matrix(matrix) if as_sparse else matrix
    approximate = scslide.sketched_leverage_scores(given, projection_size=256, seed=0)

    # They are used as sampling weights, so what matters is that they rank rows alike.
    assert np.corrcoef(exact, approximate)[0, 1] > 0.95


def test_sketched_leverage_scores_approach_the_exact_ones_with_a_larger_sketch(scslide):
    rng = np.random.default_rng(12)
    matrix = rng.normal(size=(3000, 12)) * rng.gamma(2.0, size=(3000, 1))
    exact = scslide.exact_leverage_scores(matrix)

    approximate = scslide.sketched_leverage_scores(
        matrix, sketch_rows=1500, projection_size=512, seed=0
    )

    assert np.corrcoef(exact, approximate)[0, 1] > 0.98
    assert approximate.sum() == pytest.approx(exact.sum(), rel=0.15)


def test_sample_landmarks_is_reproducible_and_skips_zero_scores(scslide):
    scores = np.array([0.0, 1.0, 1.0, 5.0, 0.0, 2.0])

    first = scslide.sample_landmarks(scores, 3, seed=1)
    second = scslide.sample_landmarks(scores, 3, seed=1)
    everything = scslide.sample_landmarks(scores, 10, seed=1)

    assert first.tolist() == second.tolist()
    assert len(set(first.tolist())) == 3
    assert everything.tolist() == [1, 2, 3, 5]


def test_nearest_label_votes_takes_the_majority_and_breaks_ties_by_closeness(scslide):
    # Sample 0 has neighbours 1 (a), 2 (b), 3 (b): majority b.
    # Sample 4 has neighbours 5 (a) then 6 (b), a tie at k=2: the nearer one wins.
    distances = np.full((7, 7), 50.0)
    np.fill_diagonal(distances, 0.0)
    for (i, j), d in {(0, 1): 1, (0, 2): 2, (0, 3): 3, (4, 5): 1, (4, 6): 2}.items():
        distances[i, j] = distances[j, i] = d
    labels = ["x", "a", "b", "b", "y", "a", "b"]

    assert scslide.nearest_label_votes(distances, labels, k=3)[0] == "b"
    assert scslide.nearest_label_votes(distances, labels, k=2)[4] == "a"
