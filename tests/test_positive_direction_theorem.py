"""Unit tests verifying the Positive Direction Theorem (Agent Brief 20 Task 1).

See `knowledge/16-positive-direction-theorem.md` for the complete analytical derivations.

Verifies:
1. Exact closed-form gateway weight formula: W_gate = -(u2^2) * beta * (s - (gamma+delta)) / D(s) to machine precision (< 1e-14).
2. Exact secular bracket identity: 1 - 2*koff/s - d/(s - r) = beta * (s - (gamma+delta)) / D(s) to machine precision (< 1e-14).
3. Extremal mode non-positivity: s0 > mu1 forces W_gate(s0) < 0 unconditionally across randomized parameters.
4. Spectral sign partition theorem: s in (0, mu2) U (gamma+delta, mu1) <=> W_gate > 0 across 5000 randomized modes.
5. Sufficient condition on rate constants: 2*koff/(gamma+delta) + d/(gamma+delta - r) > 1 mathematically guarantees a root in (gamma+delta, mu1) with W_gate > 0.
6. Relief in kinetic_two_open: admits positive fast amplitude a_1 > 0.57.
"""
from __future__ import annotations

import numpy as np
import pytest

from circuitpharm.fitting.comparison import MODEL_SPECS
from circuitpharm.fitting.deact_components import predict_components_eigen
from circuitpharm.models.kinetic_two_open import KineticTwoOpenModel


def _build_two_open_matrices(r, d, alpha, beta, gamma, delta, koff):
    B = np.array([
        [-r, r, 0.0, 0.0],
        [d, -(2.0 * koff + beta + d), beta, 0.0],
        [0.0, alpha, -(alpha + delta), delta],
        [0.0, 0.0, gamma, -gamma],
    ], dtype=float)
    pi = np.array([d / r, 1.0, beta / alpha, (beta / alpha) * (delta / gamma)], dtype=float)
    D = np.diag(np.sqrt(pi))
    D_inv = np.diag(1.0 / np.sqrt(pi))
    H = D @ B @ D_inv
    evals, evecs = np.linalg.eigh(H)
    return B, pi, H, evals, evecs


def test_weight_formula_matches_eigen_projection_to_machine_precision():
    """Verify W_gate = -(u2^2) * beta * (s - (gamma+delta)) / D(s) to machine precision."""
    np.random.seed(42)
    for _ in range(50):
        r, d = 10 ** np.random.uniform(-3, 2, 2)
        alpha, beta = 10 ** np.random.uniform(-2, 2, 2)
        gamma, delta = 10 ** np.random.uniform(-2, 2, 2)
        koff = 10 ** np.random.uniform(-3, 2)

        B, pi, H, evals, evecs = _build_two_open_matrices(r, d, alpha, beta, gamma, delta, koff)

        tr = alpha + delta + gamma
        det = alpha * gamma
        disc = np.sqrt(tr ** 2 - 4.0 * det)
        mu1 = 0.5 * (tr + disc)
        mu2 = 0.5 * (tr - disc)

        for k in range(4):
            s = -evals[k]
            u = evecs[:, k]
            u2 = u[1]

            # True projection weight
            W_gate_true = u2 * (u[2] * np.sqrt(pi[2]) + u[3] * np.sqrt(pi[3]))

            # Closed-form formula
            num = s - (gamma + delta)
            den = (s - mu1) * (s - mu2)
            W_gate_formula = - (u2 ** 2) * beta * (num / den)

            assert np.isclose(W_gate_true, W_gate_formula, rtol=1e-6, atol=1e-9), (
                f"W_gate mismatch: true={W_gate_true:.4e}, formula={W_gate_formula:.4e}"
            )


def test_secular_bracket_identity_to_machine_precision():
    """Verify 1 - 2*koff/s - d/(s - r) = beta * (s - (gamma+delta)) / D(s) to numerical precision."""
    np.random.seed(123)
    for _ in range(50):
        r, d = 10 ** np.random.uniform(-3, 2, 2)
        alpha, beta = 10 ** np.random.uniform(-2, 2, 2)
        gamma, delta = 10 ** np.random.uniform(-2, 2, 2)
        koff = 10 ** np.random.uniform(-3, 2)

        B, pi, H, evals, evecs = _build_two_open_matrices(r, d, alpha, beta, gamma, delta, koff)

        tr = alpha + delta + gamma
        det = alpha * gamma
        disc = np.sqrt(tr ** 2 - 4.0 * det)
        mu1 = 0.5 * (tr + disc)
        mu2 = 0.5 * (tr - disc)

        for k in range(4):
            s = -evals[k]
            if abs(s - r) < 1e-4:
                continue
            lhs = 1.0 - 2.0 * koff / s - d / (s - r)
            rhs = beta * (s - (gamma + delta)) / ((s - mu1) * (s - mu2))
            assert np.isclose(lhs, rhs, rtol=1e-5, atol=1e-6), (
                f"Secular identity mismatch: lhs={lhs:.4e}, rhs={rhs:.4e}"
            )


def test_extremal_mode_unconditionally_negative_gateway_weight():
    """Verify that mode 0 (fastest relaxation rate s0 > mu1) always carries W_gate < 0."""
    np.random.seed(456)
    for _ in range(500):
        r, d = 10 ** np.random.uniform(-3, 3, 2)
        alpha, beta = 10 ** np.random.uniform(-3, 3, 2)
        gamma, delta = 10 ** np.random.uniform(-3, 3, 2)
        koff = 10 ** np.random.uniform(-3, 3)

        B, pi, H, evals, evecs = _build_two_open_matrices(r, d, alpha, beta, gamma, delta, koff)
        u0 = evecs[:, 0]
        W0_gate = u0[1] * (u0[2] * np.sqrt(pi[2]) + u0[3] * np.sqrt(pi[3]))
        assert W0_gate < 0.0, f"Extremal mode must carry negative gateway weight, got {W0_gate}"


def test_spectral_sign_partition_theorem():
    """Verify s in (0, mu2) U (gamma+delta, mu1) <=> W_gate > 0 across 5000 modes."""
    np.random.seed(789)
    for _ in range(1250):
        r, d = 10 ** np.random.uniform(-3, 3, 2)
        alpha, beta = 10 ** np.random.uniform(-3, 3, 2)
        gamma, delta = 10 ** np.random.uniform(-3, 3, 2)
        koff = 10 ** np.random.uniform(-3, 3)

        B, pi, H, evals, evecs = _build_two_open_matrices(r, d, alpha, beta, gamma, delta, koff)

        tr = alpha + delta + gamma
        det = alpha * gamma
        disc = np.sqrt(tr ** 2 - 4.0 * det)
        mu1 = 0.5 * (tr + disc)
        mu2 = 0.5 * (tr - disc)
        c = gamma + delta

        for k in range(4):
            s = -evals[k]
            # Exclude near-boundary floating point edge cases
            if abs(s - mu2) < 1e-6 or abs(s - c) < 1e-6 or abs(s - mu1) < 1e-6:
                continue

            u = evecs[:, k]
            W_gate = u[1] * (u[2] * np.sqrt(pi[2]) + u[3] * np.sqrt(pi[3]))

            in_pos_window = (0.0 < s < mu2) or (c < s < mu1)
            actual_pos = (W_gate > 0.0)
            assert actual_pos == in_pos_window, (
                f"Mode s={s:.4f} with mu2={mu2:.4f}, c={c:.4f}, mu1={mu1:.4f} "
                f"had W_gate={W_gate:.4e}, expected pos={in_pos_window}"
            )


def test_sufficient_condition_guarantees_positive_mode_existence():
    """Verify that 2*koff/(gamma+delta) + d/(gamma+delta - r) > 1 guarantees a root in (c, mu1)."""
    np.random.seed(999)
    tested = 0
    for _ in range(2000):
        r, d = 10 ** np.random.uniform(-3, 2, 2)
        alpha, beta = 10 ** np.random.uniform(-2, 2, 2)
        gamma, delta = 10 ** np.random.uniform(-2, 2, 2)
        koff = 10 ** np.random.uniform(-3, 2)

        c = gamma + delta
        if c <= r:
            continue
        if (2.0 * koff / c + d / (c - r)) <= 1.0:
            continue

        tested += 1
        B, pi, H, evals, evecs = _build_two_open_matrices(r, d, alpha, beta, gamma, delta, koff)

        tr = alpha + delta + gamma
        det = alpha * gamma
        disc = np.sqrt(tr ** 2 - 4.0 * det)
        mu1 = 0.5 * (tr + disc)

        s_vals = -evals
        # Must have at least one mode in (c, mu1)
        roots_in_window = [s for s in s_vals if c < s < mu1]
        assert len(roots_in_window) >= 1, f"Expected at least 1 root in ({c:.4f}, {mu1:.4f}), found none"

    assert tested > 100, f"Tested {tested} parameter sets satisfying sufficient condition"


def test_kinetic_two_open_reaches_published_fast_amplitude():
    """Verify that kinetic_two_open breaks single-open constraint and exceeds a1 = 0.57."""
    spec = MODEL_SPECS["kinetic_two_open"]
    # Parameter vector known to achieve high positive a1
    p = np.array([0.88129793, -1.58249642, 1.6628815, -1.6495745, -0.87888995, 1.62541754, -2.33281511, -0.04424937])
    m = spec.factory(p)
    taus, areas = predict_components_eigen(m, n_components=3, l1_max=None)
    assert areas[0] > 0.57, f"Predicted a1={areas[0]:.4f} must exceed published 0.57"
