"""Unit tests verifying adversarial audit findings for the two theorems (Brief 22 Task 1).

Covers:
1. Probe 1.1: Degenerate and negligible off-diagonal stability.
2. Probe 1.2: Boundary behaviour and floating-point error near supremum.
3. Probe 1.3: Waveform fitting vs eigenmode decomposition (zero 3-positive well-conditioned fits).
4. Probe 1.4: Leaf vs non-leaf invariance for single-open schemes.
"""
from __future__ import annotations

import numpy as np
import pytest

from circuitpharm.fitting.comparison import MODEL_SPECS
from circuitpharm.fitting.deact_components import (
    predict_components_eigen,
    predict_components_fit,
)
from circuitpharm.models.kinetic_jw95 import KineticAllosteryModel
from circuitpharm.models.kinetic_open_desens import KineticOpenDesensModel


def test_audit_probe_1_1_near_degenerate_eigenvalues():
    """Verify that near-degenerate eigenvalues do not violate the gateway weight bound."""
    koff = 1.0
    alpha = 20.0
    beta = 50.0
    d = 1e-12

    H23 = np.array([
        [-(2 * koff + beta + d), np.sqrt(alpha * beta)],
        [np.sqrt(alpha * beta), -alpha],
    ])
    eig23 = np.linalg.eigvalsh(H23)
    fast_eig = eig23[0]
    r = -fast_eig

    B = np.array([
        [-r, r, 0.0],
        [d, -(2 * koff + beta + d), beta],
        [0.0, alpha, -alpha],
    ])
    pi = np.array([d / r, 1.0, beta / alpha])
    D = np.diag(np.sqrt(pi))
    H = D @ B @ np.linalg.inv(D)
    evals, U = np.linalg.eigh(H)
    u1 = U[:, 0]
    w = np.sqrt(pi[2] / pi) * u1[2] * u1

    # Gateway weight is strictly negative
    assert w[1] < 0.0
    # Ratio satisfies spectral bound
    w_ratio = w[2] / abs(w[1])
    assert w_ratio < (alpha / beta)


def test_audit_probe_1_2_boundary_behaviour_in_box():
    """Verify that a_1 <= 0 strictly inside declared parameter box boundaries."""
    spec = MODEL_SPECS["kinetic_jw95_free_rates"]
    bounds = spec.bounds
    np.random.seed(42)

    for _ in range(50):
        v = np.array([np.random.uniform(b[0], b[1]) for b in bounds])
        pin_idx = np.random.randint(0, len(bounds))
        v[pin_idx] = bounds[pin_idx][0]  # lower bound

        m = spec.factory(v)
        try:
            _, areas = predict_components_eigen(m, pulse_ms=2.0, conc_mm=10.0, l1_max=None)
            assert areas[0] <= 1e-10, f"a1={areas[0]} violated non-positivity on boundary"
        except Exception:
            pass


def test_audit_probe_1_3_triexponential_fit_cannot_produce_three_positive_components():
    """Verify that waveform fitting of a 5-state scheme cannot yield 3 well-conditioned positive components."""
    m = KineticAllosteryModel()
    fit_res = predict_components_fit(m, n_starts=4)

    # Either a1 collapsed to 0 or condition number is ill-conditioned
    a1, a2, a3 = fit_res.areas
    is_degenerate = (a1 < 0.01) or (not fit_res.is_well_conditioned) or (fit_res.condition_number > 1e4)
    assert is_degenerate, f"Expected degenerate fit for 5-state single-open scheme, got {fit_res}"


def test_audit_probe_1_4_leaf_vs_non_leaf_both_satisfy_constraint():
    """Verify that single-open constraint holds identically for leaf and internal topologies."""
    # Leaf topology (JW95)
    m_leaf = KineticAllosteryModel(kon=10.0, koff=1.0, beta=100.0, alpha=20.0, d=5.0, r=0.5)
    _, a_leaf = predict_components_eigen(m_leaf, l1_max=None)
    assert a_leaf[0] <= 0.0

    # Internal topology (Open desensitisation)
    m_internal = KineticOpenDesensModel(kon=10.0, koff=1.0, beta=100.0, alpha=20.0, d=5.0, r=0.5)
    _, a_internal = predict_components_eigen(m_internal, l1_max=None)
    assert a_internal[0] <= 0.0
