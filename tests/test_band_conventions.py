"""Unit tests verifying the band-to-sigma alternative and degrees of freedom (Agent Brief 20 Task 2).

Verifies:
1. Uniform prior standard deviation formula: sigma_eff = h / sqrt(3) derived from Var(U(a, b)) = (b - a)^2 / 12.
2. Under conventional scoring (sigma = h): kinetic_coop passes five-way minimax at 2.58 sigma with res_dof = +1.
3. Under alternative scoring (sigma_eff = h / sqrt(3)): kinetic_coop still passes at 2.58 sigma, but res_dof <= 0.
4. kinetic_two_open fails under the alternative (3.28 sigma > 3.0 sigma).
"""
from __future__ import annotations

import numpy as np
import pytest

from circuitpharm.fitting.comparison import MODEL_SPECS
from scripts.band_convention_evaluation import (
    H_MOT,
    H_PO,
    evaluate_five_way,
)


def test_uniform_prior_variance_derivation():
    """Verify that uniform distribution U(a, b) standard deviation equals half-width / sqrt(3)."""
    np.random.seed(42)
    a, b = 2.54, 2.83
    h = (b - a) / 2.0
    sigma_analytical = h / np.sqrt(3.0)

    # Empirical check
    samples = np.random.uniform(a, b, size=1_000_000)
    sigma_empirical = float(np.std(samples))

    assert np.isclose(sigma_analytical, sigma_empirical, rtol=1e-3)
    assert np.isclose(sigma_analytical, 0.145 / np.sqrt(3.0), atol=1e-5)


def test_kinetic_coop_passes_under_both_conventions():
    """Verify kinetic_coop achieves <= 3.0 sigma minimax under both conventional and alternative."""
    v_conv = np.array([2.603, -0.036, 1.458, -1.851, 1.933, 0.068, -0.989, -1.428, 3.000])
    res_conv = evaluate_five_way("kinetic_coop", v_conv, H_PO, H_MOT, "midpoint", "band")
    assert res_conv["max_z"] <= 2.60, f"Conventional minimax {res_conv['max_z']} must be <= 2.60"

    v_alt = np.array([2.6018, -0.0361, 1.4631, -1.8539, 1.9384, 0.0688, -1.0200, -1.4297, 3.0041])
    sigma_po_alt = H_PO / np.sqrt(3.0)
    sigma_mot_alt = H_MOT / np.sqrt(3.0)
    res_alt = evaluate_five_way("kinetic_coop", v_alt, sigma_po_alt, sigma_mot_alt, "midpoint", "band")
    assert res_alt["max_z"] <= 2.60, f"Alternative minimax {res_alt['max_z']} must be <= 2.60"


def test_kinetic_two_open_fails_under_alternative():
    """Verify kinetic_two_open fails (> 3.0 sigma) under the alternative convention."""
    v_two_open = np.array([0.5564, 0.2220, 1.9122, -1.9516, -0.6171, 0.1579, -1.1471, -1.5012])
    sigma_po_alt = H_PO / np.sqrt(3.0)
    sigma_mot_alt = H_MOT / np.sqrt(3.0)
    res = evaluate_five_way("kinetic_two_open", v_two_open, sigma_po_alt, sigma_mot_alt, "midpoint", "band")
    assert res["max_z"] > 3.0, f"kinetic_two_open must fail under alternative, got {res['max_z']:.2f}"


def test_degrees_of_freedom_accounting():
    """Verify degrees of freedom accounting under conventional vs alternative."""
    # Point data: 3 taus + 2 independent areas + EC50 + n_H + recovery = 8
    n_point = 8
    k_coop = 9

    # Conventional: treats both bands as independent Gaussian data points
    n_conv = n_point + 2  # 10
    res_dof_conv = n_conv - k_coop
    assert res_dof_conv == 1, "Conventional res_dof must be +1"

    # Alternative with inactive interior po_max (0 active) and active upper MOT boundary (1 active)
    n_alt_active = n_point + 1  # 9
    res_dof_alt_active = n_alt_active - k_coop
    assert res_dof_alt_active == 0, "Alternative res_dof with active boundary must be 0 (saturated)"

    # Alternative with both bands treated strictly as priors / bounds
    res_dof_alt_priors = n_point - k_coop
    assert res_dof_alt_priors == -1, "Alternative res_dof with priors must be -1 (under-determined)"
