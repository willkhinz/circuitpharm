"""Unit tests verifying adversarial audit findings for Observables and Bounds (Brief 22 Tasks 2 & 3).

Covers:
1. Task 2: Mean open time band evaluation across M-mode, H-mode, and Full modal span.
2. Task 3: Windowed Hill slope definition vs analytic EC50 derivative.
3. Task 3: Mono-liganded conductance effect on Hill steepness.
4. Task 4: Simplex covariance Mahalanobis floor = exactly (0.57 / 0.04)^2 = 203.0625.
"""
from __future__ import annotations

import numpy as np
import pytest

from circuitpharm.fitting.comparison import MODEL_SPECS
from circuitpharm.fitting.crc_components import extract_crc_parameters
from circuitpharm.fitting.data import (
    KERAMIDAS_H_MODE_BAND,
    KERAMIDAS_M_MODE_BAND,
    KERAMIDAS_MODAL_SPAN,
)
from circuitpharm.fitting.deact_components import get_barberis_amplitude_covariance
from circuitpharm.models.base import hill_slope


def test_audit_task2_mean_open_time_modal_bands():
    """Verify that kinetic_coop's predicted MOT (~3.20 ms) falls comfortably inside full modal span."""
    v_coop = np.array([2.603, -0.036, 1.458, -1.851, 1.933, 0.068, -0.989, -1.428, 3.000])
    m = MODEL_SPECS["kinetic_coop"].factory(v_coop)
    from circuitpharm.fitting.likelihood import predict_mean_open_time
    from circuitpharm.fitting.data import KERAMIDAS2008_MEAN_OPEN_TIME_M

    mot = float(predict_mean_open_time(m, KERAMIDAS2008_MEAN_OPEN_TIME_M))
    assert 2.9 < mot < 3.5, f"Expected MOT ~3.20 ms, got {mot:.2f} ms"

    # M-mode: outside upper edge [2.54, 2.83]
    assert mot > KERAMIDAS_M_MODE_BAND[1]
    z_m_mode = (mot - KERAMIDAS_M_MODE_BAND[1]) / ((KERAMIDAS_M_MODE_BAND[1] - KERAMIDAS_M_MODE_BAND[0]) / 2.0)
    assert z_m_mode < 3.0, f"M-mode excess {z_m_mode:.2f} sigma must be <= 3 sigma"

    # Full modal span [2.54, 8.63]: strictly interior!
    assert KERAMIDAS_MODAL_SPAN[0] <= mot <= KERAMIDAS_MODAL_SPAN[1]


def test_audit_task3_windowed_hill_slope_can_exceed_two():
    """Verify that the windowed Hill slope definition (1-10 uM) can exceed 2.000 for unconstrained models."""
    spec = MODEL_SPECS["kinetic_jw95"]
    # Known vector from hill_slope_ceiling producing high windowed slope
    v = np.array([1.2888524, 3.0, -1.3534856])
    m = spec.factory(v)
    crc = extract_crc_parameters(m)

    # Windowed slope exceeds 2.000
    assert crc.hill_measured > 2.000, f"Expected windowed n_H > 2.000, got {crc.hill_measured:.4f}"


def test_audit_task3_mono_liganded_conductance_depresses_hill_slope():
    """Verify that mono-liganded conductance g1 > 0 depresses the Hill slope and never lifts it above 2.0."""
    # Two sites: K1 = 100, K2 = 0.1 (cooperativity c = 1000)
    concs = np.geomspace(0.01, 1000.0, 300)

    def get_slope(g1):
        num = 2.0 * g1 * 0.1 * concs + concs**2
        den = 100.0 * 0.1 + 2.0 * 0.1 * concs + concs**2
        y = num / den
        return hill_slope(concs, y, window_um=(1.0, 10.0))

    slope_g0 = get_slope(0.0)
    slope_g1 = get_slope(0.2)
    slope_g2 = get_slope(0.5)

    assert slope_g0 > slope_g1 > slope_g2
    assert slope_g1 < 2.000
    assert slope_g2 < 2.000


def test_audit_task4_simplex_covariance_mahalanobis_floor_exact():
    """Verify that the minimum Mahalanobis distance under simplex covariance when a_1 <= 0 is (0.57/0.04)^2."""
    from scipy.optimize import minimize
    cov, cov_pinv, det_sub = get_barberis_amplitude_covariance()
    a_pub = np.array([0.57, 0.23, 0.20])

    def chi2(a):
        d = a - a_pub
        return float(d @ cov_pinv @ d)

    cons = [
        {"type": "eq", "fun": lambda a: a[0] + a[1] + a[2] - 1.0},
        {"type": "ineq", "fun": lambda a: -a[0]},  # a1 <= 0
    ]
    res = minimize(chi2, [0.0, 0.5, 0.5], constraints=cons, method="SLSQP", tol=1e-12)
    expected_floor = (0.57 / 0.04) ** 2  # 203.0625

    assert np.isclose(res.fun, expected_floor, rtol=1e-6)
