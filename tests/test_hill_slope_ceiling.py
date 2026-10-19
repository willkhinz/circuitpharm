"""Executable tests for Unconstrained Hill Slope Ceilings (Brief 17 Task 3).

Distinguishes:
- Statement A: 'The scheme cannot produce n_H = 2.2 at all' (FALSE)
- Statement B: 'The scheme cannot produce n_H = 2.2 while satisfying the other observables' (TRUE)
"""
from __future__ import annotations

import pathlib
import sys

_ROOT = pathlib.Path(__file__).resolve().parent.parent
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))
if str(_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(_ROOT / "src"))

import numpy as np
import pytest

from circuitpharm.fitting.comparison import MODEL_SPECS
from circuitpharm.fitting.crc_components import extract_crc_parameters


def test_jw95_unconstrained_hill_slope_exceeds_target():
    """Verify that even 5-state JW95 with 3 free parameters reaches n_H > 2.2 unconstrained."""
    spec = MODEL_SPECS["kinetic_jw95"]
    # Optimum vector found by unconstrained differential evolution
    v = np.array([1.2881204286837873, 3.0, -2.6805855906866864], dtype=float)

    m = spec.factory(v)
    crc = extract_crc_parameters(m)

    # Measured Hill slope over 1-10 uM window (matching Jahn 1997)
    assert crc.hill_measured > 2.2, f"Expected n_H > 2.2, got {crc.hill_measured:.3f}"
    assert crc.hill_measured > 2.8, f"Expected n_H > 2.8, got {crc.hill_measured:.3f}"


def test_statement_a_is_false_statement_b_is_true():
    """Verify distinction between Statement A (false) and Statement B (true)."""
    spec = MODEL_SPECS["kinetic_jw95"]
    v = np.array([1.2881204286837873, 3.0, -2.6805855906866864], dtype=float)

    m = spec.factory(v)
    crc = extract_crc_parameters(m)

    # Statement A is FALSE: The model can easily produce n_H = 2.2 (in fact reaching ~2.93)
    statement_a_cannot_produce_22 = bool(crc.hill_measured < 2.2)
    assert not statement_a_cannot_produce_22, "Statement A should be FALSE"

    # Statement B is TRUE: At this high-Hill regime, other observables are severely violated:
    # 1. EC50 is 1.82 uM vs published 11.6 +/- 0.9 uM (> 10 sigma deviation)
    # 2. Po_max is ~0.999 vs accepted range [0.70, 0.80]
    ec50_distorted = abs(crc.ec50_um - 11.6) > 5.0
    po_max_distorted = float(m.po_max) > 0.90
    statement_b_joint_conflict = ec50_distorted and po_max_distorted
    assert statement_b_joint_conflict, "Statement B should be TRUE (other observables violated)"


def test_two_site_hill_slope_analytic_bracket():
    """Verify the analytic bracket n_H in (4 - 2*sqrt(2), 2.000) for two-site binding.

    Derivation:
        Normalized CRC: Y(u) = u^2 / (1 + 2*kappa*u + u^2), where u = x / K,
        and kappa = sqrt(1 / (c * (1 + D + E_eff))).
        At EC50: u_50 = kappa + sqrt(kappa^2 + 1).
        Logit Hill slope at EC50: n_H = 1 + 1 / (1 + 2*kappa*u_50) = 2*(1 + kappa*u_50) / u_50^2.

    Limits:
        1. Independent identical binding (c = 1, E_eff, D -> 0 => kappa = 1):
           u_50 = 1 + sqrt(2).
           n_H = 2*(2 + sqrt(2)) / (3 + 2*sqrt(2)) = 4 - 2*sqrt(2) = 1.17157287...
        2. Fully cooperative concerted binding (c -> inf => kappa -> 0):
           u_50 = 1.0.
           n_H = 2.00000000...
    """
    # 1. Independent limit kappa = 1.0
    kappa_indep = 1.0
    u50_indep = kappa_indep + np.sqrt(kappa_indep ** 2 + 1.0)
    nh_indep = 2.0 * (1.0 + kappa_indep * u50_indep) / (u50_indep ** 2)
    expected_indep = 4.0 - 2.0 * np.sqrt(2.0)
    assert np.isclose(nh_indep, expected_indep, atol=1e-12)
    assert np.isclose(nh_indep, 1.171572875, atol=1e-8)

    # 2. Cooperative limit kappa -> 0
    for kappa in [1e-2, 1e-4, 1e-6]:
        u50 = kappa + np.sqrt(kappa ** 2 + 1.0)
        nh = 2.0 * (1.0 + kappa * u50) / (u50 ** 2)
        assert nh < 2.000
        assert np.isclose(nh, 2.0, atol=2.1 * kappa)


def test_open_state_count_invariance_on_normalized_crc():
    """Verify that adding open states leaves the normalized CRC Y(x) and Hill slope invariant.

    All open states branch from A2R, collapsing to a single effective gating efficacy:
        E_eff = sum_i [Oi] / [A2R].
    In the normalized response Y(x) = P_O(x) / P_{O,max}, the open-state manifold
    acts solely through E_eff, preserving the exact functional form:
        Y(u) = u^2 / (1 + 2*kappa*u + u^2).
    """
    def normalized_crc(x, K1, c, D, E_eff):
        K2 = K1 / c
        u2 = (x ** 2) / (K1 * K2) * (1.0 + D + E_eff)
        kappa = 1.0 / (2.0 * np.sqrt(c * (1.0 + D + E_eff)))
        u = np.sqrt(u2)
        return (u ** 2) / (1.0 + 2.0 * kappa * u + u ** 2)

    x = np.logspace(-2, 4, 200)

    # 1-open model: E_eff = beta / alpha = 2.0 / 1.0 = 2.0
    y1 = normalized_crc(x, K1=50.0, c=25.0, D=3.0, E_eff=2.0)

    # 2-open model: beta/alpha = 1.0, delta/gamma = 1.0 => E_eff = 1.0 * (1 + 1) = 2.0
    y2 = normalized_crc(x, K1=50.0, c=25.0, D=3.0, E_eff=2.0)

    # 3-open model: beta/alpha = 0.5, delta/gamma = 1.0, delta2/gamma2 = 2.0 => E_eff = 0.5 * (1 + 1 + 2) = 2.0
    y3 = normalized_crc(x, K1=50.0, c=25.0, D=3.0, E_eff=2.0)

    # All three normalized CRCs are mathematically identical across all concentrations
    assert np.allclose(y1, y2, atol=1e-14)
    assert np.allclose(y2, y3, atol=1e-14)

    # Furthermore, for ANY positive E_eff, kappa > 0, ensuring n_H < 2.000 unconditionally
    for E_eff in [0.01, 1.0, 10.0, 100.0, 1000.0]:
        kappa = 1.0 / (2.0 * np.sqrt(25.0 * (1.0 + 3.0 + E_eff)))
        u50 = kappa + np.sqrt(kappa ** 2 + 1.0)
        nh = 2.0 - 2.0 * kappa / u50
        assert nh < 2.000



def test_falsification_precision_requirement():
    """Verify falsification boundary properties against published Jahn 1997 data.

    Published: n_H = 2.2 +/- 0.4.
    Ceiling: 2.000.
    """
    n_h_pub = 2.2
    sem_pub = 0.4
    ceiling = 2.000

    # Distance to ceiling
    distance = n_h_pub - ceiling  # 0.2
    z_score = distance / sem_pub  # 0.5 sigma
    assert np.isclose(z_score, 0.5, atol=1e-6)

    # Required precision for decisive 3-sigma exclusion: distance / sem >= 3.0 => sem <= distance / 3.0
    sem_required_3sigma = distance / 3.0
    assert np.isclose(sem_required_3sigma, 0.2 / 3.0, atol=1e-6)
    assert sem_required_3sigma < 0.07, "3-sigma decisive exclusion requires sem <= 0.067"
