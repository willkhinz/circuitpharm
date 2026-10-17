"""Tests for CRC component extraction, Weighting iv likelihood scoring, and discrepancies.

AGENT BRIEF 12 (Task 1):
- CRC scored against EC50 and nH at published errors; values and errors read from
  project records, source named, no literals.
- FIT_TARGETS['ec50_um'] discrepancy surfaced with both numbers, unresolved.
- Hill slope extracted over the measured 1-10 uM window, with a full-range sensitivity check reported.
- Normalisation convention stated (fraction_of_max).
- Nine-point objective retained as a reported view.
"""
from __future__ import annotations

import numpy as np
import pytest

from circuitpharm.fitting import (
    CRC_BASELINE_RMSE,
    JAHN1997_PEAK_CRC,
    JAHN1997_PEAK_CRC_RESTRICTED,
    JAHN_EC50_SEM,
    JAHN_EC50_UM,
    JAHN_HILL,
    JAHN_HILL_SEM,
    JAHN_SLOPE_WINDOW_UM,
    CrcComponentsLogLikelihood,
    CrcParameters,
    crc_components_log_likelihood,
    extract_crc_parameters,
    get_jahn_published_records,
    surface_ec50_discrepancy,
)
from circuitpharm.gabaa_kinetics import FIT_TARGETS
from circuitpharm.models.base import (
    JAHN1997_SLOPE_WINDOW_REL,
    hill_slope,
)
from circuitpharm.models.kinetic_jw95 import KineticAllosteryModel
from circuitpharm.models.kinetic_two_open import KineticTwoOpenModel


def test_jahn_published_records_read_from_project_records():
    """Verify published values and errors are read directly from project records without literals."""
    rec = get_jahn_published_records()
    assert rec["ec50_um"] == pytest.approx(11.6)
    assert rec["ec50_sem"] == pytest.approx(0.9)
    assert rec["hill"] == pytest.approx(2.2)
    assert rec["hill_sem"] == pytest.approx(0.4)
    assert rec["measured_window_um"] == (1.0, 10.0)
    assert "Jahn et al. 1997" in rec["source"]
    assert rec["dataset_key"] == JAHN1997_PEAK_CRC.citation_label
    assert rec["restricted_dataset_key"] == JAHN1997_PEAK_CRC_RESTRICTED.citation_label


def test_ec50_discrepancy_surfaced_and_resolved():
    """Verify the discrepancy between FIT_TARGETS['ec50_um'] and the sourced value is surfaced and resolved."""
    disc = surface_ec50_discrepancy()
    assert disc["fit_targets_ec50_um"] == pytest.approx(11.6)
    assert disc["jahn_sourced_ec50_um"] == pytest.approx(11.6)
    assert disc["jahn_sourced_ec50_sem"] == pytest.approx(0.9)
    assert disc["discrepancy_um"] == pytest.approx(0.0)
    assert disc["discrepancy_sigma"] == pytest.approx(0.0)
    assert "resolved" in disc["resolution"].lower()


def test_hill_slope_measured_window_and_sensitivity():
    """Verify hill_slope extracts over the measured window (1-10 uM) and full range."""
    # Test on a pure Hill curve with published parameters
    c = np.logspace(-2, 4.5, 400)
    x = (c / 11.6) ** 2.2
    y = x / (1.0 + x)

    # Over measured window
    sl_meas = hill_slope(c, y, window_um=(1.0, 10.0))
    assert sl_meas == pytest.approx(2.2, rel=1e-3)

    # Over full curve band
    sl_full = hill_slope(c, y, band=(0.02, 0.98))
    assert sl_full == pytest.approx(2.2, rel=1e-3)

    # Over relative window
    sl_rel = hill_slope(c, y, window_rel=JAHN1997_SLOPE_WINDOW_REL)
    assert sl_rel == pytest.approx(2.2, rel=1e-3)


def test_extract_crc_parameters_normalisation_and_extraction():
    """Verify normalisation convention (fraction_of_max) and parameter extraction on models."""
    m = KineticTwoOpenModel()
    params = extract_crc_parameters(m)

    assert isinstance(params, CrcParameters)
    assert params.normalisation == "fraction_of_max"
    assert np.isfinite(params.ec50_um)
    assert params.ec50_um > 0
    assert np.isfinite(params.hill_measured)
    assert np.isfinite(params.hill_full)
    assert np.isfinite(params.hill_rel)
    assert 0.0 < params.y_max <= 1.0

    # Bending receptor curve: local measured slope differs from full-range regression
    assert abs(params.hill_measured - params.hill_full) > 0.1


def test_crc_components_log_likelihood_scoring():
    """Verify CRC log-likelihood and residuals in sigma."""
    m = KineticTwoOpenModel()
    ll = crc_components_log_likelihood(m, compute_9point_view=True)

    assert isinstance(ll, CrcComponentsLogLikelihood)
    assert ll.z_ec50 == pytest.approx((ll.ec50_um - JAHN_EC50_UM) / JAHN_EC50_SEM)
    assert ll.z_hill == pytest.approx((ll.hill_slope - JAHN_HILL) / JAHN_HILL_SEM)
    assert ll.chi2_ec50 == pytest.approx(ll.z_ec50 ** 2)
    assert ll.chi2_hill == pytest.approx(ll.z_hill ** 2)
    assert ll.chi2_total == pytest.approx(ll.chi2_ec50 + ll.chi2_hill)
    assert ll.max_z == pytest.approx(max(abs(ll.z_ec50), abs(ll.z_hill)))

    # 9-point reported view retained
    assert np.isfinite(ll.crc_rmse)
    assert np.isfinite(ll.crc_ratio)
    assert ll.crc_ratio == pytest.approx(ll.crc_rmse / CRC_BASELINE_RMSE)

    # Sensitivity checks
    assert np.isfinite(ll.hill_slope_full)
    assert np.isfinite(ll.z_hill_full)
    assert np.isfinite(ll.hill_slope_rel)
    assert np.isfinite(ll.z_hill_rel)


def test_crc_synthetic_recovery():
    """Synthetic check: perfect model matching published targets yields chi2 ~ 0."""
    # Synthetic parameters perfectly at published values
    p_exact = CrcParameters(
        ec50_um=float(JAHN_EC50_UM),
        hill_measured=float(JAHN_HILL),
        hill_rel=float(JAHN_HILL),
        hill_full=float(JAHN_HILL),
        y_max=0.75,
    )
    ll_exact = crc_components_log_likelihood(p_exact, compute_9point_view=False)
    assert ll_exact.z_ec50 == pytest.approx(0.0, abs=1e-8)
    assert ll_exact.z_hill == pytest.approx(0.0, abs=1e-8)
    assert ll_exact.chi2_total == pytest.approx(0.0, abs=1e-8)
    assert ll_exact.max_z == pytest.approx(0.0, abs=1e-8)
