"""Tests for Task B: isolating dataset pull against holdout."""
import pytest
import numpy as np

from circuitpharm.fitting.data import (
    JAHN1997_PEAK_CRC,
    BARBERIS2007_DEACT,
    BARBERIS2007_PAIRED_PULSE,
)
from circuitpharm.fitting.comparison import MODEL_SPECS, fit_one
from circuitpharm.fitting.likelihood import predict_paired_pulse_recovery


@pytest.mark.slow
def test_deactivation_flips_holdout_recovery_error_sign():
    """Verify that CRC alone under-recovers while deactivation drives severe over-recovery."""
    spec = MODEL_SPECS["kinetic_jw95"]

    # 1. CRC alone
    fit_crc = fit_one(spec, [JAHN1997_PEAK_CRC], n_starts=4, seed=42, anchor_po_max=True)
    m_crc = spec.factory(np.array([fit_crc.values[k] for k in spec.param_names]))
    r_crc = float(predict_paired_pulse_recovery(m_crc, BARBERIS2007_PAIRED_PULSE))
    resid_crc = (r_crc - float(BARBERIS2007_PAIRED_PULSE.mean_response)) / float(BARBERIS2007_PAIRED_PULSE.sem)
    assert resid_crc < -3.0, f"CRC alone should under-recover by >3 sigma, got {resid_crc:.2f} sigma"

    # 2. Deact alone
    fit_deact = fit_one(spec, [BARBERIS2007_DEACT], n_starts=4, seed=42, anchor_po_max=True)
    m_deact = spec.factory(np.array([fit_deact.values[k] for k in spec.param_names]))
    r_deact = float(predict_paired_pulse_recovery(m_deact, BARBERIS2007_PAIRED_PULSE))
    resid_deact = (r_deact - float(BARBERIS2007_PAIRED_PULSE.mean_response)) / float(BARBERIS2007_PAIRED_PULSE.sem)
    assert resid_deact > 10.0, f"Deact alone should severely over-recover by >10 sigma, got {resid_deact:.2f} sigma"

    # 3. Both
    fit_both = fit_one(spec, [JAHN1997_PEAK_CRC, BARBERIS2007_DEACT], n_starts=4, seed=42, anchor_po_max=True)
    m_both = spec.factory(np.array([fit_both.values[k] for k in spec.param_names]))
    r_both = float(predict_paired_pulse_recovery(m_both, BARBERIS2007_PAIRED_PULSE))
    resid_both = (r_both - float(BARBERIS2007_PAIRED_PULSE.mean_response)) / float(BARBERIS2007_PAIRED_PULSE.sem)
    assert resid_both > 10.0, f"Both should over-recover by >10 sigma, got {resid_both:.2f} sigma"

    # Assert that deactivation flips the sign from negative to positive
    assert resid_crc < 0 < resid_both


def test_known_optima_isolation_invariants():
    """Fast check: evaluate holdout predictions at the three known optima to verify sign flip."""
    spec = MODEL_SPECS["kinetic_jw95"]
    pp_mean = float(BARBERIS2007_PAIRED_PULSE.mean_response)  # 0.33
    pp_sem = float(BARBERIS2007_PAIRED_PULSE.sem)             # 0.03

    # 1. CRC alone optimum
    m_crc = spec.factory(np.array([2.0183, 2.2127, -3.0000]))
    r_crc = float(predict_paired_pulse_recovery(m_crc, BARBERIS2007_PAIRED_PULSE))
    sig_crc = (r_crc - pp_mean) / pp_sem
    assert -5.0 < sig_crc < -4.0, f"CRC alone error {sig_crc:.2f}σ should be ~ -4.5σ"

    # 2. Deact alone optimum
    m_deact = spec.factory(np.array([3.0579, 0.1606, 1.4117]))
    r_deact = float(predict_paired_pulse_recovery(m_deact, BARBERIS2007_PAIRED_PULSE))
    sig_deact = (r_deact - pp_mean) / pp_sem
    assert sig_deact > 15.0, f"Deact alone error {sig_deact:.2f}σ should be > +15σ"

    # 3. Both (CRC + Deactivation)
    m_both = spec.factory(np.array([0.9474, 0.1799, 1.4060]))
    r_both = float(predict_paired_pulse_recovery(m_both, BARBERIS2007_PAIRED_PULSE))
    sig_both = (r_both - pp_mean) / pp_sem
    assert sig_both > 15.0, f"Both error {sig_both:.2f}σ should be > +15σ"

    # Sign flip verified
    assert sig_crc < 0.0 < sig_deact and sig_crc < 0.0 < sig_both
