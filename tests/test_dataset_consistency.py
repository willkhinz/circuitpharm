"""Tests for the dataset consistency diagnostic (Brief 4 Task A)."""
import pytest
import numpy as np

from circuitpharm.fitting.data import (
    JAHN1997_PEAK_CRC,
    BARBERIS2007_DEACT,
    BARBERIS2007_PAIRED_PULSE,
    TRAIN,
    HOLDOUT,
    holdout_guard,
)
from circuitpharm.fitting.comparison import MODEL_SPECS
from circuitpharm.fitting.likelihood import (
    make_log_likelihood,
    _predict,
    predict_paired_pulse_recovery,
)
from circuitpharm.fitting.diagnostic import (
    diagnostic_log_likelihood,
    fit_diagnostic_consistency,
    DiagnosticFitResult,
)


def test_holdout_guard_is_not_weakened_or_bypassed():
    """Verify that holdout_guard strictly refuses held-out data on the standard fitting path."""
    with pytest.raises(ValueError, match="held-out dataset"):
        make_log_likelihood([BARBERIS2007_PAIRED_PULSE], caller="test_fitting_path", n_params=3)

    with pytest.raises(ValueError, match="held-out dataset"):
        make_log_likelihood([JAHN1997_PEAK_CRC, BARBERIS2007_PAIRED_PULSE], caller="test_fitting_path", n_params=3)

    # Train datasets pass holdout_guard without error
    holdout_guard(list(TRAIN), caller="test_train_guard")


def test_diagnostic_objective_computes_joint_likelihood():
    """Verify diagnostic_log_likelihood computes scores with declared noise models and anchor."""
    spec = MODEL_SPECS["kinetic_jw95"]
    v = np.array([3.0, 0.5, 2.5])
    model = spec.factory(v)

    dsets = [JAHN1997_PEAK_CRC, BARBERIS2007_DEACT, BARBERIS2007_PAIRED_PULSE]
    res = diagnostic_log_likelihood(model, dsets, n_params=3, anchor_po_max=True)

    # 3 datasets + 1 anchor = 4 entries in per_dataset
    assert len(res.per_dataset) == 4
    assert np.isfinite(res.total_log_likelihood)

    # Check that noise models are respected
    labels = {f.label: f for f in res.per_dataset}
    assert labels[JAHN1997_PEAK_CRC.citation_label].noise_model == "estimated"
    assert labels[BARBERIS2007_DEACT.citation_label].noise_model == "estimated"
    assert labels[BARBERIS2007_PAIRED_PULSE.citation_label].noise_model == "stated_sem"
    assert labels["po_max gating anchor"].noise_model == "stated_sem"


def test_diagnostic_fit_reports_invariants_and_verdict():
    """Verify fit_diagnostic_consistency returns a valid DiagnosticFitResult."""
    spec = MODEL_SPECS["kinetic_jw95"]
    # Fast test run with 1 start
    fit_res = fit_diagnostic_consistency(
        spec.factory,
        [JAHN1997_PEAK_CRC],
        param_names=spec.param_names,
        bounds=spec.bounds,
        starts=[np.array([2.02, 2.21, -3.00])],
        anchor_po_max=True,
    )
    assert isinstance(fit_res, DiagnosticFitResult)
    assert "log10_D" in fit_res.values
    assert np.isfinite(fit_res.total_log_likelihood)
    assert fit_res.po_max is not None
    assert "log10_D" in fit_res.at_bound  # CRC alone drives log10_D to lower bound


def test_datasets_are_not_mutually_consistent():
    """Verify that no tested parameter set satisfies all three real datasets within 3 sigma."""
    spec = MODEL_SPECS["kinetic_jw95"]
    # At the joint optimum or holdout-satisfying point:
    v_joint = np.array([0.8632, 0.0984, 1.9759])
    m_joint = spec.factory(v_joint)
    
    # Check holdout recovery residual
    pp_pred = float(predict_paired_pulse_recovery(m_joint, BARBERIS2007_PAIRED_PULSE))
    pp_resid_sigma = abs(pp_pred - float(BARBERIS2007_PAIRED_PULSE.mean_response)) / float(BARBERIS2007_PAIRED_PULSE.sem)
    assert pp_resid_sigma > 3.0, f"Recovery residual {pp_resid_sigma:.2f} should exceed 3 sigma"

    # And at the point satisfying recovery (e.g. 3.0, 0.5, 2.5), deactivation residual is huge
    m_holdout = spec.factory(np.array([3.0, 0.5, 2.5]))
    deact_pred = _predict(m_holdout, BARBERIS2007_DEACT)
    deact_rmse = np.sqrt(np.mean((deact_pred - BARBERIS2007_DEACT.normalized_current) ** 2))
    assert deact_rmse > 0.25, f"Deactivation RMSE {deact_rmse:.4f} should be large when recovery is satisfied"


def test_enlarged_specs_exist_and_construct_models():
    """Verify that kinetic_jw95_free_r and kinetic_jw95_free_rates construct valid models."""
    assert "kinetic_jw95_free_r" in MODEL_SPECS
    assert "kinetic_jw95_free_rates" in MODEL_SPECS

    spec_4 = MODEL_SPECS["kinetic_jw95_free_r"]
    assert len(spec_4.param_names) == 4
    m4 = spec_4.factory(np.array([2.0, 0.5, 1.0, -2.7]))
    assert hasattr(m4, "q_matrix")
    assert 0.0 < m4.po_max <= 1.0

    spec_6 = MODEL_SPECS["kinetic_jw95_free_rates"]
    assert len(spec_6.param_names) == 6
    m6 = spec_6.factory(np.array([2.0, 0.5, 1.0, -2.7, -0.74, -0.44]))
    assert hasattr(m6, "q_matrix")
    assert 0.0 < m6.po_max <= 1.0


def test_conflict_survives_freeing_rate_conventions():
    """Verify that freeing rate conventions (r, koff, alpha) does not reconcile Barberis observables."""
    spec = MODEL_SPECS["kinetic_jw95_free_rates"]

    # Point from 6-param Deact alone fit:
    v_deact = np.array([3.06, 0.16, 1.41, -2.70, -0.74, -0.44])
    m_deact = spec.factory(v_deact)
    pp_pred_deact = float(predict_paired_pulse_recovery(m_deact, BARBERIS2007_PAIRED_PULSE))
    sig_deact = abs(pp_pred_deact - 0.33) / 0.03
    assert sig_deact > 10.0, f"Deactivation alone miss {sig_deact:.2f} sigma should be large"

    # Joint 6-parameter fit result from Task B:
    # {'log10_kd': 3.99, 'log10_E': 0.47, 'log10_D': 2.72, 'log10_r': -1.97, 'log10_koff': 0.73, 'log10_alpha': -0.47}
    v_joint = np.array([3.9902, 0.4749, 2.7198, -1.9659, 0.7263, -0.4747])
    m_joint = spec.factory(v_joint)
    pp_pred_joint = float(predict_paired_pulse_recovery(m_joint, BARBERIS2007_PAIRED_PULSE))
    sig_joint = abs(pp_pred_joint - 0.33) / 0.03
    assert sig_joint > 3.0, f"Joint recovery miss {sig_joint:.2f} sigma must exceed 3 sigma"

    # And CRC RMSE is blown out (>50x baseline 0.004917)
    crc_pred = _predict(m_joint, JAHN1997_PEAK_CRC)
    crc_rmse = np.sqrt(np.mean((crc_pred - JAHN1997_PEAK_CRC.mean_response) ** 2))
    assert crc_rmse / 0.004917 > 50.0, f"CRC RMSE ratio {crc_rmse/0.004917:.1f}x should be blown out"

