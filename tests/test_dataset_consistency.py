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
