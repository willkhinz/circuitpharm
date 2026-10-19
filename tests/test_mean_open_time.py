"""Tests for generator-derived mean open time predictor and datasets (Brief 15 Task 1).

Covers:
  1. Generator matrix Q derivation equivalence for multi-open-state schemes.
  2. Single-open-state schemes (JW95, OpenDesens, Flip) yielding 1/alpha (or 1/(alpha+d)).
  3. Two-open-state schemes (TwoOpen, Coop) matching pi_O @ (-Q_OO)^-1 @ 1 exactly.
  4. Operational model refusal (ObservableMismatch).
  5. Dataset declarations and band definitions for Keramidas 2008 and Li 2008.
"""
import numpy as np
import pytest

from circuitpharm.fitting.data import (
    KERAMIDAS2008_MEAN_OPEN_TIME_H,
    KERAMIDAS2008_MEAN_OPEN_TIME_M,
    KERAMIDAS_H_MODE_BAND,
    KERAMIDAS_M_MODE_BAND,
    KERAMIDAS_MODAL_SPAN,
    LI2008_MEAN_OPEN_TIME,
    MeanOpenTimeDataset,
)
from circuitpharm.fitting.likelihood import (
    predict_mean_open_time,
)
from circuitpharm.models import (
    KineticAllosteryModel,
    OperationalScalarModel,
)
from circuitpharm.models.base import Observable, ObservableMismatch
from circuitpharm.models.kinetic_coop import KineticCoopModel
from circuitpharm.models.kinetic_flip import KineticFlipModel
from circuitpharm.models.kinetic_open_desens import KineticOpenDesensModel
from circuitpharm.models.kinetic_two_open import KineticTwoOpenModel


def test_single_open_state_jw95():
    """For 1-open-state JW95, mean open time is exactly 1/alpha."""
    model = KineticAllosteryModel(alpha=0.15)
    tau = predict_mean_open_time(model)
    assert np.isclose(tau, 1.0 / 0.15)


def test_single_open_state_open_desens():
    """For KineticOpenDesens, exiting open state can be closing (alpha) or desensitising (d),
    so mean open time is 1/(alpha + d)."""
    model = KineticOpenDesensModel(alpha=0.20, d=0.05)
    tau = predict_mean_open_time(model)
    assert np.isclose(tau, 1.0 / (0.20 + 0.05))


def test_single_open_state_flip():
    """For KineticFlipModel, intermediate shut flip state does not change open dwell time 1/alpha."""
    model = KineticFlipModel(alpha=0.25)
    tau = predict_mean_open_time(model)
    assert np.isclose(tau, 1.0 / 0.25)


def test_two_open_state_generator_matrix_equivalence():
    """Verify analytical formula (1/alpha)*(1 + delta/gamma) matches fundamental matrix
    inversion pi_O @ inv(-Q_OO) @ 1 to machine precision."""
    alphas = [0.05, 0.1, 0.5, 1.0]
    deltas = [0.01, 0.2, 0.8, 2.0]
    gammas = [0.05, 0.3, 1.5, 5.0]

    for a in alphas:
        for d in deltas:
            for g in gammas:
                model = KineticTwoOpenModel(alpha=a, delta=d, gamma=g)
                tau_analytic = float(predict_mean_open_time(model))
                expected = (1.0 / a) * (1.0 + d / g)
                assert np.isclose(tau_analytic, expected, rtol=1e-14)

                # Direct generator matrix sub-block inversion:
                # States: R, AR, A2R, A2O (idx 3), A2D, A2O' (idx 5)
                # Open set O = {A2O, A2O'} -> indices 3, 5
                Q = model.q_matrix(0.0)
                open_indices = [3, 4]
                Q_OO = Q[np.ix_(open_indices, open_indices)]
                # Q_OO must be [[-(a + d), d], [g, -g]]
                assert np.isclose(Q_OO[0, 0], -(a + d))
                assert np.isclose(Q_OO[0, 1], d)
                assert np.isclose(Q_OO[1, 0], g)
                assert np.isclose(Q_OO[1, 1], -g)

                inv_neg_Q = np.linalg.inv(-Q_OO)
                # Entry flux is strictly into A2O from A2R, so pi_O = [1, 0]
                pi_O = np.array([1.0, 0.0])
                tau_matrix = float(pi_O @ inv_neg_Q @ np.ones(2))
                assert np.isclose(tau_analytic, tau_matrix, rtol=1e-12)


def test_kinetic_coop_matches_two_open_formula():
    """KineticCoopModel inherits the same 2-open-state topology (A2O <-> A2O') and must evaluate identically."""
    model = KineticCoopModel(alpha=0.12, delta=0.45, gamma=0.30)
    tau = float(predict_mean_open_time(model))
    expected = (1.0 / 0.12) * (1.0 + 0.45 / 0.30)
    assert np.isclose(tau, expected, rtol=1e-14)


def test_two_open_limits_to_one_open():
    """As delta -> 0 (no transitions to secondary open state) or gamma -> inf (infinitely brief),
    mean open time continuously limits to 1/alpha."""
    m_zero_delta = KineticTwoOpenModel(alpha=0.1, delta=0.0, gamma=1.0)
    assert np.isclose(predict_mean_open_time(m_zero_delta), 10.0)

    m_fast_gamma = KineticTwoOpenModel(alpha=0.1, delta=1.0, gamma=1e8)
    assert np.isclose(predict_mean_open_time(m_fast_gamma), 10.0, rtol=1e-6)


def test_operational_model_refused():
    """Operational models have no Markov states and must raise ObservableMismatch."""
    op_model = OperationalScalarModel()
    with pytest.raises(ObservableMismatch, match="no Markov states or closing rate alpha"):
        predict_mean_open_time(op_model)


def test_dataset_declarations():
    """Verify dataset fields, bands, and provenance metadata."""
    assert KERAMIDAS2008_MEAN_OPEN_TIME_M.observable is Observable.MEAN_OPEN_TIME
    assert KERAMIDAS2008_MEAN_OPEN_TIME_M.mean_open_time_ms == 2.68
    assert KERAMIDAS2008_MEAN_OPEN_TIME_M.sem_ms is None
    assert KERAMIDAS2008_MEAN_OPEN_TIME_M.band == KERAMIDAS_M_MODE_BAND
    assert KERAMIDAS_M_MODE_BAND == (2.54, 2.83)

    assert KERAMIDAS2008_MEAN_OPEN_TIME_H.observable is Observable.MEAN_OPEN_TIME
    assert KERAMIDAS2008_MEAN_OPEN_TIME_H.mean_open_time_ms == 7.25
    assert KERAMIDAS2008_MEAN_OPEN_TIME_H.sem_ms is None
    assert KERAMIDAS2008_MEAN_OPEN_TIME_H.band == KERAMIDAS_H_MODE_BAND
    assert KERAMIDAS_H_MODE_BAND == (5.86, 8.63)

    assert KERAMIDAS_MODAL_SPAN == (2.54, 8.63)

    assert LI2008_MEAN_OPEN_TIME.observable is Observable.MEAN_OPEN_TIME
    assert LI2008_MEAN_OPEN_TIME.mean_open_time_ms == 2.96
    assert LI2008_MEAN_OPEN_TIME.sem_ms is None
    assert LI2008_MEAN_OPEN_TIME.n_cells == 4
    assert "gamma2L" in LI2008_MEAN_OPEN_TIME.preparation
