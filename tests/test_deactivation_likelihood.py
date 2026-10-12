"""The deactivation half of the likelihood, which did not exist until 2026-10-11.

WHY THESE TESTS EXIST. `fitting.data.MISSING_DATASETS` calls a macroscopic deactivation time
course "the dataset that makes an absolute rate (k_off) identifiable at all", and
`DEACTIVATION_BENCHMARK`'s own `origin` says it "exists so the kinetic-timescale half of the
likelihood can be exercised before P1-2". Neither was true: `_predict` raised
`ObservableMismatch` for every observable except PEAK and EQUILIBRIUM, `_observed` read only
`mean_response`, and NOTHING in `src/` read a deactivation trace's `time_ms` or
`normalized_current` at all. A deactivation dataset was inert however good it was -- which is
why adding a correct, attributable Barberis trace could not, on its own, change what the data
determines.

These tests pin the half that was added, and in particular the three things that could go
wrong silently:

  * the PROTOCOL must come from the dataset, never from a default -- Barberis measured the
    same receptor at 52.5 ms (2 ms pulse) and 364 ms (3 s pulse), 6.9x apart, so a default
    pulse length would silently choose which measurement a model is compared against;
  * the integration grid must not move the score, because it is a numerical setting and not
    data (a comment in `predict_deactivation` claims this, so it is asserted here);
  * a model with no states to relax between must be REFUSED rather than scored, because
    every model in this project has a `simulate_waveform` and `OperationalScalarModel`'s
    relaxes towards a Hill target at an assumed time constant -- so it would return a
    "deactivation" that measures that assumption and not a rate.
"""
import warnings

import numpy as np
import pytest

from circuitpharm import (
    ExtendedDesensitizationModel, KineticAllosteryModel, OperationalScalarModel,
)
from circuitpharm.fitting.data import ALL, DeactivationDataset
from circuitpharm.fitting.likelihood import (
    _observed, make_log_likelihood, predict_deactivation,
)
from circuitpharm.models.base import Observable, ObservableMismatch

BARBERIS = ALL["barberis2007_deact"]


def _score(datasets, model, n_params=4):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return make_log_likelihood(datasets, caller="test", n_params=n_params)(model)


# ---------------------------------------------------------------------------------------
# The prediction itself.
# ---------------------------------------------------------------------------------------

def test_a_deactivation_dataset_can_now_be_scored_at_all():
    """The regression this file exists for. Before `predict_deactivation`, this raised
    AttributeError: 'DeactivationDataset' object has no attribute 'mean_response'."""
    r = _score([BARBERIS], KineticAllosteryModel())
    assert np.isfinite(r.total_log_likelihood)
    assert r.per_dataset[0].n == len(BARBERIS.time_ms)


def test_the_prediction_is_a_decay_normalised_to_its_own_peak():
    pred = predict_deactivation(KineticAllosteryModel(), BARBERIS)
    assert pred.shape == np.asarray(BARBERIS.time_ms).shape
    assert pred[0] == pytest.approx(1.0, abs=1e-9), "a decay must start at its own peak"
    assert np.all(pred <= 1.0 + 1e-9)
    assert pred[-1] < pred[0], "it must actually decay"
    assert np.all(np.isfinite(pred))


def test_observed_reads_the_decay_column_and_dispatches_on_the_observable():
    obs = _observed(BARBERIS)
    assert np.allclose(obs, np.asarray(BARBERIS.normalized_current, float))
    # and a concentration-response dataset still reads its own column
    crc = ALL["jahn1997_peak_crc"]
    assert np.allclose(_observed(crc), np.asarray(crc.mean_response, float))


# ---------------------------------------------------------------------------------------
# The grid is a numerical setting, not data.
# ---------------------------------------------------------------------------------------

@pytest.mark.parametrize("n_pulse,n_tail", [(200, 800), (400, 1600), (800, 3200)])
def test_the_integration_grid_does_not_move_the_prediction(n_pulse, n_tail):
    """`predict_deactivation`'s comment claims the integration settings cannot move the
    score because the result is interpolated onto the dataset's own grid. Asserted, because
    an unasserted claim about a computation is exactly what this project distrusts."""
    ref = predict_deactivation(KineticAllosteryModel(), BARBERIS)
    got = predict_deactivation(KineticAllosteryModel(), BARBERIS,
                               n_pulse=n_pulse, n_tail=n_tail)
    assert np.max(np.abs(ref - got)) < 2e-3, (
        "the integration grid changed the prediction by more than interpolation error, so "
        "it is acting as data rather than as a numerical setting")


# ---------------------------------------------------------------------------------------
# The protocol comes from the dataset. This is the PEAK lesson, applied to deactivation.
# ---------------------------------------------------------------------------------------

def test_the_declared_pulse_duration_changes_the_prediction():
    """Barberis measured tau_w = 52.5 ms after a 2 ms pulse and 364 ms after a 3 s pulse --
    6.9x from the protocol alone. So the pulse length had better matter here, or the
    prediction is not reading the protocol it claims to."""
    model = KineticAllosteryModel()
    short = predict_deactivation(model, BARBERIS)
    long_ds = DeactivationDataset(
        citation_label="same trace, 3 s pulse declared",
        time_ms=BARBERIS.time_ms, normalized_current=BARBERIS.normalized_current,
        pulse_ms=3000.0, conc_mm=BARBERIS.conc_mm, n_components=3,
        fit_window_ms=BARBERIS.fit_window_ms, observable=Observable.DEACTIVATION,
        normalisation="fraction_of_max", kind="synthetic", synthetic=True,
        origin="the Barberis trace with a 3 s pulse declared, to prove the pulse is read",
    )
    long = predict_deactivation(model, long_ds)
    # the integrated decay is the weighted tau; a 3 s pulse must give a different one
    t = np.asarray(BARBERIS.time_ms, float)
    tw_short = float(np.trapezoid(short, t))
    tw_long = float(np.trapezoid(long, t))
    assert abs(tw_long - tw_short) > 1.0, (
        f"declaring a 3 s pulse instead of 2 ms moved the predicted weighted tau only from "
        f"{tw_short:.2f} to {tw_long:.2f} ms, so the declared protocol is not being used")


def test_a_deactivation_dataset_without_its_protocol_cannot_be_predicted():
    class _Bare:
        citation_label = "no protocol"
        observable = Observable.DEACTIVATION
        time_ms = np.linspace(0.0, 100.0, 10)
        normalized_current = np.linspace(1.0, 0.0, 10)
        pulse_ms = None
        conc_mm = 10.0
        fit_window_ms = 100.0
    with pytest.raises(ValueError, match="pulse_ms"):
        predict_deactivation(KineticAllosteryModel(), _Bare())


# ---------------------------------------------------------------------------------------
# Which models may be scored on a decay shape.
# ---------------------------------------------------------------------------------------

def test_a_model_with_no_states_is_refused_not_scored():
    """OperationalScalarModel HAS a simulate_waveform -- it relaxes towards an instantaneous
    Hill target -- so a `hasattr(model, 'simulate_waveform')` guard would let it through and
    score an assumed relaxation time constant as if it were a rate. The guard is on
    `q_matrix` for that reason, and this test is why."""
    m = OperationalScalarModel()
    assert hasattr(m, "simulate_waveform"), (
        "if this is ever false the guard may be weakened to simulate_waveform; until then "
        "it must stay on q_matrix")
    assert not hasattr(m, "q_matrix")
    with pytest.raises(ObservableMismatch, match="q_matrix"):
        predict_deactivation(m, BARBERIS)
    with pytest.raises(ObservableMismatch):
        _score([BARBERIS], m, n_params=3)


@pytest.mark.parametrize("cls", [KineticAllosteryModel, ExtendedDesensitizationModel])
def test_both_state_based_schemes_can_be_scored_on_a_decay(cls):
    r = _score([BARBERIS], cls())
    assert np.isfinite(r.total_log_likelihood)


def test_a_charge_dataset_is_still_refused_and_says_why():
    """CHARGE is an integral over a window, not a decay shape, and still has no door. The
    message must not now claim the likelihood handles it."""
    with pytest.raises(ObservableMismatch) as e:
        _score([ALL["deactivation_charge"]], KineticAllosteryModel())
    msg = str(e.value)
    assert "CHARGE" in msg and "DEACTIVATION" in msg


def test_the_barberis_dataset_is_the_one_being_scored():
    """Guards the dataset's identity, so a later edit cannot repoint this suite at a
    synthetic trace and still pass."""
    assert BARBERIS.kind == "parametric" and BARBERIS.synthetic is True
    assert BARBERIS.source_key == "barberis2007"
    assert BARBERIS.observable is Observable.DEACTIVATION
    assert BARBERIS.pulse_ms == 2.0 and BARBERIS.conc_mm == 10.0
    assert BARBERIS.n_components == 3
