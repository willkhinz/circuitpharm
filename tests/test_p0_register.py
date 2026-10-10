"""Regression tests for the P0 defect register (knowledge/09-strategic-roadmap.md §3).

Each test names its register item. Several are ANCHOR tests in the roadmap's §2.7 sense --
they assert a VALUE against something the code was not fitted to, with the tolerance
justified in the docstring -- and those are marked. A test here that only asserted
`isfinite` or a shape would not be evidence that the defect is fixed.

TESTS MARKED "PERMANENT RECORD" MUST BE AMENDED, NOT DELETED, when a later phase changes
the behaviour they pin. They are the only durable statement of why that phase exists.
"""
from __future__ import annotations

import warnings

import numpy as np
import pytest

from circuitpharm.models import (
    ExtendedDesensitizationModel,
    KineticAllosteryModel,
    OperationalScalarModel,
)
from circuitpharm.models.base import (
    DEFAULT_E_CL_MV,
    DecayFit,
    Observable,
    ObservableMismatch,
    ObservedQuantity,
    fit_biexponential_decay,
    resolve_initial_state,
)
from circuitpharm.protocols.waveforms import (
    synaptic_transient,
    synaptic_transient_biexp_clearance,
)

MARKOV_MODELS = [
    pytest.param(KineticAllosteryModel(), 5, id="kinetic_jw95"),
    pytest.param(ExtendedDesensitizationModel(), 6, id="extended_desens"),
]


def _grid(n=1200, span=60.0):
    return np.linspace(0.0, span, n)


# ======================================================================== P0-1
# A scalar initial state crashed both Markov models, though models/base.py documents
# "state vector OR initial open probability". np.asarray(0.1) is 0-d, so len() raised
# before the length check; the float form raised in len() directly.

@pytest.mark.parametrize("model,n_states", MARKOV_MODELS)
def test_p0_1_scalar_initial_state_is_honoured(model, n_states):
    """ANCHOR. A scalar must be USED, not replaced by the resting distribution.

    Tolerance 1e-9 on the first sample: the solver is asked to start exactly there, so
    anything looser would hide the previous behaviour, which silently substituted the
    resting state and discarded the caller's instruction.
    """
    t = _grid(400, 10.0)
    g = np.full_like(t, 0.4)
    res = model.simulate_waveform(t, g, initial_state=0.25)
    assert res.p_open[0] == pytest.approx(0.25, abs=1e-9)
    assert res.states is not None
    assert res.states[0].sum() == pytest.approx(1.0, abs=1e-9)


@pytest.mark.parametrize("model,n_states", MARKOV_MODELS)
def test_p0_1_full_vector_accepted_and_validated(model, n_states):
    rest = model.state_distribution(0.0)
    t = _grid(200, 5.0)
    g = np.full_like(t, 0.4)
    res = model.simulate_waveform(t, g, initial_state=rest)
    assert res.p_open[0] == pytest.approx(float(rest[3]), abs=1e-9)

    with pytest.raises(ValueError, match="length"):
        model.simulate_waveform(t, g, initial_state=np.ones(n_states - 1) / (n_states - 1))
    bad = np.zeros(n_states)
    bad[0] = 1.3
    with pytest.raises(ValueError, match="sums to"):
        model.simulate_waveform(t, g, initial_state=bad)
    neg = np.zeros(n_states)
    neg[0], neg[1] = 1.5, -0.5
    with pytest.raises(ValueError, match="negative"):
        model.simulate_waveform(t, g, initial_state=neg)


def test_p0_1_scalar_out_of_range_raises():
    m = KineticAllosteryModel()
    t = _grid(100, 5.0)
    g = np.full_like(t, 0.4)
    for bad in (-0.1, 1.5):
        with pytest.raises(ValueError, match=r"\[0, 1\]"):
            m.simulate_waveform(t, g, initial_state=bad)


def test_p0_1_resolver_preserves_total_probability():
    """The scalar branch redistributes the remaining mass; it must still sum to 1."""
    rest = KineticAllosteryModel().state_distribution(0.0)
    for p in (0.0, 0.01, 0.5, 1.0):
        out = resolve_initial_state(p, 5, rest)
        assert out.sum() == pytest.approx(1.0, abs=1e-12)
        assert out[3] == pytest.approx(p, abs=1e-12)
        assert np.all(out >= 0.0)


def test_p0_1_operational_rejects_a_state_vector():
    """The mirror defect: Model A read element 0 of a state vector as an open probability.

    For a Markov resting distribution element 0 is P_R ~ 1.0, so a resting state started
    the operational trace fully open.
    """
    a = OperationalScalarModel()
    t = _grid(200, 20.0)
    g = synaptic_transient(t, peak_um=1000.0, rise_ms=0.1, clear_ms=1.0)
    rest = KineticAllosteryModel().state_distribution(0.0)
    with pytest.raises(ValueError, match="ONE state variable"):
        a.simulate_waveform(t, g, initial_state=rest)
    # the scalar form still works and is honoured
    res = a.simulate_waveform(t, g, initial_state=0.1)
    assert res.p_open[0] == pytest.approx(0.1, abs=1e-6)


# ======================================================================== P0-2
# rise_ms == clear_ms made exp(-t/c) - exp(-t/r) identically zero, the normaliser was
# clamped to 1e-6, and every sample came out 0 uM for a requested 1000 uM.

def test_p0_2_equal_time_constants_give_the_requested_peak():
    """ANCHOR. The alpha-function limit must reach `peak_um`.

    Tolerance 0.1%: the peak of (t/r)exp(-t/r) sits exactly at t = r, and the grid below
    samples that point, so the only error is float rounding. The defect produced 0.0,
    which no tolerance admits.
    """
    t = np.linspace(0.0, 50.0, 5001)          # 0.01 ms spacing; t = 1.0 is on the grid
    g = synaptic_transient(t, peak_um=1000.0, rise_ms=1.0, clear_ms=1.0)
    assert g.max() == pytest.approx(1000.0, rel=1e-3)
    assert np.argmax(g) == int(round(1.0 / 0.01))   # peaks at t = rise_ms


@pytest.mark.parametrize("ratio", [1.0, 1.0 + 1e-13, 1.000001, 1.001, 1.1, 10.0])
def test_p0_2_amplitude_is_continuous_across_the_limit(ratio):
    """ANCHOR. No discontinuity as clear_ms/rise_ms passes through 1.

    0.5% because the sampled peak of the two-exponential form lands slightly off the
    analytic maximum on a finite grid; the defect was a 100% error at exactly 1.0.
    """
    t = np.linspace(0.0, 60.0, 6001)
    g = synaptic_transient(t, peak_um=1000.0, rise_ms=1.0, clear_ms=1.0 * ratio)
    assert g.max() == pytest.approx(1000.0, rel=5e-3)
    assert np.all(g >= -1e-9)


def test_p0_2_guard_is_relative_not_an_equality_test():
    """The near-equal case was ALWAYS fine; a `c == r` guard would miss the real risk.

    Permanent record: if someone "simplifies" the tolerance back to equality, the
    1 + 1e-13 case above starts returning zeros again while this one still passes.
    """
    t = np.linspace(0.0, 60.0, 6001)
    near = synaptic_transient(t, peak_um=1000.0, rise_ms=0.999, clear_ms=1.0)
    exact = synaptic_transient(t, peak_um=1000.0, rise_ms=1.0, clear_ms=1.0)
    assert near.max() == pytest.approx(1000.0, rel=1e-3)
    # shapes agree closely either side of the limit
    assert float(np.abs(near - exact).max()) < 2.0      # of a 1000 uM amplitude


def test_p0_2_biexponential_clearance_reduces_to_one_component():
    """§2.6: the new two-clearance form must not change the single-clearance result.

    Compared as SHAPES (each normalised to its own peak) because the mixture normalises to
    its sampled maximum while `synaptic_transient` normalises analytically -- a ~0.02%
    amplitude difference on this grid, documented in the function.
    """
    t = np.linspace(0.0, 60.0, 4000)
    mixed = synaptic_transient_biexp_clearance(
        t, peak_um=1000.0, rise_ms=0.1, clear_fast_ms=1.0, clear_slow_ms=1.0,
        weight_fast=1.0)
    single = synaptic_transient(t, peak_um=1000.0, rise_ms=0.1, clear_ms=1.0)
    assert float(np.abs(mixed / mixed.max() - single / single.max()).max()) < 1e-9


def test_p0_2_biexponential_clearance_slows_the_decay():
    """A slow component must lengthen the tail; otherwise the weights do nothing."""
    t = np.linspace(0.0, 120.0, 6000)
    fast = synaptic_transient_biexp_clearance(
        t, peak_um=1.0, rise_ms=0.1, clear_fast_ms=1.0, clear_slow_ms=1.0, weight_fast=1.0)
    slow = synaptic_transient_biexp_clearance(
        t, peak_um=1.0, rise_ms=0.1, clear_fast_ms=1.0, clear_slow_ms=30.0,
        weight_fast=0.7)
    _trapz = getattr(np, "trapezoid", None) or getattr(np, "trapz", None)
    assert float(_trapz(slow, t)) > 3.0 * float(_trapz(fast, t))
    with pytest.raises(ValueError, match="faster than"):
        synaptic_transient_biexp_clearance(t, clear_fast_ms=10.0, clear_slow_ms=1.0)


# ======================================================================== P0-3
# A failed integration was replaced by a plausible-looking trace: a tiled constant in
# Models B and C, and in Model A the instantaneous-Hill target, which deletes the
# relaxation the model exists to express.

@pytest.mark.parametrize("model", [
    pytest.param(KineticAllosteryModel(), id="kinetic_jw95"),
    pytest.param(ExtendedDesensitizationModel(), id="extended_desens"),
    pytest.param(OperationalScalarModel(), id="operational"),
])
def test_p0_3_solver_failure_raises(model, monkeypatch):
    """A non-converged solve must raise, with the solver's own message in the text."""
    import scipy.integrate

    class _Failed:
        success = False
        message = "SENTINEL_SOLVER_MESSAGE"
        y = np.zeros((6, 3))

    monkeypatch.setattr(scipy.integrate, "solve_ivp", lambda *a, **k: _Failed())
    for mod_name in ("circuitpharm.models.kinetic_jw95",
                     "circuitpharm.models.extended_desens",
                     "circuitpharm.models.operational"):
        import importlib
        monkeypatch.setattr(importlib.import_module(mod_name), "solve_ivp",
                            lambda *a, **k: _Failed())

    t = _grid(300, 20.0)
    g = synaptic_transient(t, peak_um=1000.0, rise_ms=0.1, clear_ms=1.0)
    with pytest.raises(RuntimeError, match="SENTINEL_SOLVER_MESSAGE"):
        model.simulate_waveform(t, g)


def test_p0_3_simplex_drift_raises(monkeypatch):
    """Renormalising a drifted trajectory hides an integration failure behind a number."""
    import circuitpharm.models.kinetic_jw95 as kj

    class _Drifted:
        success = True
        message = ""

        def __init__(self, n):
            # probabilities summing to 1.5 -- far outside solver tolerance
            self.y = np.tile(np.array([[0.3], [0.3], [0.3], [0.3], [0.3]]), (1, n))

    t = _grid(50, 5.0)
    g = np.full_like(t, 0.4)
    monkeypatch.setattr(kj, "solve_ivp", lambda *a, **k: _Drifted(len(t)))
    with pytest.raises(RuntimeError, match="simplex"):
        KineticAllosteryModel().simulate_waveform(t, g)


# ======================================================================== P0-4
# decay_tau_ms fell back to the literal 15.0, which IS
# gabaa_kinetics.FIT_TARGETS["tau_ms"], so a failed estimate was indistinguishable from a
# perfect reproduction of the anchor.

def test_p0_4_biexponential_recovery():
    """ANCHOR. Recover a known two-component decay.

    Tolerances from the roadmap: weighted tau within 2%, each component within 5%. The
    synthetic trace is 0.7*exp(-t/15) + 0.3*exp(-t/70), so the weighted tau is
    (0.7*15 + 0.3*70)/1.0 = 31.5 ms. These are generous enough for a two-exponential fit
    on a finite window and far tighter than the 15.0 the defect returned.
    """
    t = np.linspace(0.0, 400.0, 4000)
    y = 0.70 * np.exp(-t / 15.0) + 0.30 * np.exp(-t / 70.0)
    fit = fit_biexponential_decay(t, y)
    assert fit.ok
    assert fit.r_squared > 0.999
    assert fit.tau_weighted_ms == pytest.approx(31.5, rel=0.02)
    assert fit.tau_fast_ms == pytest.approx(15.0, rel=0.05)
    assert fit.tau_slow_ms == pytest.approx(70.0, rel=0.05)
    assert fit.weight_fast == pytest.approx(0.70, rel=0.05)


def test_p0_4_single_exponential_is_recovered_too():
    """ANCHOR. A mono-exponential tail must give tau_weighted == its own tau.

    3% on a pure exponential fitted with a two-component form: the second amplitude goes
    to ~0 and the weighted mean collapses onto the single tau.
    """
    t = np.linspace(0.0, 300.0, 3000)
    fit = fit_biexponential_decay(t, np.exp(-t / 25.0))
    assert fit.ok
    assert fit.tau_weighted_ms == pytest.approx(25.0, rel=0.03)


def test_p0_4_failure_is_nan_never_the_fit_target():
    """PERMANENT RECORD. 15.0 must never be a fallback value.

    `gabaa_kinetics.FIT_TARGETS["tau_ms"]` is 15.0, so returning it on failure made a
    broken estimate look like a perfect anchor reproduction. A flat trace, a too-short
    trace and a NaN trace must all give all-NaN with r_squared 0.0.
    """
    from circuitpharm.gabaa_kinetics import FIT_TARGETS
    assert FIT_TARGETS["tau_ms"] == 15.0, "the sentinel this test guards against moved"

    t = np.linspace(0.0, 100.0, 500)
    for y in (np.ones_like(t),                      # flat
              np.zeros_like(t),                     # empty
              np.full_like(t, np.nan),              # NaN
              np.exp(-t / 10.0)[:4]):               # too short
        tt = t if y.shape == t.shape else t[:4]
        fit = fit_biexponential_decay(tt, y)
        assert not fit.ok
        assert fit.r_squared == 0.0
        assert fit.n_points == 0
        for field in (fit.tau_fast_ms, fit.tau_slow_ms, fit.tau_weighted_ms,
                      fit.weight_fast):
            assert np.isnan(field)
            assert field != 15.0


@pytest.mark.parametrize("model", [
    pytest.param(KineticAllosteryModel(), id="kinetic_jw95"),
    pytest.param(ExtendedDesensitizationModel(), id="extended_desens"),
    pytest.param(OperationalScalarModel(), id="operational"),
])
def test_p0_4_single_sample_trace_gives_nan_decay(model):
    t = np.array([0.0])
    g = np.array([0.4])
    res = model.simulate_waveform(t, g)
    assert np.isnan(res.decay_tau_ms)
    assert res.decay is not None and not res.decay.ok


# ======================================================================== P0-11
# current_pA defaulted to v_hold == e_rev == -70, so it returned an identically zero
# trace, and extract_electrophys_metrics worked around it with a second convention.

def test_p0_11_zero_driving_force_raises_rather_than_returning_zeros():
    t = _grid(400, 30.0)
    g = synaptic_transient(t, peak_um=1000.0, rise_ms=0.1, clear_ms=1.0)
    res = KineticAllosteryModel().simulate_waveform(t, g)
    with pytest.raises(ValueError, match="driving force"):
        res.current_pA(1.0, v_hold_mv=-70.0, e_cl_mv=-70.0)


def test_p0_11_current_scales_with_conductance_and_driving_force():
    """ANCHOR. I = g * P_open * (V_hold - E_Cl), exactly.

    Exact to float precision: this is a multiplication, so any tolerance beyond rounding
    would admit a wrong formula.
    """
    t = _grid(400, 30.0)
    g = synaptic_transient(t, peak_um=1000.0, rise_ms=0.1, clear_ms=1.0)
    res = KineticAllosteryModel().simulate_waveform(t, g)
    i1 = res.current_pA(1.0, v_hold_mv=-60.0, e_cl_mv=-75.0)
    i2 = res.current_pA(2.0, v_hold_mv=-60.0, e_cl_mv=-75.0)
    i3 = res.current_pA(1.0, v_hold_mv=-45.0, e_cl_mv=-75.0)
    assert np.allclose(i2, 2.0 * i1, rtol=0, atol=0)
    assert np.allclose(i3, 2.0 * i1, rtol=1e-12)
    assert np.allclose(i1, res.p_open * 15.0, rtol=1e-12)


def test_p0_11_chloride_reversal_matches_the_circuit_layer():
    """PERMANENT RECORD. One constant, two layers -- §5.3.

    A second copy of the GABA-A reversal is exactly how the SYNAPTIC_PULSE divergence
    happened, so the receptor layer's default is asserted equal to the circuit's.
    """
    from circuitpharm.cpg import E_REV
    assert DEFAULT_E_CL_MV == E_REV["gabaa"]


# =============================================== P1 contract (built alongside P0)

def test_observed_quantity_rejects_cross_observable_arithmetic():
    peak = ObservedQuantity(0.75, Observable.PEAK)
    eq = ObservedQuantity(0.17, Observable.EQUILIBRIUM)
    with pytest.raises(ObservableMismatch):
        peak / eq
    with pytest.raises(ObservableMismatch):
        peak - eq
    with pytest.raises(ObservableMismatch):
        peak > eq
    # same observable is fine
    assert ObservedQuantity(1.0, Observable.PEAK) / peak == pytest.approx(1 / 0.75)


def test_charge_quantity_requires_a_window():
    with pytest.raises(ValueError, match="window_ms"):
        ObservedQuantity(1.0, Observable.CHARGE)
    q = ObservedQuantity(1.0, Observable.CHARGE, window_ms=50.0)
    with pytest.raises(ObservableMismatch, match="different windows"):
        q / ObservedQuantity(1.0, Observable.CHARGE, window_ms=25.0)


def test_model_a_refuses_an_equilibrium_curve():
    """Model A has no desensitisation, so it has no equilibrium distinct from its peak."""
    a = OperationalScalarModel()
    with pytest.raises(ObservableMismatch, match="PEAK-native"):
        a.dose_response(np.array([1.0, 10.0]), observable=Observable.EQUILIBRIUM)
    assert a.native_observable is Observable.PEAK
    assert KineticAllosteryModel().native_observable is Observable.EQUILIBRIUM
    assert ExtendedDesensitizationModel().native_observable is Observable.EQUILIBRIUM


def test_peak_and_equilibrium_curves_differ_for_the_markov_models():
    """ANCHOR. The two observables are not relabellings of each other.

    At the provisional parameters the EQUILIBRIUM EC50 is ~6.3 uM and the PEAK EC50 is
    several times higher, because the stationary distribution has desensitisation
    competing throughout the rise while the peak is reached before it dominates. Asserting
    a factor > 2 rather than a precise value: the point is that they are different
    measurements, and the ratio itself is a prediction that P4 will pin at fitted
    parameters.
    """
    m = KineticAllosteryModel()
    c = np.logspace(-1, 4, 40)
    eq = m.dose_response(c)
    pk = m.peak_dose_response(c)

    def ec50(y):
        half = 0.5 * y.max()
        i = int(np.argmax(y >= half))
        if i == 0:
            return float(c[0])
        x0, x1, y0, y1 = c[i - 1], c[i], y[i - 1], y[i]
        f = (half - y0) / max(1e-18, y1 - y0)
        return float(np.exp(np.log(x0) + f * (np.log(x1) - np.log(x0))))

    assert pk.max() > eq.max()          # peak exceeds the desensitised stationary level
    assert ec50(pk) > 2.0 * ec50(eq)


# ======================================================================== P0-12
# apply_pam clamped its factor with max(f, 1e-6), admitting NAMs and 0.0 (which turns
# koff into 1e6x its value), while decomposition.py separately clamped to >= 1.0.

@pytest.mark.parametrize("model", [
    pytest.param(KineticAllosteryModel(), id="kinetic_jw95"),
    pytest.param(ExtendedDesensitizationModel(), id="extended_desens"),
])
def test_p0_12_negative_modulation_is_rejected_not_clamped(model):
    for bad in (0.0, 0.5, 0.999):
        with pytest.raises(ValueError, match="below 1.0"):
            model.apply_pam(bad) if isinstance(model, ExtendedDesensitizationModel) \
                else model.apply_pam(affinity_factor=bad)


def test_p0_12_desensitisation_denominator_is_guarded():
    """pam_desens_factor is fittable, and the old expression could zero the denominator."""
    m = ExtendedDesensitizationModel(pam_desens_factor=-10.0)
    # 1 + 0.1*(pf-1)*(-10) == 0 at pf == 2.0 -> infinite d_fast
    with pytest.raises(ValueError, match="denominator"):
        m.apply_pam(2.0)
    # and a legitimate value still works
    assert ExtendedDesensitizationModel().apply_pam(2.5).d_fast > 0.0


def test_p0_12_affinity_factor_is_the_equilibrium_ec50_shift():
    """ANCHOR. pam_factor divides K_d exactly, so it IS the equilibrium fold-shift.

    Exact to 1e-9: K_d = koff/kon and the modulation divides koff, so the shift is
    algebraic. This pins the semantics that P6-2 warns must not be confused with
    `gabaa_kinetics.calibrate_pam`'s PEAK-calibrated multiplier (2.945 for a nominal 2.5).
    """
    m = KineticAllosteryModel()
    assert m.apply_pam(affinity_factor=2.5).kd_um == pytest.approx(m.kd_um / 2.5, rel=1e-9)
