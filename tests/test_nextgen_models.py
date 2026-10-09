"""Comprehensive bug checks and unit tests for next-generation GABA-A models and protocols."""
from __future__ import annotations

import numpy as np
import pytest

from circuitpharm.models import (
    ReceptorModel,
    WaveformResult,
    OperationalScalarModel,
    KineticAllosteryModel,
    ExtendedDesensitizationModel,
)
from circuitpharm.protocols import (
    synaptic_transient,
    pulse_train,
    ambient_with_spillover,
    extract_electrophys_metrics,
    decompose_modulation_gain,
    compare_phasic_tonic_mechanisms,
    evaluate_model_fit,
    find_discriminating_protocol,
)
from circuitpharm.fitting import (
    DOSE_RESPONSE_BENCHMARK,
    DEACTIVATION_BENCHMARK,
    compute_kinetic_objective,
    compute_profile_likelihood,
    compute_fisher_information_matrix,
    run_ensemble_mcmc,
)
from circuitpharm.dynamic_range import evaluate_dynamic_range
from circuitpharm.gabaa_kinetics import Scheme


# ==============================================================================
# 1. Model Protocol & Individual Implementations
# ==============================================================================

def test_operational_scalar_model_protocol_and_behavior():
    model = OperationalScalarModel(ec50_um=20.0, hill_n=1.5, po_max=0.75, s_max=2.5)
    assert isinstance(model, ReceptorModel)
    assert model.name == "OperationalScalar"
    assert "ec50_um" in model.param_names

    # Zero agonist -> zero open prob
    assert model.steady_state(0.0) == 0.0

    # At EC50 -> exactly half of po_max
    p_half = model.steady_state(20.0, pam_factor=1.0)
    assert abs(p_half - 0.375) < 1e-4

    # PAM leftward shift capped at s_max (2.5x)
    p_pam_2x = model.steady_state(10.0, pam_factor=2.0)
    p_pam_5x = model.steady_state(10.0, pam_factor=5.0)  # should cap at 2.5x
    p_pam_2_5x = model.steady_state(10.0, pam_factor=2.5)
    assert abs(p_pam_5x - p_pam_2_5x) < 1e-9
    assert p_pam_2x < p_pam_2_5x

    # Dose response monotonic
    concs = np.logspace(-1, 3, 20)
    dr = model.dose_response(concs)
    assert np.all(np.diff(dr) >= 0.0)
    assert dr[-1] <= 0.75

    # Waveform simulation
    t = np.linspace(0.0, 50.0, 100)
    gaba_t = np.where((t >= 5.0) & (t <= 6.0), 1000.0, 0.0)
    res = model.simulate_waveform(t, gaba_t)
    assert isinstance(res, WaveformResult)
    assert res.peak_p_open > 0.5
    assert res.charge_integral > 0.0
    assert np.all(res.p_open >= 0.0) and np.all(res.p_open <= 1.0)


def test_kinetic_jw95_model_protocol_and_behavior():
    model = KineticAllosteryModel()
    assert isinstance(model, ReceptorModel)
    assert model.name == "KineticAllostery_JW95"

    # Microscopic constants
    assert abs(model.kd_um - (model.koff / model.kon)) < 1e-9
    assert abs(model.gating_efficacy - (model.beta / model.alpha)) < 1e-9
    assert abs(model.desens_ratio - (model.d / model.r)) < 1e-9

    # State distribution simplex constraint
    dist = model.state_distribution(0.40)
    assert len(dist) == 5
    assert abs(np.sum(dist) - 1.0) < 1e-6
    assert np.all(dist >= 0.0)

    # Steady state matches open probability in distribution
    po = model.steady_state(0.40)
    assert abs(po - dist[3]) < 1e-9

    # Infinite agonist limit matches po_inf
    po_high = model.steady_state(1e7)
    assert abs(po_high - model.po_inf) < 1e-4

    # Scheme compatibility wrapper
    s = Scheme(kon=model.kon, koff=model.koff, beta=model.beta, alpha=model.alpha, d=model.d, r=model.r)
    model_from_s = s.to_model()
    assert abs(model_from_s.steady_state(0.40) - po) < 1e-9

    # Dynamic waveform integration
    t = np.linspace(0.0, 30.0, 100)
    gaba = synaptic_transient(t, peak_um=1000.0, rise_ms=0.1, clear_ms=1.0)
    res = model.simulate_waveform(t, gaba)
    assert res.peak_p_open > 0.4
    assert res.states is not None
    assert res.states.shape == (100, 5)
    # Check probability conservation at all time points
    sums = np.sum(res.states, axis=1)
    assert np.all(np.abs(sums - 1.0) < 1e-4)


def test_extended_desensitization_model_protocol_and_behavior():
    model = ExtendedDesensitizationModel()
    assert isinstance(model, ReceptorModel)
    assert model.name == "ExtendedDesensitization_DualD"

    # Dual-pathway ratios
    assert abs(model.total_desens_ratio - (model.desens_fast_ratio + model.desens_slow_ratio)) < 1e-9

    # State distribution simplex constraint (6 states)
    dist = model.state_distribution(1.0)
    assert len(dist) == 6
    assert abs(np.sum(dist) - 1.0) < 1e-6
    assert np.all(dist >= 0.0)

    # Infinite agonist limit
    po_high = model.steady_state(1e7)
    assert abs(po_high - model.po_inf) < 1e-4

    # Transient waveform simulation
    t = np.linspace(0.0, 50.0, 100)
    gaba = synaptic_transient(t, peak_um=2000.0, rise_ms=0.1, clear_ms=1.0)
    res = model.simulate_waveform(t, gaba)
    assert res.peak_p_open > 0.3
    assert res.states is not None
    assert res.states.shape == (100, 6)
    assert np.all(np.abs(np.sum(res.states, axis=1) - 1.0) < 1e-4)


# ==============================================================================
# 2. Waveforms & Protocols
# ==============================================================================

def test_waveform_generators():
    t = np.linspace(0.0, 100.0, 1000)

    # Synaptic transient
    tr = synaptic_transient(t, peak_um=1000.0, rise_ms=0.2, clear_ms=1.5, t_start=5.0)
    assert tr[0] == 0.0
    assert abs(np.max(tr) - 1000.0) < 5.0
    assert tr[-1] < 1.0

    # Pulse train
    train = pulse_train(t, freq_hz=50.0, n_pulses=3, peak_um=800.0, clear_ms=1.0, start_ms=10.0)
    # Peak should approach 800 uM
    assert np.max(train) >= 700.0
    # Between pulses it should drop
    assert np.min(train[t > 15.0]) < 100.0

    # Ambient spillover
    amb = ambient_with_spillover(t, ambient_um=0.40, spillover_peak_um=3.0, event_times_ms=[20.0, 60.0])
    assert amb[0] == 0.40
    assert np.max(amb) > 3.0


def test_electrophys_metrics_extraction():
    # 400 samples, not 100: the biexponential decay fit (P0-4) needs at least 6 points
    # inside the 90%->10% window, and a 1 ms clearance on a 50 ms span at 100 samples
    # does not supply them. A coarse grid now yields an honest NaN rather than 15.0.
    t = np.linspace(0.0, 50.0, 400)
    gaba = synaptic_transient(t, peak_um=1000.0, rise_ms=0.1, clear_ms=1.0)
    model = KineticAllosteryModel()
    res = model.simulate_waveform(t, gaba)
    # ONE driving-force convention, from models.base (P0-11); there is no
    # `driving_force_mv` argument any more.
    metrics = extract_electrophys_metrics(res, g_max_ns=2.0, v_hold_mv=-60.0,
                                          e_cl_mv=-75.0)

    assert metrics["peak_current_pA"] > 0.0
    assert metrics["charge_integral_fC"] > 0.0
    assert metrics["driving_force_mv"] == pytest.approx(15.0)
    assert metrics["charge_window_ms"] == pytest.approx(50.0)
    # decay may legitimately not fit on a given grid; if it does, it must be positive
    assert np.isnan(metrics["decay_tau_ms"]) or metrics["decay_tau_ms"] > 0.0


# ==============================================================================
# 3. Mechanistic Decomposition
# ==============================================================================

def test_mechanistic_decomposition_attribution():
    model = KineticAllosteryModel()
    decomp = compare_phasic_tonic_mechanisms(model, tonic_gaba_um=0.40, phasic_gaba_um=1000.0, pam_factor=2.50)

    tonic = decomp["tonic"]
    phasic = decomp["phasic"]

    # Shapley percentages must sum to 100%
    assert abs(tonic.shapley_pct_occupancy + tonic.shapley_pct_gating + tonic.shapley_pct_headroom - 100.0) < 1e-3
    assert abs(phasic.shapley_pct_occupancy + phasic.shapley_pct_gating + phasic.shapley_pct_headroom - 100.0) < 1e-3

    # Tonic condition must exhibit higher total fold-gain than phasic
    assert tonic.total_fold_gain > phasic.total_fold_gain
    assert tonic.total_fold_gain > 2.0
    assert phasic.total_fold_gain < 1.5


# ==============================================================================
# 4. Model Comparison & Optimal Experimental Design
# ==============================================================================

def test_model_comparison_and_oed():
    m_a = OperationalScalarModel()
    m_b = KineticAllosteryModel()
    m_c = ExtendedDesensitizationModel()

    # Generate synthetic observations with noise
    concs = np.array([0.1, 1.0, 10.0, 30.0, 100.0, 1000.0])
    y_true = m_b.dose_response(concs)

    comp = evaluate_model_fit([m_a, m_b, m_c], concs, y_true, measurement_noise_std=0.02)
    assert len(comp) == 3
    # Akaike weights sum to 1.0
    total_w = sum(c.akaike_weight for c in comp)
    assert abs(total_w - 1.0) < 1e-4

    # OED: discover discriminating protocol
    protocol = find_discriminating_protocol([m_a, m_b, m_c], noise_sigma=0.02)
    assert protocol.discrimination_score > 0.0
    assert len(protocol.predicted_responses) == 3
    assert len(protocol.falsification_boundaries) == 3


# ==============================================================================
# 5. Parameter Estimation & Identifiability
# ==============================================================================

def test_fitting_objective_and_identifiability():
    model = KineticAllosteryModel()
    x0 = np.array([model.kon, model.koff, model.beta, model.alpha, model.d, model.r])

    # Objective evaluates finite positive chi2
    cost = compute_kinetic_objective(x0)
    assert np.isfinite(cost)
    assert cost >= 0.0

    # Profile likelihood on kon
    prof = compute_profile_likelihood("kon", model, scan_factors=(0.8, 1.0, 1.2))
    assert prof.param_name == "kon"
    assert len(prof.profile_costs) == 3
    assert prof.ci_95_bounds[0] <= prof.ci_95_bounds[1]

    # Fisher Information Matrix
    fim, cond = compute_fisher_information_matrix(model)
    assert fim.shape == (6, 6)
    assert np.isfinite(cond)
    assert cond > 0.0


def test_lightweight_ensemble_mcmc():
    model = KineticAllosteryModel()
    x0 = np.array([model.kon, model.koff, model.beta, model.alpha, model.d, model.r])

    # Short MCMC run to verify sampling mechanics
    res = run_ensemble_mcmc(x0, n_walkers=14, n_steps=15, burn_in=5, seed=123)
    assert isinstance(res.flat_samples, np.ndarray)
    assert res.flat_samples.shape[1] == 6
    assert len(res.posterior_means) == 6
    assert 0.0 <= res.acceptance_fraction <= 1.0


# ==============================================================================
# 6. Three-Tiered Dynamic Range Hierarchy
# ==============================================================================

def test_dynamic_range_hierarchy():
    model = KineticAllosteryModel()
    eval_dr = evaluate_dynamic_range(model, ambient_gaba_um=0.40, cleft_peak_um=1000.0)

    # Tier 1: Asymptotic headroom must be large (> 100x)
    assert eval_dr.theoretical_asymptotic_headroom > 100.0

    # Tier 2: Reachable gains at finite s_max
    assert 2.50 in eval_dr.reachable_gains_by_smax
    assert eval_dr.reachable_gains_by_smax[2.50] > 1.5
    # Order monotonic with s_max
    assert eval_dr.reachable_gains_by_smax[1.25] < eval_dr.reachable_gains_by_smax[2.50]

    # Tier 3: Physiological charge ratio under waveforms
    assert eval_dr.physiological_charge_ratio > 1.0  # tonic advantage preserved
    assert not eval_dr.headroom_collapses
