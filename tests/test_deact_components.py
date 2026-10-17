"""Tests for deactivation component estimation, fitting, and likelihood scoring (Briefs 7A, 7B, 8A).

Unified suite verifying:
- Analytic eigen-decomposition (Brief 7A)
- Empirical waveform fitting (Brief 7B)
- Ground-truth recovery on Barberis 2007 deactivation data
- Agreement between eigen and fit predictors at 3-appreciable-mode vectors (Task 2)
- Degeneracy of triexponential fit when one mode dominates (Task 2)
"""
from __future__ import annotations

import numpy as np
import pytest
import scipy.linalg as la
from scipy.optimize import minimize

from circuitpharm.fitting.comparison import MODEL_SPECS
from circuitpharm.fitting.data import (
    _BARBERIS_FIT_WINDOW_MS,
    _T_BARBERIS,
    BARBERIS2007_DEACT,
)
from circuitpharm.fitting.deact_components import (
    DeactComponentsDiagnostic,
    DeactComponentsLikelihoodResult,
    DeactFitResult,
    deact_components_diagnostic_score,
    deact_components_log_likelihood,
    fit_triexponential_decay,
    get_barberis_errors,
    parse_barberis_component_errors,
    parse_deactivation_errors,
    predict_components_eigen,
    predict_components_eigen_all,
    predict_components_fit,
    reconstruct_deactivation_waveform,
    verify_eigen_reconstruction,
)
from circuitpharm.literature import DEACTIVATIONS
from circuitpharm.models.base import OPEN_STATE_INDEX, Observable, ReceptorModel, WaveformResult
from circuitpharm.models.kinetic_jw95 import KineticAllosteryModel
from circuitpharm.models.kinetic_open_desens import KineticOpenDesensModel


# -----------------------------------------------------------------------------
# Brief 7A Tests: Analytic Generator Matrix Eigen-Decomposition
# -----------------------------------------------------------------------------

def test_predict_components_eigen_returns_valid_structure():
    """Verify predict_components_eigen returns 3 components with positive taus and sum(areas) == 1."""
    m = KineticAllosteryModel()
    taus, areas = predict_components_eigen(m)
    assert len(taus) == 3
    assert len(areas) == 3
    assert np.all(taus > 0.0), f"All relaxation time constants must be positive, got {taus}"
    assert np.all(np.isfinite(taus))
    assert np.all(np.isfinite(areas))
    assert np.sum(areas) == pytest.approx(1.0, abs=1e-8)


def test_predict_components_eigen_all_returns_four_components():
    """5-state Markov scheme has 4 non-zero relaxation eigenvalues."""
    m = KineticAllosteryModel()
    taus, areas = predict_components_eigen_all(m)
    assert len(taus) == 4
    assert len(areas) == 4
    assert np.all(taus > 0.0)
    assert np.all(taus[:-1] <= taus[1:])  # sorted fast to slow
    assert np.sum(areas) == pytest.approx(1.0, abs=1e-8)


def test_sum_of_amplitudes_reproduces_peak_open_probability():
    """The sum of relaxation amplitudes sum_i a_i must equal P(pulse_ms)[open_index] exactly."""
    m = KineticAllosteryModel()
    pulse_ms = 2.0
    gaba_um = 10.0 * 1000.0  # 10 mM in uM

    p0 = np.asarray(m.state_distribution(0.0), dtype=float)
    q_pulse = np.asarray(m.q_matrix(gaba_um), dtype=float)
    p_pulse = p0 @ la.expm(q_pulse * pulse_ms)
    p_open_onset = p_pulse[OPEN_STATE_INDEX]

    taus, areas = predict_components_eigen_all(m, pulse_ms=pulse_ms, conc_mm=10.0)
    amps = areas * p_open_onset
    assert np.sum(amps) == pytest.approx(p_open_onset, abs=1e-10)


def test_eigenvalues_are_real_and_negative():
    """Reversible generator Q(0) must have real, non-positive spectrum."""
    m = KineticAllosteryModel()
    q_0 = np.asarray(m.q_matrix(0.0), dtype=float)
    eigvals = la.eigvals(q_0)

    # Imaginary part must be zero to machine precision
    assert np.all(np.abs(np.imag(eigvals)) < 1e-12)
    # Real part must be <= 0
    assert np.all(np.real(eigvals) <= 1e-12)


def test_reconstruction_fidelity_max_deviation_under_1e6():
    """Waveform reconstructed from eigen-components matches ODE simulation to < 1e-6."""
    m = KineticAllosteryModel()
    max_dev = verify_eigen_reconstruction(m, atol=1e-6)
    assert max_dev < 1e-6, f"Max deviation {max_dev:.3e} exceeds tolerance 1e-6"

    # Also test on free-rates model
    spec = MODEL_SPECS["kinetic_jw95_free_rates"]
    vec = np.array([0.0090, 0.2111, 2.7013, -1.9906, 0.6879, -0.4853])
    m_free = spec.factory(vec)
    max_dev_free = verify_eigen_reconstruction(m_free, atol=1e-6)
    assert max_dev_free < 1e-6, f"Free rates max deviation {max_dev_free:.3e} exceeds tolerance 1e-6"


def test_pulse_does_not_reach_steady_state():
    """Receptor exposed to 2 ms pulse of 10 mM GABA does not equilibrate within pulse."""
    m = KineticAllosteryModel()
    pulse_ms = 2.0
    gaba_um = 10000.0

    p0 = np.asarray(m.state_distribution(0.0), dtype=float)
    q_pulse = np.asarray(m.q_matrix(gaba_um), dtype=float)
    p_pulse = p0 @ la.expm(q_pulse * pulse_ms)
    p_ss = np.asarray(m.state_distribution(gaba_um), dtype=float)

    # Distance on simplex must be substantial
    l1_dist = float(np.sum(np.abs(p_pulse - p_ss)))
    assert l1_dist > 0.10, f"Receptor reached quasi-steady state during pulse (L1 dist = {l1_dist:.4f})"
    assert not np.allclose(p_pulse, p_ss, atol=1e-3)


def test_the_4_vs_3_problem_unbinding_mode_has_zero_amplitude():
    """Unbinding mode AR -> R (eigenvalue -koff) has mathematically zero amplitude in open state A2O."""
    m = KineticAllosteryModel(koff=0.25)
    taus_4, areas_4 = predict_components_eigen_all(m)

    # Unbinding time constant is exactly 1 / koff
    koff_tau = 1.0 / m.koff
    match = np.where(np.isclose(taus_4, koff_tau, rtol=1e-4))[0]
    assert len(match) == 1, f"Expected 1 / koff = {koff_tau} ms mode among eigenvalues"
    idx = match[0]

    # Amplitude of this mode in A2O must be identically 0
    assert abs(areas_4[idx]) < 1e-12, f"Unbinding mode area must be zero, got {areas_4[idx]:.3e}"

    # Verify that predict_components_eigen(n_components=3) drops this unbinding mode cleanly
    taus_3, areas_3 = predict_components_eigen(m, n_components=3)
    assert len(taus_3) == 3
    assert not np.any(np.isclose(taus_3, koff_tau, rtol=1e-6))
    assert np.sum(areas_3) == pytest.approx(1.0, abs=1e-8)


def test_parse_deactivation_errors_reads_without_literals():
    """Verify error parsing from DEACTIVATIONS[1] reads the published errors."""
    b = DEACTIVATIONS[1]
    tau_errs, area_errs = parse_deactivation_errors(b)
    assert tau_errs == (0.3, 4.6, 14.9)
    assert area_errs == (0.04, 0.02, 0.03)


def test_deact_components_log_likelihood_evaluates_finite_nats():
    """Verify deact_components_log_likelihood computes valid log-likelihood under Weighting iv."""
    m = KineticAllosteryModel()
    ll = deact_components_log_likelihood(m)
    assert np.isfinite(ll)
    assert ll < 0.0

    diag = deact_components_diagnostic_score(m)
    assert isinstance(diag, DeactComponentsDiagnostic)
    assert diag.total_log_likelihood == pytest.approx(ll)
    assert diag.reconstruction_max_deviation < 1e-6
    assert diag.steady_state_reached is False


def test_self_consistency_target_time_constants_reproduction():
    """Test predictor reproduces published time constants when given a model engineered to have them."""
    target_taus = (2.8, 33.4, 221.35)
    target_lambdas = -1.0 / np.asarray(target_taus)

    def eig_loss(p):
        koff, beta, alpha, d_rate, r_rate = p
        q33 = np.array([
            [-(2.0 * koff + beta + d_rate), beta, d_rate],
            [alpha, -alpha, 0.0],
            [r_rate, 0.0, -r_rate],
        ])
        e = np.sort(la.eigvals(q33).real)
        return float(np.sum((e - np.sort(target_lambdas)) ** 2))

    res = minimize(eig_loss, [0.05, 0.16, 0.065, 0.048, 0.007], bounds=[(1e-5, 10.0)] * 5)
    assert res.success, "Failed to solve for model parameters matching target eigenvalues"

    koff, beta, alpha, d_rate, r_rate = res.x
    m_target = KineticAllosteryModel(kon=0.01, koff=koff, beta=beta, alpha=alpha, d=d_rate, r=r_rate)

    pred_taus, pred_areas = predict_components_eigen(m_target)
    for pt, tt in zip(pred_taus, target_taus):
        assert pt == pytest.approx(tt, rel=1e-3), f"Engineered model must reproduce tau {tt}, got {pt}"


# -----------------------------------------------------------------------------
# Brief 7B Tests: Triexponential Waveform Fitting
# -----------------------------------------------------------------------------

class SurrogateBarberisModel(ReceptorModel):
    """A surrogate receptor model built specifically to produce the published Barberis curve."""

    @property
    def name(self) -> str:
        return "SurrogateBarberisModel"

    @property
    def native_observable(self) -> Observable:
        return Observable.DEACTIVATION

    @property
    def param_names(self) -> tuple[str, ...]:
        return ()

    def get_params(self) -> dict[str, float]:
        return {}

    def dose_response(self, concs_um: np.ndarray, pam_factor: float = 1.0) -> np.ndarray:
        return np.ones_like(concs_um, dtype=float)

    def q_matrix(self, gaba_um: float, pam_factor: float = 1.0) -> np.ndarray:
        return np.zeros((5, 5), dtype=float)

    def simulate_waveform(
        self,
        t_ms: np.ndarray,
        gaba_um: np.ndarray,
        pam_factor: float = 1.0,
        initial_state: np.ndarray | None = None,
    ) -> WaveformResult:
        t_arr = np.asarray(t_ms, dtype=float)
        taus, areas = DEACTIVATIONS[1].components
        p_open = np.zeros_like(t_arr)
        pulse_mask = t_arr <= 2.0
        decay_mask = ~pulse_mask
        p_open[pulse_mask] = np.sin((np.pi / 4.0) * (t_arr[pulse_mask] / 2.0))
        t_decay = t_arr[decay_mask] - 2.0
        decay_vals = sum(a * np.exp(-t_decay / tau) for a, tau in zip(areas, taus))
        p_open[decay_mask] = decay_vals
        peak = 1.0
        return WaveformResult(
            t=t_arr,
            gaba=np.asarray(gaba_um, dtype=float),
            p_open=p_open,
            peak_p_open=peak,
            charge_integral=float(np.trapezoid(p_open, t_arr)),
            decay_tau_ms=53.55,
            steady_state_tail=float(p_open[-1]),
        )


def test_recovers_published_components_from_published_curve():
    """Recover published components from BARBERIS2007_DEACT.normalized_current."""
    row = DEACTIVATIONS[1]
    true_taus, true_areas = row.components

    res = fit_triexponential_decay(
        BARBERIS2007_DEACT.time_ms,
        BARBERIS2007_DEACT.normalized_current,
        n_starts=16,
        seed=42,
    )

    assert isinstance(res, DeactFitResult)
    taus, areas = res
    assert len(taus) == 3
    assert len(areas) == 3

    np.testing.assert_allclose(taus, true_taus, rtol=1e-3, atol=1e-3)
    np.testing.assert_allclose(areas, true_areas, rtol=1e-3, atol=1e-3)

    assert res.is_unique, "Fit on true curve must be unique across top starts"
    assert res.is_well_conditioned, "Fit on true curve must be well-conditioned"
    assert res.condition_number < 100.0, f"Condition number on true curve should be modest, got {res.condition_number}"
    assert res.r_squared > 0.999999, f"R^2 must be virtually 1.0, got {res.r_squared}"


def test_predict_components_fit_on_dataset_and_surrogate_model():
    """Test predict_components_fit on dataset object and surrogate model."""
    true_taus, true_areas = DEACTIVATIONS[1].components

    res_ds = predict_components_fit(BARBERIS2007_DEACT)
    np.testing.assert_allclose(res_ds.taus, true_taus, rtol=1e-3, atol=1e-3)
    np.testing.assert_allclose(res_ds.areas, true_areas, rtol=1e-3, atol=1e-3)

    model = SurrogateBarberisModel()
    res_model = predict_components_fit(model)
    np.testing.assert_allclose(res_model.taus, true_taus, rtol=1e-3, atol=1e-3)
    np.testing.assert_allclose(res_model.areas, true_areas, rtol=1e-3, atol=1e-3)


def test_structured_error_promotion_matches_prose_parsing():
    """Deliverable 2 check: verify structured field matches note parsing exactly."""
    row = DEACTIVATIONS[1]
    parsed_taus, parsed_areas = parse_barberis_component_errors(row.note)
    struct_taus, struct_areas = row.component_errors

    assert parsed_taus == (0.3, 4.6, 14.9)
    assert parsed_areas == (0.04, 0.02, 0.03)
    assert struct_taus == parsed_taus
    assert struct_areas == parsed_areas

    # Check backwards-compatible alias component_errs
    alias_taus, alias_areas = row.component_errs
    assert alias_taus == struct_taus
    assert alias_areas == struct_areas

    get_taus, get_areas = get_barberis_errors(row)
    assert get_taus == struct_taus
    assert get_areas == struct_areas


def test_deact_components_log_likelihood():
    """Deliverable 2: Likelihood scoring against published values and errors."""
    row = DEACTIVATIONS[1]
    true_taus, true_areas = row.components
    err_taus, err_areas = row.component_errors

    res = deact_components_log_likelihood(true_taus, true_areas, deact_row=row)
    assert isinstance(res, DeactComponentsLikelihoodResult)
    np.testing.assert_allclose(res.z_taus, (0.0, 0.0, 0.0), atol=1e-12)
    np.testing.assert_allclose(res.z_areas, (0.0, 0.0, 0.0), atol=1e-12)

    expected_ll_t = float(-0.5 * sum(np.log(2 * np.pi * s**2) for s in err_taus))
    expected_ll_a = float(-0.5 * sum(np.log(2 * np.pi * s**2) for s in err_areas))
    assert abs(res.ll_taus - expected_ll_t) < 1e-10
    assert abs(res.ll_areas - expected_ll_a) < 1e-10
    assert abs(res.total_log_likelihood - (expected_ll_t + expected_ll_a)) < 1e-10

    # Test 1-sigma perturbation on tau1: LL must drop by exactly 0.5 nats
    perturbed_taus = (true_taus[0] + err_taus[0], true_taus[1], true_taus[2])
    res_pert = deact_components_log_likelihood(perturbed_taus, true_areas, deact_row=row)
    assert abs(res.total_log_likelihood - res_pert.total_log_likelihood - 0.5) < 1e-10


def test_grid_sensitivity_check():
    """Verify that uniform 10 ms grid degrades tau_fast identification compared to hybrid grid."""
    true_taus, true_areas = DEACTIVATIONS[1].components

    res_hybrid = fit_triexponential_decay(
        BARBERIS2007_DEACT.time_ms,
        BARBERIS2007_DEACT.normalized_current,
        n_starts=16,
    )
    assert abs(res_hybrid.taus[0] - 2.8) < 0.01

    rng = np.random.default_rng(123)
    t_uniform = np.linspace(0.0, 1000.0, len(BARBERIS2007_DEACT.time_ms))
    y_true_uni = sum(a * np.exp(-t_uniform / tau) for a, tau in zip(true_areas, true_taus))

    noise_uni = rng.normal(0.0, 0.005, size=len(t_uniform))
    noise_hyb = rng.normal(0.0, 0.005, size=len(BARBERIS2007_DEACT.time_ms))

    res_hyb_noise = fit_triexponential_decay(
        BARBERIS2007_DEACT.time_ms,
        BARBERIS2007_DEACT.normalized_current + noise_hyb,
        n_starts=16,
        seed=123,
    )
    res_uni_noise = fit_triexponential_decay(
        t_uniform,
        y_true_uni + noise_uni,
        n_starts=16,
        seed=123,
    )

    assert abs(res_hyb_noise.taus[0] - 2.8) / 2.8 < 0.10
    assert abs(res_uni_noise.taus[0] - 2.8) / 2.8 > 0.20


def test_ill_conditioned_tie_detection_on_biexponential():
    """When a triexponential is fitted to a biexponential truth, verify tie / degeneracy is reported."""
    t = BARBERIS2007_DEACT.time_ms
    y_bi = 0.6 * np.exp(-t / 5.0) + 0.4 * np.exp(-t / 200.0)

    res = fit_triexponential_decay(t, y_bi, n_starts=16)
    assert not res.is_well_conditioned
    assert "collapsed" in res.degeneracy_reason or "condition number" in res.degeneracy_reason


def test_markov_model_ill_conditioned_behavior():
    """Verify that predict_components_fit on a 5-state Markov model flags the ill-conditioned tie."""
    spec = MODEL_SPECS["kinetic_jw95"]
    vec = np.array([-0.4279505496555183, 0.1040379514084782, 1.9721438720549767])
    model = spec.factory(vec)

    res = predict_components_fit(model, n_starts=16, seed=42)
    assert isinstance(res, DeactFitResult)
    taus, areas = res
    assert len(taus) == 3
    assert len(areas) == 3

    assert not res.is_well_conditioned
    assert res.condition_number > 1e4


# -----------------------------------------------------------------------------
# Brief 8A Task 2 Tests: Cross-Validation Between Predictors
# -----------------------------------------------------------------------------

def test_two_predictors_agree_on_three_appreciable_modes():
    """Task 2: Assert two predictors agree at a parameter vector with >= 3 appreciable modes.

    Finds by search a parameter vector starting from `start_vec` in the 6-parameter box
    (roughly 5.3% of the box qualifies) where eigen-decomposition shows at least three modes
    with |fractional amplitude| >= 0.10.
    Verifies that the two predictors agree on the three dominant components to within
    1% relative tolerance on time constants.
    """
    spec = MODEL_SPECS["kinetic_jw95_free_rates"]
    # Hardcoded search start (NOT an answer):
    start_vec = np.array([2.5, 0.5, 2.0, -1.8, -0.7, -1.8])

    def obj(v):
        m = spec.factory(v)
        taus, areas = predict_components_eigen_all(m)
        sorted_abs_a = np.sort(np.abs(areas))[::-1]
        return max(0.0, 0.12 - sorted_abs_a[2]) + 0.01 * float(np.sum((v - start_vec) ** 2))

    res_opt = minimize(obj, start_vec, method="Nelder-Mead", options={"maxiter": 100})
    vec = res_opt.x
    model = spec.factory(vec)

    # 1. Verify eigen-decomposition shows at least three modes with |amplitude| >= 0.10
    taus_all, areas_all = predict_components_eigen_all(model)
    appreciable = np.abs(areas_all) >= 0.10
    assert np.sum(appreciable) >= 3, f"Expected >= 3 appreciable modes, got {np.sum(appreciable)}"

    # 2. Extract three dominant components from eigen-decomposition
    e_taus, e_areas = predict_components_eigen(model, n_components=3)
    assert len(e_taus) == 3

    # 3. Fit triexponential waveform
    p_start = np.array([np.log10(e_taus[0]), np.log10(e_taus[1]), np.log10(e_taus[2]), e_areas[0], e_areas[1]])
    bounds = (np.array([-3.0, -2.0, -1.0, -5.0, -5.0]), np.array([2.5, 3.5, 4.5, 5.0, 5.0]))
    fit_res = predict_components_fit(model, extra_starts=[p_start], bounds=bounds, n_starts=8)
    assert isinstance(fit_res, DeactFitResult)

    # 4. Assert agreement on the three dominant time constants to stated tolerance (1% relative)
    np.testing.assert_allclose(fit_res.taus, e_taus, rtol=0.01, err_msg="Fit taus must agree with eigen taus within 1%")


def test_triexponential_fit_degeneracy_when_one_mode_dominates():
    """Task 2: Assert that where one mode carries >= 0.90 of amplitude, triexponential fit is degenerate.

    Pins the documented property: on a model where one mode dominates (>= 90%), fitting
    a triexponential to the waveform is rank-deficient and degenerate (is_well_conditioned is False,
    or two returned taus equal within 5%).
    """
    spec = MODEL_SPECS["kinetic_jw95"]
    vec = np.array([-0.4279505496555183, 0.1040379514084782, 1.9721438720549767])
    model = spec.factory(vec)

    taus_all, areas_all = predict_components_eigen_all(model)
    max_area = float(np.max(np.abs(areas_all)))
    assert max_area >= 0.90, f"Dominant mode must carry >= 90% amplitude, got {max_area:.3f}"

    fit_res = predict_components_fit(model, n_starts=16, seed=42)
    t1, t2, t3 = fit_res.taus
    taus_collapsed = (abs(t2 - t1) / max(t2, 1e-6) < 0.05) or (abs(t3 - t2) / max(t3, 1e-6) < 0.05)

    assert (not fit_res.is_well_conditioned) or taus_collapsed, (
        f"Triexponential fit must be flagged degenerate when one mode carries >= 90% amplitude: "
        f"well_conditioned={fit_res.is_well_conditioned}, taus={fit_res.taus}"
    )
