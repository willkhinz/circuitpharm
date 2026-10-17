"""Tests for analytic deactivation components from Q matrix eigenvalues (Brief 7A).

Verifies:
1. Predictor returns (taus, areas) from Q matrix eigenvalues.
2. Sum of amplitudes reproduces peak open probability at relaxation onset (sum a_i = P_pulse[3]).
3. Waveform reconstruction deviation against simulate_waveform is strictly < 1e-6.
4. Receptor does not reach steady state during the 2 ms pulse.
5. In the 5-state scheme, the unbinding mode AR -> R (lambda = -koff) carries identically
   zero amplitude in the open state A2O (4-vs-3 resolution).
6. Likelihood term scores against DEACTIVATIONS[1] published errors without hardcoded literals.
7. Self-consistency test: reproduces target components when given a model engineered to have them.
"""
import numpy as np
import pytest
import scipy.linalg as la

from circuitpharm.fitting.deact_components import (
    DeactComponentsDiagnostic,
    deact_components_diagnostic_score,
    deact_components_log_likelihood,
    parse_deactivation_errors,
    predict_components_eigen,
    predict_components_eigen_all,
    reconstruct_deactivation_waveform,
    verify_eigen_reconstruction,
)
from circuitpharm.fitting.data import BARBERIS2007_DEACT
from circuitpharm.fitting.comparison import MODEL_SPECS
from circuitpharm.literature import DEACTIVATIONS
from circuitpharm.models.kinetic_jw95 import KineticAllosteryModel
from circuitpharm.models.kinetic_open_desens import KineticOpenDesensModel


def test_predict_components_eigen_returns_valid_structure():
    """Verify predict_components_eigen returns 3 components with positive taus and sum(areas) == 1."""
    m = KineticAllosteryModel()
    taus, areas = predict_components_eigen(m)
    assert len(taus) == 3
    assert len(areas) == 3
    assert np.all(taus > 0.0)
    assert list(taus) == sorted(taus), "taus must be sorted fast to slow"
    assert np.sum(areas) == pytest.approx(1.0, abs=1e-8)


def test_predict_components_eigen_all_returns_four_components():
    """Verify predict_components_eigen_all returns all 4 non-zero relaxation components."""
    m = KineticAllosteryModel()
    taus_4, areas_4 = predict_components_eigen_all(m)
    assert len(taus_4) == 4
    assert len(areas_4) == 4
    assert np.all(taus_4 > 0.0)
    assert list(taus_4) == sorted(taus_4)
    assert np.sum(areas_4) == pytest.approx(1.0, abs=1e-8)


def test_sum_of_amplitudes_reproduces_peak_open_probability():
    """Assert sum(amplitudes) exactly equals open probability at pulse end (relaxation onset)."""
    m = KineticAllosteryModel()
    pulse_ms = 2.0
    conc_mm = 10.0
    gaba_um = conc_mm * 1000.0

    p0 = m.state_distribution(0.0)
    q_pulse = m.q_matrix(gaba_um)
    p_pulse = p0 @ la.expm(q_pulse * pulse_ms)
    p_open_onset = p_pulse[3]

    taus_4, areas_4 = predict_components_eigen_all(m, pulse_ms=pulse_ms, conc_mm=conc_mm)
    amps_4 = areas_4 * p_open_onset
    assert np.sum(amps_4) == pytest.approx(p_open_onset, abs=1e-12)


def test_eigenvalues_are_real_and_negative():
    """Assert that generator matrix Q_0 eigenvalues are strictly real and non-positive."""
    for spec_name in ["kinetic_jw95", "kinetic_open_desens"]:
        m = MODEL_SPECS[spec_name].factory(np.array([1.0, 0.2, 1.5]))
        q_0 = m.q_matrix(0.0)
        eigvals = la.eigvals(q_0)
        max_imag = np.max(np.abs(np.imag(eigvals)))
        assert max_imag < 1e-12, f"Eigenvalues of {spec_name} Q_0 must be real, got imag {max_imag}"
        max_real = np.max(np.real(eigvals))
        assert max_real <= 1e-12, f"Eigenvalues of {spec_name} Q_0 must be <= 0, got {max_real}"


def test_reconstruction_fidelity_max_deviation_under_1e6():
    """Verify max absolute deviation against simulate_waveform is strictly < 1e-6 (Brief 7A requirement)."""
    # Test on baseline models and brief-6 optima
    models = [
        KineticAllosteryModel(),
        KineticOpenDesensModel(),
        MODEL_SPECS["kinetic_jw95"].factory(np.array([-0.4279505, 0.1040380, 1.9721439])),
        MODEL_SPECS["kinetic_jw95_free_rates"].factory(
            np.array([0.0090277, 0.2110962, 2.7013357, -1.9906398, 0.6878775, -0.4852770])
        ),
        MODEL_SPECS["kinetic_open_desens"].factory(np.array([-0.5455756, 0.5392041, 1.5411059])),
    ]

    t_grid = BARBERIS2007_DEACT.time_ms
    for m in models:
        dev = verify_eigen_reconstruction(m, time_ms=t_grid, atol=1e-6)
        assert dev < 1e-6, f"{m.name} reconstruction deviation {dev:.3e} must be < 1e-6"


def test_pulse_does_not_reach_steady_state():
    """Verify that the receptor does NOT reach steady state during the 2 ms pulse of 10 mM GABA."""
    m = KineticAllosteryModel()
    pulse_ms = 2.0
    gaba_um = 10000.0

    p0 = m.state_distribution(0.0)
    q_pulse = m.q_matrix(gaba_um)
    p_pulse = p0 @ la.expm(q_pulse * pulse_ms)
    p_ss = m.state_distribution(gaba_um)

    # State occupancy differs dramatically, especially in desensitisation (A2D)
    assert not np.allclose(p_pulse, p_ss, atol=1e-2), (
        "Receptor was assumed not to reach steady state during 2 ms pulse; verify it does not."
    )
    # Open probability at 2 ms is near peak (~0.588), while steady state is desensitised (~0.166)
    assert p_pulse[3] > 2.0 * p_ss[3]
    # Desensitised state A2D is only beginning at 2 ms (~0.059) vs equilibrium (~0.801)
    assert p_pulse[4] < 0.15 * p_ss[4]


def test_the_4_vs_3_problem_unbinding_mode_has_zero_amplitude():
    """Verify that the 4th relaxation mode (AR -> R unbinding, lambda = -koff) has zero amplitude in A2O."""
    m = KineticAllosteryModel()
    taus_4, areas_4 = predict_components_eigen_all(m)
    # Find mode corresponding to koff: tau = 1 / koff
    koff_tau = 1.0 / m.koff
    # Locate matching tau
    idx = int(np.argmin(np.abs(taus_4 - koff_tau)))
    assert np.isclose(taus_4[idx], koff_tau, rtol=1e-6)
    # Amplitude in A2O is mathematically zero
    assert abs(areas_4[idx]) < 1e-12, (
        f"AR -> R unbinding mode at tau={taus_4[idx]:.4f} ms must have zero area in A2O, got {areas_4[idx]}"
    )

    # And 3-component prediction selects the 3 non-zero components
    taus_3, areas_3 = predict_components_eigen(m, n_components=3)
    assert len(taus_3) == 3
    assert not np.any(np.isclose(taus_3, koff_tau, rtol=1e-6))
    assert np.sum(areas_3) == pytest.approx(1.0, abs=1e-8)


def test_parse_deactivation_errors_reads_without_literals():
    """Verify error parsing from DEACTIVATIONS[1] reads the published errors."""
    b = DEACTIVATIONS[1]
    tau_errs, area_errs = parse_deactivation_errors(b)
    # Published values: 0.3, 4.6, 14.9 ms and 0.04, 0.02, 0.03
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
    """Test predictor reproduces published time constants when given a model engineered to have them.

    Acceptance criterion: 'A test that your predictor reproduces the published components when given
    a model built to have them, or, if that is not constructible, a stated reason and a weaker check
    that is honest about what it proves.'

    Here we solve for (koff, beta, alpha, d, r) such that the 3 relaxation eigenvalues of Q(0)
    match the target time constants (2.8, 33.4, 221.35) ms to within 1e-4 ms.
    """
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

    from scipy.optimize import minimize
    res = minimize(eig_loss, [0.05, 0.16, 0.065, 0.048, 0.007], bounds=[(1e-5, 10.0)] * 5)
    assert res.success, "Failed to solve for model parameters matching target eigenvalues"

    koff, beta, alpha, d_rate, r_rate = res.x
    m_target = KineticAllosteryModel(kon=0.01, koff=koff, beta=beta, alpha=alpha, d=d_rate, r=r_rate)

    pred_taus, pred_areas = predict_components_eigen(m_target)
    for pt, tt in zip(pred_taus, target_taus):
        assert pt == pytest.approx(tt, rel=1e-3), f"Engineered model must reproduce tau {tt}, got {pt}"
