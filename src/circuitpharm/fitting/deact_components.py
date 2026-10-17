"""Analytic prediction of deactivation components from generator matrix eigenvalues.

Agent Brief 7A:
Predicts deactivation relaxation time constants and component amplitudes analytically
from the Q matrix eigenvalues and eigenvectors at zero agonist, for a step application
protocol (Barberis et al. 2007: 2 ms pulse of 10 mM GABA in alpha1beta2gamma2S).

The three steps:
1. Initial condition:
   Starting from equilibrium at zero agonist (P(0) = state_distribution(0.0)),
   the receptor is exposed to a square pulse of agonist (conc_mm, default 10.0 mM = 10,000 uM)
   for duration pulse_ms (default 2.0 ms). Because the generator Q_pulse is constant during the pulse,
   the occupancy at the end of the pulse is obtained exactly via the matrix exponential:
       P(pulse_ms) = P(0) @ expm(Q_pulse * pulse_ms)
   The pulse is brief and saturating; the receptor does not reach steady state during it.

2. Relaxation:
   At zero agonist (t >= pulse_ms), the generator matrix Q_0 = q_matrix(0.0) governs relaxation:
       dP/dt = P Q_0
   We compute the eigen-decomposition Q_0 = V diag(lambda) V^{-1}.
   The zero eigenvalue corresponds to the stationary absorbing state (all channels in R).
   The non-zero eigenvalues lambda_i < 0 yield relaxation time constants:
       tau_i = -1 / lambda_i

3. Observable:
   Current is proportional to open probability (occupancy of A2O, open_index=3).
   The open probability relaxation is:
       p_open(t) = sum_i a_i exp(lambda_i t) = sum_i a_i exp(-t / tau_i)
   where amplitude a_i is the projection of the initial state onto eigenvector i and coordinate 3:
       a_i = (P(pulse_ms) @ V)_i * (V^{-1})_{i, open_index}
   The sum of all amplitudes sum_i a_i reproduces the open probability at relaxation onset,
   p_open(0) = P(pulse_ms)[open_index].
   Normalised component areas are:
       area_i = a_i / sum_j a_j = a_i / P(pulse_ms)[open_index]
   such that sum_i area_i = 1.0, exactly matching the normalisation of BARBERIS2007_DEACT.

4. The 4-vs-3 problem:
   A 5-state Markov scheme has 4 non-zero relaxation eigenvalues. However, because agonist is zero,
   no transitions occur from R -> AR or AR -> A2R. The unbinding mode AR -> R (eigenvalue lambda = -koff)
   has identically zero projection onto the {A2R, A2O, A2D} coupled subsystem.
   Therefore, the amplitude of the fourth mode in the open state A2O is mathematically zero (a_4 = 0).
   Reporting the 3 largest components by absolute amplitude retains 100% of the observable signal
   and discards a mode with zero observable current.
"""
from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Sequence

import numpy as np
import scipy.linalg as la

from ..literature import DEACTIVATIONS, Deactivation
from ..models.base import OPEN_STATE_INDEX, ReceptorModel


@dataclass(frozen=True)
class DeactComponentsDiagnostic:
    """Detailed score breakdown for deactivation components likelihood (Weighting iv)."""
    total_log_likelihood: float
    tau_log_likelihood: float
    area_log_likelihood: float
    predicted_taus: tuple[float, ...]
    predicted_areas: tuple[float, ...]
    target_taus: tuple[float, ...]
    target_areas: tuple[float, ...]
    tau_residuals_sigma: tuple[float, ...]
    area_residuals_sigma: tuple[float, ...]
    tau_w_predicted: float
    tau_w_residual_sigma: float
    reconstruction_max_deviation: float
    steady_state_reached: bool


def predict_components_eigen_all(
    model: ReceptorModel,
    *,
    pulse_ms: float = 2.0,
    conc_mm: float = 10.0,
    open_index: int = OPEN_STATE_INDEX,
) -> tuple[np.ndarray, np.ndarray]:
    """Compute all non-zero relaxation components (taus, areas) from the generator matrix.

    Parameters
    ----------
    model : ReceptorModel
        Markov kinetic model providing `q_matrix` and `state_distribution`.
    pulse_ms : float
        Pulse duration in ms (default 2.0 ms for Barberis 2007 brief pulse).
    conc_mm : float
        GABA concentration during pulse in mM (default 10.0 mM).
    open_index : int
        State index for open channel A2O (default 3).

    Returns
    -------
    taus : np.ndarray
        Array of all 4 relaxation time constants in ms, sorted ascending.
    areas : np.ndarray
        Array of all 4 normalised component areas (summing to 1.0).
    """
    if not hasattr(model, "q_matrix") or not hasattr(model, "state_distribution"):
        raise TypeError(
            f"{type(model).__name__} does not provide q_matrix and state_distribution; "
            f"analytic deactivation components require a discrete Markov generator."
        )

    # 1. Initial condition: 2 ms pulse of 10 mM GABA starting from equilibrium at 0 agonist
    gaba_pulse_um = float(conc_mm) * 1000.0
    p0 = np.asarray(model.state_distribution(0.0), dtype=float)
    q_pulse = np.asarray(model.q_matrix(gaba_pulse_um), dtype=float)
    p_pulse = p0 @ la.expm(q_pulse * float(pulse_ms))

    # Numerical hygiene: ensure state probabilities remain on simplex
    p_pulse = np.clip(p_pulse, 0.0, 1.0)
    p_pulse = p_pulse / np.sum(p_pulse)
    p_open_onset = float(p_pulse[open_index])

    # 2. Relaxation: decompose Q at zero agonist
    q_0 = np.asarray(model.q_matrix(0.0), dtype=float)
    eigvals, v = la.eig(q_0)

    # Check: eigenvalues must be real and non-positive within numerical tolerance
    max_imag_val = float(np.max(np.abs(np.imag(eigvals))))
    if max_imag_val > 1e-7:
        raise ValueError(
            f"Q_0 eigenvalues contain significant imaginary parts (max |imag| = {max_imag_val:.2e}); "
            f"reversible Markov scheme generators must have real spectrum. Check conditioning."
        )

    eigvals_real = np.real(eigvals)
    max_pos = float(np.max(eigvals_real))
    if max_pos > 1e-10:
        raise ValueError(
            f"Q_0 has positive eigenvalue ({max_pos:.2e}); generator must be negative semi-definite."
        )

    # Identify non-zero eigenvalues (stationary mode has lambda = 0)
    zero_tol = 1e-10
    non_zero = np.abs(eigvals_real) > zero_tol
    if np.sum(non_zero) != len(eigvals_real) - 1:
        # Stationary direction should be unique (absorbing state R)
        non_zero = np.argsort(eigvals_real)[:-1]

    lambda_nz = eigvals_real[non_zero]
    taus = -1.0 / lambda_nz

    # 3. Observable: projection onto open state coordinate
    v_nz = np.real(v)
    v_inv = la.inv(v_nz)

    # a_i = (P(pulse_ms) @ V)_i * (V^{-1})_{i, open_index}
    proj_pulse = p_pulse @ v_nz
    amps = proj_pulse[non_zero] * v_inv[non_zero, open_index]

    # Check: sum of amplitudes must reproduce peak open probability at onset
    sum_amps = float(np.sum(amps))
    if not np.isclose(sum_amps, p_open_onset, atol=1e-8, rtol=1e-6):
        raise ValueError(
            f"Sum of eigen-amplitudes ({sum_amps:.8e}) does not reproduce open probability "
            f"at pulse end ({p_open_onset:.8e})."
        )

    # Normalise areas to peak (p_open at relaxation onset)
    if p_open_onset > 0.0:
        areas = amps / p_open_onset
    else:
        areas = np.zeros_like(amps)

    # Sort fast-to-slow
    order = np.argsort(taus)
    return taus[order], areas[order]


def predict_components_eigen(
    model: ReceptorModel,
    *,
    pulse_ms: float = 2.0,
    conc_mm: float = 10.0,
    n_components: int = 3,
    open_index: int = OPEN_STATE_INDEX,
    return_all: bool = False,
) -> tuple[np.ndarray, np.ndarray]:
    """Predict deactivation (taus, areas) from generator matrix eigenvalues.

    By default returns the 3 physically active relaxation components, resolving the 4-vs-3 problem
    by dropping the unbinding mode AR -> R which carries mathematically zero amplitude in A2O.

    Parameters
    ----------
    model : ReceptorModel
        Kinetic model.
    pulse_ms : float
        Pulse duration in ms (default 2.0).
    conc_mm : float
        GABA concentration in mM (default 10.0).
    n_components : int
        Number of components to return (default 3 to match Barberis 2007).
    open_index : int
        Index of the open state (default 3).
    return_all : bool
        If True, returns all 4 non-zero components without reduction.

    Returns
    -------
    taus : np.ndarray
        Component time constants in ms, sorted ascending.
    areas : np.ndarray
        Component normalised areas (fractions).
    """
    taus_4, areas_4 = predict_components_eigen_all(
        model, pulse_ms=pulse_ms, conc_mm=conc_mm, open_index=open_index
    )
    if return_all:
        return taus_4, areas_4

    if n_components >= len(taus_4):
        return taus_4, areas_4

    # 4-vs-3 resolution:
    # Select the n_components with largest absolute amplitude/area.
    # In this topology, the unbinding mode AR -> R has amplitude ~ 0.0, so it is cleanly discarded.
    top_indices = np.argsort(-np.abs(areas_4))[:n_components]
    sort_taus = np.argsort(taus_4[top_indices])
    chosen = top_indices[sort_taus]
    return taus_4[chosen], areas_4[chosen]


def reconstruct_deactivation_waveform(
    taus: Sequence[float],
    amplitudes: Sequence[float],
    time_ms: np.ndarray,
) -> np.ndarray:
    """Reconstruct deactivation current waveform as a sum of exponentials.

    y(t) = sum_i a_i * exp(-t / tau_i)
    """
    t = np.asarray(time_ms, dtype=float)
    y = np.zeros_like(t)
    for a, tau in zip(amplitudes, taus):
        y += float(a) * np.exp(-t / float(tau))
    return y


def verify_eigen_reconstruction(
    model: ReceptorModel,
    time_ms: np.ndarray | None = None,
    *,
    pulse_ms: float = 2.0,
    conc_mm: float = 10.0,
    open_index: int = OPEN_STATE_INDEX,
    atol: float = 1e-6,
) -> float:
    """Reconstruct waveform from eigen-components and assert agreement with simulate_waveform.

    Parameters
    ----------
    model : ReceptorModel
        Kinetic model.
    time_ms : np.ndarray, optional
        Time grid for deactivation relaxation. Defaults to BARBERIS2007_DEACT.time_ms.
    pulse_ms : float
        Pulse duration in ms (default 2.0).
    conc_mm : float
        GABA concentration in mM (default 10.0).
    open_index : int
        Index of the open state (default 3).
    atol : float
        Maximum allowed absolute deviation (default 1e-6).

    Returns
    -------
    max_dev : float
        Maximum absolute deviation across the grid.
    """
    if time_ms is None:
        from ..fitting.data import BARBERIS2007_DEACT
        t_grid = np.asarray(BARBERIS2007_DEACT.time_ms, dtype=float)
    else:
        t_grid = np.asarray(time_ms, dtype=float)

    # 1. State at end of 2 ms pulse
    gaba_pulse_um = float(conc_mm) * 1000.0
    p0 = np.asarray(model.state_distribution(0.0), dtype=float)
    q_pulse = np.asarray(model.q_matrix(gaba_pulse_um), dtype=float)
    p_pulse = p0 @ la.expm(q_pulse * float(pulse_ms))
    p_pulse = np.clip(p_pulse, 0.0, 1.0)
    p_pulse = p_pulse / np.sum(p_pulse)

    # 2. All 4 components and unnormalized amplitudes
    taus_4, areas_4 = predict_components_eigen_all(
        model, pulse_ms=pulse_ms, conc_mm=conc_mm, open_index=open_index
    )
    p_open_onset = float(p_pulse[open_index])
    amps_4 = areas_4 * p_open_onset

    # 3. Reconstruct analytical waveform
    recon_p_open = reconstruct_deactivation_waveform(taus_4, amps_4, t_grid)

    # 4. Integrate with simulate_waveform starting from p_pulse at zero agonist
    wf = model.simulate_waveform(t_grid, np.zeros_like(t_grid), initial_state=p_pulse)
    sim_p_open = np.asarray(wf.p_open, dtype=float)

    max_dev = float(np.max(np.abs(recon_p_open - sim_p_open)))
    if max_dev > atol:
        raise AssertionError(
            f"Eigen-reconstruction deviation {max_dev:.3e} exceeds tolerance {atol:.3e}. "
            f"Check generator conditioning or stiff ODE solver tolerances."
        )
    return max_dev


def parse_deactivation_errors(
    record: Deactivation | None = None,
) -> tuple[tuple[float, ...], tuple[float, ...]]:
    """Extract published (tau_errors, area_errors) from Deactivation record.

    Reads `component_errs` structured field if populated; otherwise extracts from `note` prose.
    Never uses hardcoded numeric literals.
    """
    if record is None:
        record = DEACTIVATIONS[1]

    if getattr(record, "component_errs", None):
        c_errs = record.component_errs
        if len(c_errs) == 2 and len(c_errs[0]) > 0:
            return tuple(float(x) for x in c_errs[0]), tuple(float(x) for x in c_errs[1])

    # Regex parse from prose note:
    note = record.note
    # Example: "component errors are 0.3, 4.6, 14.9 ms and area errors 0.04, 0.02, 0.03."
    m_tau = re.search(r"component errors are ([0-9.,\s]+) ms", note)
    m_area = re.search(r"area errors ([0-9.,\s]+)", note)

    if not (m_tau and m_area):
        raise ValueError(
            f"Could not parse component and area errors from note of {record.label!r}: {note!r}"
        )

    tau_errs = tuple(float(x.strip()) for x in m_tau.group(1).split(",") if x.strip())
    # Strip any trailing punctuation from area string
    raw_area = m_area.group(1).rstrip(".;")
    area_errs = tuple(float(x.strip()) for x in raw_area.split(",") if x.strip())
    return tau_errs, area_errs


def deact_components_log_likelihood(
    model: ReceptorModel,
    *,
    deact_record: Deactivation | None = None,
    pulse_ms: float = 2.0,
    conc_mm: float = 10.0,
    score_areas: bool = True,
    score_tau_w: bool = False,
) -> float:
    """Evaluate Gaussian log-likelihood of predicted deactivation components against published record.

    Weighting (iv): Scores predicted components against published parameters at their published errors,
    avoiding the point-count leverage and noise profiling of reconstructed curves.

    Parameters
    ----------
    model : ReceptorModel
        Model to evaluate.
    deact_record : Deactivation, optional
        Literature record providing published values and errors. Defaults to DEACTIVATIONS[1]
        (Barberis 2007, brief pulse).
    pulse_ms : float
        Pulse length in ms.
    conc_mm : float
        GABA concentration in mM.
    score_areas : bool
        If True, scores both taus and areas (6 observations total).
    score_tau_w : bool
        If True, also scores derived weighted tau against published tau_w.

    Returns
    -------
    log_likelihood : float
        Total log-likelihood in nats.
    """
    if deact_record is None:
        deact_record = DEACTIVATIONS[1]

    target_taus, target_areas = deact_record.components
    target_taus = np.asarray(target_taus, dtype=float)
    target_areas = np.asarray(target_areas, dtype=float)

    tau_errs, area_errs = parse_deactivation_errors(deact_record)
    tau_errs = np.asarray(tau_errs, dtype=float)
    area_errs = np.asarray(area_errs, dtype=float)

    n_comp = len(target_taus)
    pred_taus, pred_areas = predict_components_eigen(
        model, pulse_ms=pulse_ms, conc_mm=conc_mm, n_components=n_comp
    )

    # Sort fast-to-slow to align components
    tau_diff = pred_taus - target_taus
    ll_tau = -0.5 * float(np.sum((tau_diff / tau_errs) ** 2 + np.log(2.0 * np.pi * (tau_errs ** 2))))

    ll_total = ll_tau
    if score_areas:
        area_diff = pred_areas - target_areas
        ll_area = -0.5 * float(np.sum((area_diff / area_errs) ** 2 + np.log(2.0 * np.pi * (area_errs ** 2))))
        ll_total += ll_area

    if score_tau_w:
        tau_w_pred = float(np.sum(pred_areas * pred_taus) / np.sum(pred_areas))
        tau_w_targ = float(deact_record.tau_w_ms)
        tau_w_err = float(deact_record.err_ms)
        tau_w_diff = tau_w_pred - tau_w_targ
        ll_tau_w = -0.5 * ((tau_w_diff / tau_w_err) ** 2 + np.log(2.0 * np.pi * (tau_w_err ** 2)))
        ll_total += ll_tau_w

    return float(ll_total)


def deact_components_diagnostic_score(
    model: ReceptorModel,
    *,
    deact_record: Deactivation | None = None,
    pulse_ms: float = 2.0,
    conc_mm: float = 10.0,
) -> DeactComponentsDiagnostic:
    """Full diagnostic evaluation of deactivation components with residuals and verification."""
    if deact_record is None:
        deact_record = DEACTIVATIONS[1]

    target_taus, target_areas = deact_record.components
    target_taus = tuple(float(x) for x in target_taus)
    target_areas = tuple(float(x) for x in target_areas)

    tau_errs, area_errs = parse_deactivation_errors(deact_record)
    tau_errs_arr = np.asarray(tau_errs, dtype=float)
    area_errs_arr = np.asarray(area_errs, dtype=float)

    pred_taus, pred_areas = predict_components_eigen(
        model, pulse_ms=pulse_ms, conc_mm=conc_mm, n_components=len(target_taus)
    )

    tau_diff = pred_taus - np.asarray(target_taus)
    area_diff = pred_areas - np.asarray(target_areas)

    tau_sig = tuple(float(x) for x in (tau_diff / tau_errs_arr))
    area_sig = tuple(float(x) for x in (area_diff / area_errs_arr))

    ll_tau = -0.5 * float(np.sum((tau_diff / tau_errs_arr) ** 2 + np.log(2.0 * np.pi * (tau_errs_arr ** 2))))
    ll_area = -0.5 * float(np.sum((area_diff / area_errs_arr) ** 2 + np.log(2.0 * np.pi * (area_errs_arr ** 2))))

    tau_w_pred = float(np.sum(pred_areas * pred_taus) / np.sum(pred_areas))
    tau_w_sig = (tau_w_pred - float(deact_record.tau_w_ms)) / float(deact_record.err_ms)

    # Reconstruction deviation check
    max_dev = verify_eigen_reconstruction(model, pulse_ms=pulse_ms, conc_mm=conc_mm, atol=1e-5)

    # Check whether steady state was reached during pulse
    gaba_pulse_um = float(conc_mm) * 1000.0
    p0 = np.asarray(model.state_distribution(0.0), dtype=float)
    q_pulse = np.asarray(model.q_matrix(gaba_pulse_um), dtype=float)
    p_pulse = p0 @ la.expm(q_pulse * float(pulse_ms))
    p_ss = np.asarray(model.state_distribution(gaba_pulse_um), dtype=float)
    steady_state_reached = bool(np.allclose(p_pulse, p_ss, atol=1e-3, rtol=1e-2))

    return DeactComponentsDiagnostic(
        total_log_likelihood=ll_tau + ll_area,
        tau_log_likelihood=ll_tau,
        area_log_likelihood=ll_area,
        predicted_taus=tuple(float(x) for x in pred_taus),
        predicted_areas=tuple(float(x) for x in pred_areas),
        target_taus=target_taus,
        target_areas=target_areas,
        tau_residuals_sigma=tau_sig,
        area_residuals_sigma=area_sig,
        tau_w_predicted=tau_w_pred,
        tau_w_residual_sigma=float(tau_w_sig),
        reconstruction_max_deviation=max_dev,
        steady_state_reached=steady_state_reached,
    )
