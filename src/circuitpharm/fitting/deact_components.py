"""Deactivation component prediction, fitting, and likelihood scoring (Weighting iv).

Unified module (Brief 8A):
Reconciles the analytic eigen-decomposition predictor (Brief 7A) and the
triexponential waveform-fitting predictor (Brief 7B) for the Barberis et al. 2007
protocol (2 ms pulse of 10 mM GABA in alpha1beta2gamma2S).

Predictors:
1. `predict_components_eigen`: Analytic prediction from the infinitesimal generator
   eigenvalues and eigenvectors at zero agonist (Q_0). Resolves the 4-vs-3 problem
   by identifying the zero-amplitude unbinding mode AR -> R.
2. `predict_components_eigen_all`: Computes all 4 relaxation components without reduction.
3. `predict_components_fit`: Triexponential non-linear least squares fit to simulated
   decay waveforms, matching the experimental methodology of Barberis et al. 2007.
4. `fit_triexponential_decay`: Multi-start bounded optimizer with Latin Hypercube sampling,
   Jacobian SVD condition analysis, and rank-deficiency / tie detection.
5. `deact_components_log_likelihood`: Unified Gaussian log-likelihood term scoring predicted
   (taus, areas) against published parameters at published errors (Weighting iv), reading from
   `circuitpharm.literature.DEACTIVATIONS[1]` without literals.
"""
from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Any, Sequence

import numpy as np
import scipy.linalg as la
from scipy.optimize import least_squares
from scipy.stats import qmc

from ..literature import DEACTIVATIONS, Deactivation
from ..models.base import OPEN_STATE_INDEX, ReceptorModel
from .data import _BARBERIS_FIT_WINDOW_MS, _T_BARBERIS, BARBERIS2007_DEACT, DeactivationDataset


class DeactFitResult(tuple):
    """Result of triexponential deactivation decay fit.

    Inherits from tuple with length 2: (taus, areas), allowing unpacking:
        taus, areas = result
    while also exposing rich diagnostic attributes:
        result.taus: (tau1, tau2, tau3) in ms, ordered fastest-first
        result.areas: (a1, a2, a3) fractional amplitudes, matching taus
        result.spread_taus: max - min across top optima (decades or ms)
        result.spread_areas: max - min across top optima
        result.is_unique: whether top optima agree within tolerance
        result.is_well_conditioned: whether components are distinct and identifiable
        result.condition_number: condition number of the Jacobian at the optimum
        result.singular_values: singular values of the Jacobian at the optimum
        result.cost: half sum of squared residuals
        result.r_squared: coefficient of determination
        result.degeneracy_reason: diagnostic note if fit is ill-conditioned or tied
    """

    def __new__(
        cls,
        taus: tuple[float, float, float],
        areas: tuple[float, float, float],
        *,
        spread_taus: tuple[float, float, float] = (0.0, 0.0, 0.0),
        spread_areas: tuple[float, float, float] = (0.0, 0.0, 0.0),
        is_unique: bool = True,
        is_well_conditioned: bool = True,
        condition_number: float = 1.0,
        singular_values: tuple[float, ...] = (),
        cost: float = 0.0,
        r_squared: float = 1.0,
        degeneracy_reason: str = "",
        top_optima: list[dict[str, Any]] | None = None,
    ):
        instance = super().__new__(cls, (taus, areas))
        instance._taus = taus
        instance._areas = areas
        instance._spread_taus = spread_taus
        instance._spread_areas = spread_areas
        instance._is_unique = is_unique
        instance._is_well_conditioned = is_well_conditioned
        instance._condition_number = condition_number
        instance._singular_values = singular_values
        instance._cost = cost
        instance._r_squared = r_squared
        instance._degeneracy_reason = degeneracy_reason
        instance._top_optima = top_optima or []
        return instance

    @property
    def taus(self) -> tuple[float, float, float]:
        return self._taus

    @property
    def areas(self) -> tuple[float, float, float]:
        return self._areas

    @property
    def spread_taus(self) -> tuple[float, float, float]:
        return self._spread_taus

    @property
    def spread_areas(self) -> tuple[float, float, float]:
        return self._spread_areas

    @property
    def is_unique(self) -> bool:
        return self._is_unique

    @property
    def is_well_conditioned(self) -> bool:
        return self._is_well_conditioned

    @property
    def condition_number(self) -> float:
        return self._condition_number

    @property
    def singular_values(self) -> tuple[float, ...]:
        return self._singular_values

    @property
    def cost(self) -> float:
        return self._cost

    @property
    def r_squared(self) -> float:
        return self._r_squared

    @property
    def degeneracy_reason(self) -> str:
        return self._degeneracy_reason

    @property
    def top_optima(self) -> list[dict[str, Any]]:
        return self._top_optima

    @property
    def tau_weighted_ms(self) -> float:
        """Derived amplitude-weighted deactivation time constant: sum(a_i * tau_i)."""
        return float(sum(a * t for a, t in zip(self._areas, self._taus)))

    def __repr__(self) -> str:
        cond_str = "well-conditioned" if self._is_well_conditioned else f"ill-conditioned ({self._degeneracy_reason})"
        return (
            f"DeactFitResult(taus={tuple(round(t, 3) for t in self._taus)}, "
            f"areas={tuple(round(a, 3) for a in self._areas)}, "
            f"tau_w={self.tau_weighted_ms:.2f} ms, "
            f"cond={self._condition_number:.2e}, {cond_str})"
        )


class DeactComponentsLikelihoodResult(float):
    """Log-likelihood evaluation of predicted components against literature.

    Subclasses float so that it functions as a plain float log-likelihood in nats:
        ll = deact_components_log_likelihood(model)
        assert np.isfinite(ll)
        assert ll < 0.0
    while also exposing rich diagnostic attributes (matching Brief 7B):
        res.total_log_likelihood: float
        res.ll_taus: float
        res.ll_areas: float
        res.z_taus: tuple[float, ...]
        res.z_areas: tuple[float, ...]
        res.taus_pred: tuple[float, ...]
        res.taus_pub: tuple[float, ...]
        res.taus_err: tuple[float, ...]
        res.areas_pred: tuple[float, ...]
        res.areas_pub: tuple[float, ...]
        res.areas_err: tuple[float, ...]
        res.tau_w_pred: float
        res.tau_w_pub: float
        res.tau_w_err: float
        res.z_tau_w: float
    """

    def __new__(
        cls,
        total_log_likelihood: float,
        ll_taus: float,
        ll_areas: float,
        z_taus: tuple[float, ...],
        z_areas: tuple[float, ...],
        taus_pred: tuple[float, ...],
        taus_pub: tuple[float, ...],
        taus_err: tuple[float, ...],
        areas_pred: tuple[float, ...],
        areas_pub: tuple[float, ...],
        areas_err: tuple[float, ...],
        tau_w_pred: float,
        tau_w_pub: float,
        tau_w_err: float,
        z_tau_w: float,
    ):
        instance = super().__new__(cls, float(total_log_likelihood))
        instance.total_log_likelihood = float(total_log_likelihood)
        instance.ll_taus = float(ll_taus)
        instance.ll_areas = float(ll_areas)
        instance.z_taus = tuple(float(x) for x in z_taus)
        instance.z_areas = tuple(float(x) for x in z_areas)
        instance.taus_pred = tuple(float(x) for x in taus_pred)
        instance.taus_pub = tuple(float(x) for x in taus_pub)
        instance.taus_err = tuple(float(x) for x in taus_err)
        instance.areas_pred = tuple(float(x) for x in areas_pred)
        instance.areas_pub = tuple(float(x) for x in areas_pub)
        instance.areas_err = tuple(float(x) for x in areas_err)
        instance.tau_w_pred = float(tau_w_pred)
        instance.tau_w_pub = float(tau_w_pub)
        instance.tau_w_err = float(tau_w_err)
        instance.z_tau_w = float(z_tau_w)
        return instance

    def __repr__(self) -> str:
        return (
            f"DeactComponentsLikelihoodResult(total_ll={self.total_log_likelihood:.4f}, "
            f"ll_taus={self.ll_taus:.4f}, ll_areas={self.ll_areas:.4f}, "
            f"z_taus={tuple(round(z, 2) for z in self.z_taus)}, "
            f"z_areas={tuple(round(z, 2) for z in self.z_areas)})"
        )


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


def parse_barberis_component_errors(
    note: str,
) -> tuple[tuple[float, float, float], tuple[float, float, float]]:
    """Parse component errors and area errors from a Deactivation record note prose.

    Matches e.g.:
      'Triple exponential; component errors are 0.3, 4.6, 14.9 ms and area errors 0.04, 0.02, 0.03.'
    """
    pattern = (
        r"component errors are\s+([\d]+(?:\.[\d]+)?),\s*([\d]+(?:\.[\d]+)?),\s*([\d]+(?:\.[\d]+)?)\s*ms\s+"
        r"and area errors\s+([\d]+(?:\.[\d]+)?),\s*([\d]+(?:\.[\d]+)?),\s*([\d]+(?:\.[\d]+)?)"
    )
    m = re.search(pattern, note)
    if not m:
        m_tau = re.search(r"component errors are ([0-9.,\s]+) ms", note)
        m_area = re.search(r"area errors ([0-9.,\s]+)", note)
        if not (m_tau and m_area):
            raise ValueError(f"Could not parse component errors from note: {note!r}")
        tau_errs = tuple(float(x.strip()) for x in m_tau.group(1).split(",") if x.strip())
        raw_area = m_area.group(1).rstrip(".;")
        area_errs = tuple(float(x.strip()) for x in raw_area.split(",") if x.strip())
        return tau_errs, area_errs
    tau_errs = (float(m.group(1)), float(m.group(2)), float(m.group(3)))
    area_errs = (float(m.group(4)), float(m.group(5)), float(m.group(6)))
    return tau_errs, area_errs


def get_barberis_errors(
    record: Deactivation | None = None,
) -> tuple[tuple[float, ...], tuple[float, ...]]:
    """Retrieve published component errors from structured field or prose fallback."""
    if record is None:
        record = DEACTIVATIONS[1]
    for attr in ("component_errors", "component_errs"):
        val = getattr(record, attr, None)
        if val and len(val) == 2 and len(val[0]) > 0:
            return tuple(float(x) for x in val[0]), tuple(float(x) for x in val[1])
    return parse_barberis_component_errors(record.note)


def parse_deactivation_errors(
    record: Deactivation | None = None,
) -> tuple[tuple[float, ...], tuple[float, ...]]:
    """Alias for get_barberis_errors (Brief 7A name)."""
    return get_barberis_errors(record)


# -----------------------------------------------------------------------------
# 1. Analytic Eigen-Decomposition Predictor (Brief 7A)
# -----------------------------------------------------------------------------

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
        Array of all relaxation time constants in ms, sorted ascending.
    areas : np.ndarray
        Array of all normalised component areas (summing to 1.0).
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
        If True, returns all non-zero components without reduction.

    Returns
    -------
    taus : np.ndarray
        Component time constants in ms, sorted ascending.
    areas : np.ndarray
        Component normalised areas (fractions).
    """
    taus_all, areas_all = predict_components_eigen_all(
        model, pulse_ms=pulse_ms, conc_mm=conc_mm, open_index=open_index
    )
    if return_all:
        return taus_all, areas_all

    if n_components >= len(taus_all):
        return taus_all, areas_all

    top_indices = np.argsort(-np.abs(areas_all))[:n_components]
    sort_taus = np.argsort(taus_all[top_indices])
    chosen = top_indices[sort_taus]
    return taus_all[chosen], areas_all[chosen]


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

    # 2. All components and unnormalized amplitudes
    taus_all, areas_all = predict_components_eigen_all(
        model, pulse_ms=pulse_ms, conc_mm=conc_mm, open_index=open_index
    )
    p_open_onset = float(p_pulse[open_index])
    amps_all = areas_all * p_open_onset

    # 3. Reconstruct analytical waveform
    recon_p_open = reconstruct_deactivation_waveform(taus_all, amps_all, t_grid)

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


# -----------------------------------------------------------------------------
# 2. Triexponential Fitting Predictor (Brief 7B)
# -----------------------------------------------------------------------------

def fit_triexponential_decay(
    t_ms: np.ndarray,
    y: np.ndarray,
    *,
    n_starts: int = 16,
    seed: int = 42,
    spread_tol: float = 0.05,
    extra_starts: Sequence[np.ndarray] | None = None,
    bounds: tuple[np.ndarray, np.ndarray] | None = None,
) -> DeactFitResult:
    """Fit a normalized triexponential decay sum(a_i * exp(-t / tau_i)) to a time series.

    The model parameterisation enforces sum(a_i) = 1.0 (with a_3 = 1 - a_1 - a_2) and
    searches in log10(tau) to cover timescales spanning from 0.01 ms to 10,000 ms.

    Convergence and conditioning checks:
    - Multi-start optimisation via Latin Hypercube + warm starts.
    - Top-5 optima are tracked to detect ties and spread.
    - SVD of the Jacobian evaluates condition number and rank deficiency.
    - If two time constants collapse within 5% or an area is near 0, the fit is
      flagged as ill-conditioned (rank deficient / non-identifiable tie).
    """
    t_arr = np.asarray(t_ms, dtype=float)
    y_arr = np.asarray(y, dtype=float)
    if t_arr.shape != y_arr.shape or t_arr.ndim != 1 or t_arr.size < 10:
        raise ValueError(f"Invalid input arrays: shapes {t_arr.shape} vs {y_arr.shape}")
    if not (np.all(np.isfinite(t_arr)) and np.all(np.isfinite(y_arr))):
        raise ValueError("Non-finite values encountered in deactivation decay trace.")

    # Objective: residuals for p = [log10_tau1, log10_tau2, log10_tau3, a1, a2]
    # a3 = 1.0 - a1 - a2
    def resid(p):
        t1, t2, t3 = 10.0 ** p[0], 10.0 ** p[1], 10.0 ** p[2]
        a1, a2 = p[3], p[4]
        a3 = 1.0 - a1 - a2
        y_mod = a1 * np.exp(-t_arr / t1) + a2 * np.exp(-t_arr / t2) + a3 * np.exp(-t_arr / t3)
        return y_mod - y_arr

    if bounds is None:
        lo_b = np.array([-2.0, -1.0, 0.5, 0.0, 0.0])
        hi_b = np.array([1.5, 2.8, 4.2, 1.0, 1.0])
        bounds = (lo_b, hi_b)
    else:
        lo_b = np.asarray(bounds[0], dtype=float)
        hi_b = np.asarray(bounds[1], dtype=float)
        bounds = (lo_b, hi_b)

    sampler = qmc.LatinHypercube(d=5, seed=seed)
    sampled = lo_b + sampler.random(n_starts) * (hi_b - lo_b)
    for row in sampled:
        row[:3] = np.sort(row[:3])
        if row[3] >= 0 and row[4] >= 0 and row[3] + row[4] > 0.95:
            scale = 0.85 / (row[3] + row[4])
            row[3] *= scale
            row[4] *= scale

    warm_starts = [
        np.array([np.log10(2.8), np.log10(33.4), np.log10(221.35), 0.57, 0.23]),
        np.array([np.log10(5.0), np.log10(50.0), np.log10(300.0), 0.33, 0.33]),
        np.array([np.log10(1.0), np.log10(15.0), np.log10(500.0), 0.40, 0.30]),
    ]
    if extra_starts:
        warm_starts.extend(extra_starts)

    all_starts = np.vstack([warm_starts, sampled])

    optima = []
    for s0 in all_starts:
        s0_clamped = np.clip(s0, lo_b + 1e-6, hi_b - 1e-6)
        try:
            res = least_squares(
                resid,
                s0_clamped,
                bounds=bounds,
                method="trf",
                ftol=1e-12,
                xtol=1e-12,
                gtol=1e-12,
                max_nfev=15000,
            )
        except Exception:
            continue

        if not res.success or not np.all(np.isfinite(res.x)):
            continue

        p = res.x
        t1, t2, t3 = 10.0 ** p[0], 10.0 ** p[1], 10.0 ** p[2]
        a1, a2 = p[3], p[4]
        a3 = 1.0 - a1 - a2

        tau_raw = np.array([t1, t2, t3])
        area_raw = np.array([a1, a2, a3])
        perm = np.argsort(tau_raw)
        taus_ordered = tuple(float(tau_raw[i]) for i in perm)
        areas_ordered = tuple(float(area_raw[i]) for i in perm)

        s_vals = np.linalg.svd(res.jac, compute_uv=False)
        cond = float(s_vals[0] / s_vals[-1]) if s_vals[-1] > 1e-14 else float("inf")

        optima.append({
            "cost": float(res.cost),
            "taus": taus_ordered,
            "areas": areas_ordered,
            "cond": cond,
            "singular_values": tuple(float(x) for x in s_vals),
            "res": res,
        })

    if not optima:
        raise RuntimeError("None of the triexponential starts converged.")

    optima.sort(key=lambda item: item["cost"])
    best = optima[0]

    top = optima[: min(5, len(optima))]
    tau_mat = np.array([item["taus"] for item in top])
    area_mat = np.array([item["areas"] for item in top])

    spread_t = tuple(float(np.ptp(tau_mat[:, i])) for i in range(3))
    spread_a = tuple(float(np.ptp(area_mat[:, i])) for i in range(3))

    log_tau_mat = np.log10(np.maximum(tau_mat, 1e-12))
    log_spread_t = np.ptp(log_tau_mat, axis=0)
    is_unique = bool(np.all(log_spread_t <= spread_tol) and np.all(np.array(spread_a) <= spread_tol))

    t1, t2, t3 = best["taus"]
    a1, a2, a3 = best["areas"]
    merge_12 = abs(t2 - t1) / max(t2, 1e-6) < 0.05
    merge_23 = abs(t3 - t2) / max(t3, 1e-6) < 0.05
    zero_area = any(abs(a) < 0.01 for a in (a1, a2, a3))
    high_cond = best["cond"] > 1e4

    is_well_conditioned = not (merge_12 or merge_23 or zero_area or high_cond)

    deg_reasons = []
    if merge_12:
        deg_reasons.append(f"tau1 and tau2 collapsed ({t1:.2f} ms vs {t2:.2f} ms; a1+a2 unidentifiable tie)")
    if merge_23:
        deg_reasons.append(f"tau2 and tau3 collapsed ({t2:.2f} ms vs {t3:.2f} ms)")
    if zero_area:
        deg_reasons.append(f"vanishing amplitude component ({a1:.3f}, {a2:.3f}, {a3:.3f})")
    if high_cond:
        deg_reasons.append(f"Jacobian condition number {best['cond']:.2e} indicates rank deficiency")

    degeneracy_str = "; ".join(deg_reasons) if deg_reasons else "None"

    ss_res = 2.0 * best["cost"]
    ss_tot = float(np.sum((y_arr - np.mean(y_arr)) ** 2))
    r2 = 1.0 - (ss_res / ss_tot) if ss_tot > 1e-30 else 0.0

    return DeactFitResult(
        taus=best["taus"],
        areas=best["areas"],
        spread_taus=spread_t,
        spread_areas=spread_a,
        is_unique=is_unique,
        is_well_conditioned=is_well_conditioned,
        condition_number=best["cond"],
        singular_values=best["singular_values"],
        cost=best["cost"],
        r_squared=float(r2),
        degeneracy_reason=degeneracy_str,
        top_optima=top,
    )


def predict_components_fit(
    model: Any,
    *,
    pulse_ms: float = 2.0,
    conc_mm: float = 10.0,
    fit_window_ms: float = _BARBERIS_FIT_WINDOW_MS,
    grid: np.ndarray | None = None,
    n_starts: int = 16,
    seed: int = 42,
    uniform_grid: bool = False,
    n_pulse: int = 200,
    n_tail: int = 500,
    extra_starts: Sequence[np.ndarray] | None = None,
    bounds: tuple[np.ndarray, np.ndarray] | None = None,
) -> DeactFitResult:
    """Predict triexponential deactivation components (taus, areas) by fitting simulated waveform.

    Implements the experimental procedure of Barberis et al. 2007:
    1. Protocol simulation: 2 ms pulse of 10 mM GABA, followed by washout to 0 mM over fit_window_ms
       (default 1000.0 ms), computed via model.simulate_waveform.
    2. Peak extraction & normalisation: current decay is re-zeroed at peak open probability and
       divided by peak open probability so y(0) = 1.0, matching BARBERIS2007_DEACT.normalized_current.
    3. Grid sampling: evaluated on the dataset's hybrid log/linear grid (_T_BARBERIS) by default,
       or on a uniform grid if uniform_grid=True (sensitivity check).
    4. Triexponential fit: multi-start non-linear least squares fit with convergence and
       rank-deficiency detection, returning (taus, areas) ordered fastest-first.
    """
    if grid is not None:
        t_eval = np.asarray(grid, dtype=float)
    elif uniform_grid:
        t_eval = np.linspace(0.0, float(fit_window_ms), len(_T_BARBERIS))
    else:
        t_eval = _T_BARBERIS

    if isinstance(model, np.ndarray):
        y_eval = np.asarray(model, dtype=float)
        return fit_triexponential_decay(t_eval, y_eval, n_starts=n_starts, seed=seed, extra_starts=extra_starts, bounds=bounds)

    if isinstance(model, DeactivationDataset) or hasattr(model, "normalized_current"):
        ds_t = getattr(model, "time_ms", t_eval)
        ds_y = getattr(model, "normalized_current")
        if np.array_equal(ds_t, t_eval):
            y_eval = np.asarray(ds_y, dtype=float)
        else:
            y_eval = np.interp(t_eval, ds_t, ds_y)
        return fit_triexponential_decay(t_eval, y_eval, n_starts=n_starts, seed=seed, extra_starts=extra_starts, bounds=bounds)

    if not hasattr(model, "q_matrix"):
        from ..models.base import ObservableMismatch
        raise ObservableMismatch(
            f"predict_components_fit requires a Markov model with state transitions (`q_matrix`), "
            f"got {type(model).__name__}."
        )

    # 1. Simulate the protocol
    t_sim = np.unique(np.concatenate([
        np.linspace(0.0, float(pulse_ms), n_pulse),
        float(pulse_ms) + np.geomspace(1e-4, float(fit_window_ms), n_tail),
    ]))
    gaba_um = np.where(t_sim <= float(pulse_ms), float(conc_mm) * 1000.0, 0.0)

    res = model.simulate_waveform(t_sim, gaba_um)
    p_open = np.asarray(res.p_open, dtype=float)

    # 2. Extract and normalise decay
    i_peak = int(np.argmax(p_open))
    peak = float(p_open[i_peak])
    if not (np.isfinite(peak) and peak > 0.0):
        raise ValueError(
            f"Model peak open probability is non-positive or non-finite ({peak!r}); decay cannot be normalised."
        )

    t_decay = t_sim[i_peak:] - t_sim[i_peak]
    y_decay = p_open[i_peak:] / peak

    # 3. Interpolate onto the evaluation grid
    y_grid = np.interp(t_eval, t_decay, y_decay)

    # 4. Fit triexponential
    return fit_triexponential_decay(t_eval, y_grid, n_starts=n_starts, seed=seed, extra_starts=extra_starts, bounds=bounds)


# -----------------------------------------------------------------------------
# 3. Unified Likelihood & Diagnostics (Weighting iv)
# -----------------------------------------------------------------------------

def deact_components_log_likelihood(
    *args,
    deact_record: Deactivation | None = None,
    deact_row: Deactivation | None = None,
    pulse_ms: float = 2.0,
    conc_mm: float = 10.0,
    score_areas: bool = True,
    score_tau_w: bool = False,
    **kwargs,
) -> DeactComponentsLikelihoodResult:
    """Score predicted deactivation components against published literature values and errors.

    Weighting (iv):
    Evaluates an exact Gaussian log-likelihood for predicted components against
    published values and errors in `circuitpharm.literature.DEACTIVATIONS[1]`
    (Barberis 2007 Table 1):
      tau_1 = 2.8 +/- 0.3 ms, area_1 = 0.57 +/- 0.04
      tau_2 = 33.4 +/- 4.6 ms, area_2 = 0.23 +/- 0.02
      tau_3 = 221.35 +/- 14.9 ms, area_3 = 0.20 +/- 0.03
      tau_w = 52.5 +/- 2.9 ms (derived)

    Accepts either:
      - `deact_components_log_likelihood(model, ...)`: computes analytic eigen-components
      - `deact_components_log_likelihood(taus, areas, ...)`: scores explicit component sequences
      - `deact_components_log_likelihood(fit_result, ...)`: scores fitted component results
    """
    if len(args) == 2:
        taus_in, areas_in = args[0], args[1]
    elif len(args) == 1:
        arg = args[0]
        if hasattr(arg, "q_matrix") or isinstance(arg, ReceptorModel):
            taus_in, areas_in = predict_components_eigen(
                arg, pulse_ms=pulse_ms, conc_mm=conc_mm, n_components=3
            )
        elif isinstance(arg, (tuple, list)) and len(arg) == 2:
            taus_in, areas_in = arg[0], arg[1]
        else:
            raise TypeError(
                f"Expected ReceptorModel, (taus, areas) tuple, or 2 positional arguments; got {type(arg).__name__}"
            )
    else:
        raise TypeError(
            f"deact_components_log_likelihood takes 1 or 2 positional arguments, got {len(args)}"
        )

    record = deact_record if deact_record is not None else (deact_row if deact_row is not None else DEACTIVATIONS[1])

    taus_pred_arr = np.asarray(taus_in, dtype=float)
    areas_pred_arr = np.asarray(areas_in, dtype=float)

    # Published values & errors
    target_taus, target_areas = record.components
    pub_taus = np.asarray(target_taus, dtype=float)
    pub_areas = np.asarray(target_areas, dtype=float)

    err_taus, err_areas = get_barberis_errors(record)
    err_taus_arr = np.asarray(err_taus, dtype=float)
    err_areas_arr = np.asarray(err_areas, dtype=float)

    # Align number of components to target
    n_comp = len(pub_taus)
    if len(taus_pred_arr) != n_comp or len(areas_pred_arr) != n_comp:
        raise ValueError(
            f"Expected {n_comp} taus and areas, got {len(taus_pred_arr)} taus and {len(areas_pred_arr)} areas."
        )

    z_tau = (taus_pred_arr - pub_taus) / err_taus_arr
    ll_tau_points = -0.5 * (z_tau ** 2 + np.log(2.0 * np.pi * err_taus_arr ** 2))
    ll_taus = float(np.sum(ll_tau_points))

    z_area = (areas_pred_arr - pub_areas) / err_areas_arr
    ll_area_points = -0.5 * (z_area ** 2 + np.log(2.0 * np.pi * err_areas_arr ** 2))
    ll_areas = float(np.sum(ll_area_points)) if score_areas else 0.0

    total_ll = ll_taus + ll_areas

    sum_areas = float(np.sum(areas_pred_arr))
    tau_w_pred = float(np.sum(areas_pred_arr * taus_pred_arr) / sum_areas) if sum_areas != 0 else float(np.sum(areas_pred_arr * taus_pred_arr))
    tau_w_pub = float(record.tau_w_ms)
    tau_w_err = float(record.err_ms)
    z_tw = float((tau_w_pred - tau_w_pub) / tau_w_err)

    if score_tau_w:
        ll_tau_w = -0.5 * (z_tw ** 2 + np.log(2.0 * np.pi * tau_w_err ** 2))
        total_ll += float(ll_tau_w)

    return DeactComponentsLikelihoodResult(
        total_log_likelihood=total_ll,
        ll_taus=ll_taus,
        ll_areas=ll_areas,
        z_taus=tuple(float(x) for x in z_tau),
        z_areas=tuple(float(x) for x in z_area),
        taus_pred=tuple(float(x) for x in taus_pred_arr),
        taus_pub=tuple(float(x) for x in pub_taus),
        taus_err=tuple(float(x) for x in err_taus_arr),
        areas_pred=tuple(float(x) for x in areas_pred_arr),
        areas_pub=tuple(float(x) for x in pub_areas),
        areas_err=tuple(float(x) for x in err_areas_arr),
        tau_w_pred=tau_w_pred,
        tau_w_pub=tau_w_pub,
        tau_w_err=tau_w_err,
        z_tau_w=z_tw,
    )


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

    tau_errs, area_errs = get_barberis_errors(deact_record)
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

    sum_a = float(np.sum(pred_areas))
    tau_w_pred = float(np.sum(pred_areas * pred_taus) / sum_a) if sum_a != 0 else float(np.sum(pred_areas * pred_taus))
    tau_w_sig = (tau_w_pred - float(deact_record.tau_w_ms)) / float(deact_record.err_ms)

    max_dev = verify_eigen_reconstruction(model, pulse_ms=pulse_ms, conc_mm=conc_mm, atol=1e-5)

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
