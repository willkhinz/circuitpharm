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
        warm_start_won: bool = True,
        cost_warm: float = 0.0,
        cost_random: float = 0.0,
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
        instance._warm_start_won = warm_start_won
        instance._cost_warm = cost_warm
        instance._cost_random = cost_random
        return instance

    @property
    def warm_start_won(self) -> bool:
        """Whether the chosen optimum came from a warm start rather than random sampling."""
        return self._warm_start_won

    @property
    def cost_warm(self) -> float:
        """Best least-squares residual cost achieved across warm starts."""
        return self._cost_warm

    @property
    def cost_random(self) -> float:
        """Best least-squares residual cost achieved across random starts."""
        return self._cost_random

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
        chi2_taus: float = 0.0,
        chi2_areas: float = 0.0,
        chi2_total: float = 0.0,
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
        instance.chi2_taus = float(chi2_taus)
        instance.chi2_areas = float(chi2_areas)
        instance.chi2_total = float(chi2_total)
        return instance

    def __repr__(self) -> str:
        return (
            f"DeactComponentsLikelihoodResult(total_ll={self.total_log_likelihood:.4f}, "
            f"ll_taus={self.ll_taus:.4f}, ll_areas={self.ll_areas:.4f}, "
            f"chi2_total={self.chi2_total:.2f}, "
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
    chi2_taus: float = 0.0
    chi2_areas: float = 0.0
    chi2_total: float = 0.0


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


def get_barberis_amplitude_covariance(
    record: Deactivation | None = None,
    *,
    err_areas: Sequence[float] | None = None,
) -> tuple[np.ndarray, np.ndarray, float]:
    """Compute rank-2 covariance matrix, pseudo-inverse, and sub-determinant for fractional amplitudes.

    Brief 10:
    The published fractional amplitudes sum to 1.0 by construction, so their variation
    is strictly constrained to the 2D simplex hyperplane (sum(a_i) = 1). The covariance
    matrix Sigma must satisfy Sigma @ 1 = 0 (zero variance on the sum constraint).

    Given marginal standard errors s = (s_1, s_2, s_3), the diagonal elements are
    Sigma_ii = s_i^2. The condition that row sums equal zero uniquely determines the
    off-diagonal covariances:
        Sigma_ij + Sigma_ik = -s_i^2  for distinct i, j, k
    with unique algebraic solution:
        Sigma_ij = 0.5 * (s_k^2 - s_i^2 - s_j^2).

    Parameters
    ----------
    record : Deactivation, optional
        Deactivation record (defaults to DEACTIVATIONS[1], Barberis et al. 2007).
    err_areas : Sequence[float], optional
        Explicit component area standard errors. If None, extracted from `record`.

    Returns
    -------
    cov : np.ndarray, shape (3, 3)
        Singular rank-2 covariance matrix satisfying cov @ 1 = 0 and diag(cov) = s^2.
    cov_pinv : np.ndarray, shape (3, 3)
        Moore-Penrose pseudo-inverse of cov, satisfying cov_pinv @ 1 = 0.
    det_sub : float
        Determinant of any 2x2 principal submatrix Sigma_{(-k, -k)}.
        Equal to det*(Sigma) / 3, identically across all k in {0, 1, 2}.
    """
    if err_areas is None:
        if record is None:
            record = DEACTIVATIONS[1]
        _, err_areas = get_barberis_errors(record)

    s = np.asarray(err_areas, dtype=float)
    if s.ndim != 1 or len(s) != 3:
        raise ValueError(
            f"Fractional amplitude covariance requires exactly 3 components, got {len(s)}."
        )

    v = s ** 2
    cov = np.zeros((3, 3), dtype=float)
    cov[0, 0] = v[0]
    cov[1, 1] = v[1]
    cov[2, 2] = v[2]
    cov[0, 1] = cov[1, 0] = 0.5 * (v[2] - v[0] - v[1])
    cov[0, 2] = cov[2, 0] = 0.5 * (v[1] - v[0] - v[2])
    cov[1, 2] = cov[2, 1] = 0.5 * (v[0] - v[1] - v[2])

    cov_pinv = np.linalg.pinv(cov, hermitian=True)
    det_sub = float(cov[0, 0] * cov[1, 1] - cov[0, 1] ** 2)

    return cov, cov_pinv, det_sub


def compute_amplitude_chi2(
    areas_pred: Sequence[float],
    pub_areas: Sequence[float],
    err_areas: Sequence[float],
    *,
    method: str = "full_covariance",
    dropped_index: int | None = None,
) -> float:
    """Compute chi^2 for constrained fractional amplitudes.

    Parameters
    ----------
    areas_pred : Sequence[float]
        Predicted fractional amplitudes (must sum to ~1.0).
    pub_areas : Sequence[float]
        Published fractional amplitudes (must sum to ~1.0).
    err_areas : Sequence[float]
        Published marginal standard errors (SEMs).
    method : {"full_covariance", "drop_one", "drop_one_diagonal", "legacy_unconstrained"}
        - "full_covariance" (default): Quadratic form d^T Sigma^+ d using the Moore-Penrose
          pseudo-inverse of the rank-2 degenerate covariance matrix (Brief 10).
        - "drop_one": Explicitly drops component `dropped_index` (default: 2) and evaluates
          d_{-k}^T (Sigma_{-k, -k})^{-1} d_{-k}, which is mathematically identical to full_covariance.
        - "drop_one_diagonal": Drops component `dropped_index` and scores remaining components
          with naive diagonal variances (demonstrating non-invariance).
        - "legacy_unconstrained": Naive unconstrained sum of 3 independent z-scores (Brief 9),
          which has an arithmetic floor of 452.99 when a_1 = 0.
    dropped_index : int, optional
        Component index to drop for "drop_one" or "drop_one_diagonal" (0, 1, or 2).

    Returns
    -------
    chi2 : float
        Evaluated chi^2 metric.
    """
    a_pred = np.asarray(areas_pred, dtype=float)
    a_pub = np.asarray(pub_areas, dtype=float)
    s_pub = np.asarray(err_areas, dtype=float)
    d = a_pred - a_pub

    if method == "full_covariance":
        _, cov_pinv, _ = get_barberis_amplitude_covariance(err_areas=s_pub)
        return float(d @ cov_pinv @ d)
    elif method == "drop_one":
        k = 2 if dropped_index is None else int(dropped_index)
        keep = [i for i in range(3) if i != k]
        cov, _, _ = get_barberis_amplitude_covariance(err_areas=s_pub)
        sub_cov = cov[np.ix_(keep, keep)]
        sub_inv = np.linalg.inv(sub_cov)
        d_sub = d[keep]
        return float(d_sub @ sub_inv @ d_sub)
    elif method == "drop_one_diagonal":
        k = 2 if dropped_index is None else int(dropped_index)
        keep = [i for i in range(3) if i != k]
        z_sub = d[keep] / s_pub[keep]
        return float(np.sum(z_sub ** 2))
    elif method == "legacy_unconstrained":
        z = d / s_pub
        return float(np.sum(z ** 2))
    else:
        raise ValueError(
            f"Unknown amplitude comparison method: {method!r}. "
            f"Must be one of 'full_covariance', 'drop_one', 'drop_one_diagonal', 'legacy_unconstrained'."
        )


# -----------------------------------------------------------------------------
# 1. Analytic Eigen-Decomposition Predictor (Brief 7A)
# -----------------------------------------------------------------------------

def predict_components_eigen_all(
    model: ReceptorModel,
    *,
    pulse_ms: float = 2.0,
    conc_mm: float = 10.0,
    open_index: int | Sequence[int] | None = None,
    return_raw: bool = False,
) -> tuple[np.ndarray, np.ndarray] | tuple[np.ndarray, np.ndarray, np.ndarray, float]:
    """Compute all non-zero relaxation components (taus, areas) from the generator matrix.

    Parameters
    ----------
    model : ReceptorModel
        Markov kinetic model providing `q_matrix` and `state_distribution`.
    pulse_ms : float
        Pulse duration in ms (default 2.0 ms for Barberis 2007 brief pulse).
    conc_mm : float
        GABA concentration during pulse in mM (default 10.0 mM).
    open_index : int, sequence of ints, or None
        State index or indices for conducting open channel(s). If None, inspects
        model.open_indices or model.open_index, falling back to OPEN_STATE_INDEX (3).

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

    if open_index is None:
        if hasattr(model, "open_indices"):
            open_indices = tuple(int(i) for i in model.open_indices)
        elif hasattr(model, "open_index"):
            idx = model.open_index
            open_indices = tuple(int(i) for i in idx) if isinstance(idx, (tuple, list, np.ndarray)) else (int(idx),)
        else:
            open_indices = (OPEN_STATE_INDEX,)
    elif isinstance(open_index, (int, np.integer)):
        open_indices = (int(open_index),)
    else:
        open_indices = tuple(int(i) for i in open_index)

    # 1. Initial condition: 2 ms pulse of 10 mM GABA starting from equilibrium at 0 agonist
    gaba_pulse_um = float(conc_mm) * 1000.0
    p0 = np.asarray(model.state_distribution(0.0), dtype=float)
    q_pulse = np.asarray(model.q_matrix(gaba_pulse_um), dtype=float)
    p_pulse = p0 @ la.expm(q_pulse * float(pulse_ms))

    # Numerical hygiene: ensure state probabilities remain on simplex
    p_pulse = np.clip(p_pulse, 0.0, 1.0)
    p_pulse = p_pulse / np.sum(p_pulse)
    p_open_onset = float(np.sum(p_pulse[list(open_indices)]))

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

    # 3. Observable: projection onto open state coordinate(s)
    v_nz = np.real(v)
    v_inv = la.inv(v_nz)

    # For multiple open states, project onto sum of basis coordinates
    # a_i = (P(pulse_ms) @ V)_i * sum_{j in open_indices} (V^{-1})_{i, j}
    proj_pulse = p_pulse @ v_nz
    summed_v_inv = np.sum(v_inv[non_zero, :][:, list(open_indices)], axis=1)
    amps = proj_pulse[non_zero] * summed_v_inv

    # Check: sum of amplitudes must reproduce peak open probability at onset
    sum_amps = float(np.sum(amps))
    if not np.isclose(sum_amps, p_open_onset, atol=1e-8, rtol=1e-6):
        raise ValueError(
            f"Sum of eigen-amplitudes ({sum_amps:.8e}) does not reproduce open probability "
            f"at pulse end ({p_open_onset:.8e})."
        )

    # Reject unphysical or ill-conditioned near-zero open probability / amplitude sum (Brief 9B §3)
    min_amplitude_sum = 1e-4
    if p_open_onset < min_amplitude_sum or abs(sum_amps) < min_amplitude_sum:
        raise ValueError(
            f"Deactivation amplitude sum is near zero or ill-conditioned "
            f"(P_o(0) = {p_open_onset:.4e}, sum(amps) = {sum_amps:.4e} < {min_amplitude_sum:.1e}); "
            f"relaxation waveform carries negligible open probability."
        )

    # Normalise areas to peak (p_open at relaxation onset)
    areas = amps / p_open_onset

    # Sort fast-to-slow
    order = np.argsort(taus)
    if return_raw:
        return taus[order], areas[order], amps[order], p_open_onset
    return taus[order], areas[order]


def predict_components_eigen(
    model: ReceptorModel,
    *,
    pulse_ms: float = 2.0,
    conc_mm: float = 10.0,
    n_components: int = 3,
    open_index: int | Sequence[int] | None = None,
    return_all: bool = False,
    return_raw: bool = False,
    min_relative_ratio: float = 0.10,
) -> tuple[np.ndarray, np.ndarray] | tuple[np.ndarray, np.ndarray, np.ndarray, float]:
    """Predict deactivation (taus, areas) from generator matrix eigenvalues.

    By default returns the 3 physically active relaxation components, resolving the 4-vs-3 problem
    by dropping the unbinding mode AR -> R which carries mathematically zero amplitude in A2O.

    Rejects parameter vectors where the raw amplitude sum of the selected components is near
    zero or small relative to P_o(0) (threshold min_relative_ratio, default 0.10; Brief 9B §3).

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
    open_index : int, sequence of ints, or None
        Index or indices of the conducting open state(s). Defaults to model open state(s).
    return_all : bool
        If True, returns all non-zero components without reduction.
    return_raw : bool
        If True, returns (taus, areas, raw_amplitudes, p_open_onset).
    min_relative_ratio : float
        Minimum required ratio |sum(raw_chosen)| / P_o(0) (default 0.10).

    Returns
    -------
    taus : np.ndarray
        Component time constants in ms, sorted ascending.
    areas : np.ndarray
        Component normalised areas (fractions).
    """
    taus_all, areas_all, amps_all, p_open_onset = predict_components_eigen_all(
        model, pulse_ms=pulse_ms, conc_mm=conc_mm, open_index=open_index, return_raw=True
    )
    if return_all or n_components >= len(taus_all):
        if return_raw:
            return taus_all, areas_all, amps_all, p_open_onset
        return taus_all, areas_all

    top_indices = np.argsort(-np.abs(areas_all))[:n_components]
    sort_taus = np.argsort(taus_all[top_indices])
    chosen = top_indices[sort_taus]

    chosen_taus = taus_all[chosen]
    chosen_areas = areas_all[chosen]
    chosen_amps = amps_all[chosen]

    # Conditioning check: raw amplitude sum of selected components must not cancel to near zero (Brief 9B §3)
    raw_sum_chosen = float(np.sum(chosen_amps))
    if p_open_onset > 0.0 and (abs(raw_sum_chosen) / p_open_onset) < min_relative_ratio:
        raise ValueError(
            f"Ill-conditioned amplitudes: raw amplitude sum of top {n_components} components "
            f"({raw_sum_chosen:.4e}) is small relative to P_o(0) ({p_open_onset:.4e}, "
            f"ratio {abs(raw_sum_chosen)/p_open_onset:.4f} < {min_relative_ratio}); vector rejected."
        )

    if return_raw:
        return chosen_taus, chosen_areas, chosen_amps, p_open_onset
    return chosen_taus, chosen_areas


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
    open_index: int | Sequence[int] | None = None,
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
    open_index : int, sequence of ints, or None
        Index or indices of the conducting open state(s). Defaults to model open state(s).
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

    # Resolve open indices
    if open_index is None:
        if hasattr(model, "open_indices"):
            open_indices = tuple(int(i) for i in model.open_indices)
        elif hasattr(model, "open_index"):
            idx = model.open_index
            open_indices = tuple(int(i) for i in idx) if isinstance(idx, (tuple, list, np.ndarray)) else (int(idx),)
        else:
            open_indices = (OPEN_STATE_INDEX,)
    elif isinstance(open_index, (int, np.integer)):
        open_indices = (int(open_index),)
    else:
        open_indices = tuple(int(i) for i in open_index)

    # 2. All components and unnormalized amplitudes
    taus_all, areas_all = predict_components_eigen_all(
        model, pulse_ms=pulse_ms, conc_mm=conc_mm, open_index=open_indices
    )
    p_open_onset = float(np.sum(p_pulse[list(open_indices)]))
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
    starts: Sequence[np.ndarray] | None = None,
    extra_starts: Sequence[np.ndarray] | None = None,
    bounds: tuple[np.ndarray, np.ndarray] | None = None,
) -> DeactFitResult:
    """Fit a normalized triexponential decay sum(a_i * exp(-t / tau_i)) to a time series.

    The model parameterisation enforces sum(a_i) = 1.0 (with a_3 = 1 - a_1 - a_2) and
    searches in log10(tau) to cover timescales spanning from 0.01 ms to 10,000 ms.

    Warm-start and multi-start discipline (Brief 9C §2):
    - Accepts explicit starting parameter vectors via `starts` (e.g. from eigen decomposition).
    - Retains random Latin Hypercube starts (`n_starts`).
    - Compares warm starts vs random starts: asserts warm-started optimum is no worse than
      the best random-start optimum on least-squares cost. If random beats warm, random is
      reported as the fit, preventing warm-start bias.
    - Tracks top optima, condition number, rank deficiency, and collapsed taus.
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

    # 1. Warm starts (explicit or default literature/heuristic)
    warm_candidates = []
    if starts is not None:
        for s in starts:
            warm_candidates.append(np.asarray(s, dtype=float))
    else:
        warm_candidates = [
            np.array([np.log10(2.8), np.log10(33.4), np.log10(221.35), 0.57, 0.23]),
            np.array([np.log10(5.0), np.log10(50.0), np.log10(300.0), 0.33, 0.33]),
            np.array([np.log10(1.0), np.log10(15.0), np.log10(500.0), 0.40, 0.30]),
        ]
    if extra_starts:
        for s in extra_starts:
            warm_candidates.append(np.asarray(s, dtype=float))

    # 2. Random Latin Hypercube starts
    sampled = []
    if n_starts > 0:
        sampler = qmc.LatinHypercube(d=5, seed=seed)
        sampled_raw = lo_b + sampler.random(n_starts) * (hi_b - lo_b)
        for row in sampled_raw:
            row[:3] = np.sort(row[:3])
            if row[3] >= 0 and row[4] >= 0 and row[3] + row[4] > 0.95:
                scale = 0.85 / (row[3] + row[4])
                row[3] *= scale
                row[4] *= scale
            sampled.append(row)

    tagged_starts = [(s, "warm") for s in warm_candidates] + [(s, "random") for s in sampled]

    optima = []
    for s0, origin in tagged_starts:
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
            "origin": origin,
        })

    if not optima:
        raise RuntimeError("None of the triexponential starts converged.")

    # Guard: Compare best warm start against best random start (Brief 9C §2)
    warm_opts = [o for o in optima if o["origin"] == "warm"]
    rand_opts = [o for o in optima if o["origin"] == "random"]
    best_warm = min(warm_opts, key=lambda x: x["cost"]) if warm_opts else None
    best_rand = min(rand_opts, key=lambda x: x["cost"]) if rand_opts else None

    cost_warm = best_warm["cost"] if best_warm else float("inf")
    cost_rand = best_rand["cost"] if best_rand else float("inf")

    if cost_rand < cost_warm - 1e-9:
        warm_start_won = False
        best = best_rand
    else:
        warm_start_won = True
        best = best_warm if best_warm else best_rand

    top = sorted(optima, key=lambda item: item["cost"])[: min(5, len(optima))]
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
        warm_start_won=warm_start_won,
        cost_warm=cost_warm,
        cost_random=cost_rand,
    )


def predict_components_fit(
    model: Any,
    *,
    pulse_ms: float = 2.0,
    conc_mm: float = 10.0,
    fit_window_ms: float = _BARBERIS_FIT_WINDOW_MS,
    grid: np.ndarray | None = None,
    n_starts: int = 4,
    seed: int = 42,
    uniform_grid: bool = False,
    n_pulse: int = 200,
    n_tail: int = 500,
    starts: Sequence[np.ndarray] | None = None,
    extra_starts: Sequence[np.ndarray] | None = None,
    bounds: tuple[np.ndarray, np.ndarray] | None = None,
    warm_start_from_eigen: bool = True,
) -> DeactFitResult:
    """Predict triexponential deactivation components (taus, areas) by fitting simulated waveform.

    Implements the experimental procedure of Barberis et al. 2007:
    1. Protocol simulation: 2 ms pulse of 10 mM GABA, followed by washout to 0 mM over fit_window_ms
       (default 1000.0 ms), computed via model.simulate_waveform.
    2. Peak extraction & normalisation: current decay is re-zeroed at peak open probability and
       divided by peak open probability so y(0) = 1.0, matching BARBERIS2007_DEACT.normalized_current.
    3. Grid sampling: evaluated on the dataset's hybrid log/linear grid (_T_BARBERIS) by default,
       or on a uniform grid if uniform_grid=True (sensitivity check).
    4. Triexponential fit: warm-started from the generator's dominant eigen components (Brief 9C §2)
       alongside random Latin Hypercube starts, with rank-deficiency / tie detection.
    """
    if grid is not None:
        t_eval = np.asarray(grid, dtype=float)
    elif uniform_grid:
        t_eval = np.linspace(0.0, float(fit_window_ms), len(_T_BARBERIS))
    else:
        t_eval = _T_BARBERIS

    if isinstance(model, np.ndarray):
        y_eval = np.asarray(model, dtype=float)
        return fit_triexponential_decay(t_eval, y_eval, n_starts=n_starts, seed=seed, starts=starts, extra_starts=extra_starts, bounds=bounds)

    if isinstance(model, DeactivationDataset) or hasattr(model, "normalized_current"):
        ds_t = getattr(model, "time_ms", t_eval)
        ds_y = getattr(model, "normalized_current")
        if np.array_equal(ds_t, t_eval):
            y_eval = np.asarray(ds_y, dtype=float)
        else:
            y_eval = np.interp(t_eval, ds_t, ds_y)
        return fit_triexponential_decay(t_eval, y_eval, n_starts=n_starts, seed=seed, starts=starts, extra_starts=extra_starts, bounds=bounds)

    if not hasattr(model, "q_matrix"):
        from ..models.base import ObservableMismatch
        raise ObservableMismatch(
            f"predict_components_fit requires a Markov model with state transitions (`q_matrix`), "
            f"got {type(model).__name__}."
        )

    # 1. Warm-start from eigen decomposition if enabled (Brief 9C §2)
    fit_starts = list(starts) if starts is not None else []
    if not fit_starts and warm_start_from_eigen:
        try:
            e_taus, e_areas = predict_components_eigen(
                model, pulse_ms=pulse_ms, conc_mm=conc_mm, n_components=3
            )
            lo_b = np.array([-2.0, -1.0, 0.5, 0.0, 0.0]) if bounds is None else bounds[0]
            hi_b = np.array([1.5, 2.8, 4.2, 1.0, 1.0]) if bounds is None else bounds[1]
            p_eigen = np.array([
                np.log10(max(float(e_taus[0]), 10.0 ** lo_b[0])),
                np.log10(max(float(e_taus[1]), 10.0 ** lo_b[1])),
                np.log10(max(float(e_taus[2]), 10.0 ** lo_b[2])),
                float(e_areas[0]),
                float(e_areas[1]),
            ])
            p_eigen = np.clip(p_eigen, lo_b + 1e-4, hi_b - 1e-4)
            fit_starts.append(p_eigen)
        except Exception:
            pass

    # 2. Simulate the protocol
    t_sim = np.unique(np.concatenate([
        np.linspace(0.0, float(pulse_ms), n_pulse),
        float(pulse_ms) + np.geomspace(1e-4, float(fit_window_ms), n_tail),
    ]))
    gaba_um = np.where(t_sim <= float(pulse_ms), float(conc_mm) * 1000.0, 0.0)

    res = model.simulate_waveform(t_sim, gaba_um)
    p_open = np.asarray(res.p_open, dtype=float)

    # 3. Extract and normalise decay
    i_peak = int(np.argmax(p_open))
    peak = float(p_open[i_peak])
    if not (np.isfinite(peak) and peak > 0.0):
        raise ValueError(
            f"Model peak open probability is non-positive or non-finite ({peak!r}); decay cannot be normalised."
        )

    t_decay = t_sim[i_peak:] - t_sim[i_peak]
    y_decay = p_open[i_peak:] / peak

    # 4. Interpolate onto the evaluation grid
    y_grid = np.interp(t_eval, t_decay, y_decay)

    # 5. Fit triexponential
    return fit_triexponential_decay(
        t_eval, y_grid, n_starts=n_starts, seed=seed,
        starts=fit_starts if fit_starts else None,
        extra_starts=extra_starts, bounds=bounds,
    )


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

    # Reject near-zero amplitude sums (Brief 9B §3)
    sum_areas = float(np.sum(areas_pred_arr))
    if abs(sum_areas) < 1e-4:
        raise ValueError(
            f"Ill-conditioned amplitudes: sum of predicted component areas ({sum_areas:.4e}) "
            f"is near zero (< 1e-4); vector rejected."
        )

    z_tau = (taus_pred_arr - pub_taus) / err_taus_arr
    ll_tau_points = -0.5 * (z_tau ** 2 + np.log(2.0 * np.pi * err_taus_arr ** 2))
    ll_taus = float(np.sum(ll_tau_points))
    chi2_taus = float(np.sum(z_tau ** 2))

    # Marginal z-scores for reporting
    z_area = (areas_pred_arr - pub_areas) / err_areas_arr

    if kwargs.get("legacy_unconstrained", False):
        chi2_areas = float(np.sum(z_area ** 2))
        ll_area_points = -0.5 * (z_area ** 2 + np.log(2.0 * np.pi * err_areas_arr ** 2))
        ll_areas = float(np.sum(ll_area_points)) if score_areas else 0.0
    else:
        cov, cov_pinv, det_sub = get_barberis_amplitude_covariance(record)
        d_area = areas_pred_arr - pub_areas
        chi2_areas = float(d_area @ cov_pinv @ d_area)
        # Rank-2 degenerate Gaussian density on the simplex hyperplane (sum(a_i) = 1)
        # Normalization factor for 2 DOF is (2*pi)^2 * det(Sigma_sub)
        log_norm_areas = float(np.log((2.0 * np.pi) ** 2 * det_sub))
        ll_areas = float(-0.5 * (chi2_areas + log_norm_areas)) if score_areas else 0.0

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
        chi2_taus=chi2_taus,
        chi2_areas=chi2_areas,
        chi2_total=chi2_taus + chi2_areas,
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

    chi2_taus = float(np.sum((tau_diff / tau_errs_arr) ** 2))
    ll_tau = -0.5 * float(chi2_taus + np.sum(np.log(2.0 * np.pi * (tau_errs_arr ** 2))))

    cov, cov_pinv, det_sub = get_barberis_amplitude_covariance(deact_record)
    chi2_areas = float(area_diff @ cov_pinv @ area_diff)
    log_norm_areas = float(np.log((2.0 * np.pi) ** 2 * det_sub))
    ll_area = float(-0.5 * (chi2_areas + log_norm_areas))

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
        chi2_taus=chi2_taus,
        chi2_areas=chi2_areas,
        chi2_total=chi2_taus + chi2_areas,
    )
