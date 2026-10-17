"""Predicted macroscopic deactivation components by fitting simulated waveforms.

Agent brief 7B (Brief 7, weighting iv):
Barberis et al. 2007 (Eur J Neurosci 25(9):2726-2740, PMC1950087, Table 1) reported
triexponential deactivation parameters (taus and fractional areas) fitted to current
decay following a 2 ms pulse of 10 mM GABA.

This module provides:
1. `predict_components_fit`: Pure function simulating the declared protocol using the
   existing waveform machinery (`simulate_waveform`), normalising the decay as
   `BARBERIS2007_DEACT.normalized_current` is normalised, and fitting a triexponential
   to extract predicted (taus, areas) ordered fastest-first.
2. `fit_triexponential_decay`: Multi-start bounded non-linear least squares solver
   reporting top-optima spread, ties, Jacobian condition number, and rank deficiency.
3. `deact_components_log_likelihood`: Likelihood term scoring predicted components
   against `circuitpharm.literature.DEACTIVATIONS[1]` published values and errors.
"""
from __future__ import annotations

from dataclasses import dataclass
import re
from typing import Any, Sequence

import numpy as np
from scipy.optimize import least_squares
from scipy.stats import qmc

from ..literature import DEACTIVATIONS, Deactivation
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


@dataclass(frozen=True)
class DeactComponentsLikelihoodResult:
    """Log-likelihood evaluation of predicted components against literature."""

    total_log_likelihood: float
    ll_taus: float
    ll_areas: float
    z_taus: tuple[float, float, float]
    z_areas: tuple[float, float, float]
    taus_pred: tuple[float, float, float]
    taus_pub: tuple[float, float, float]
    taus_err: tuple[float, float, float]
    areas_pred: tuple[float, float, float]
    areas_pub: tuple[float, float, float]
    areas_err: tuple[float, float, float]
    tau_w_pred: float
    tau_w_pub: float
    tau_w_err: float
    z_tau_w: float


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
        raise ValueError(f"Could not parse component errors from note: {note!r}")
    tau_errs = (float(m.group(1)), float(m.group(2)), float(m.group(3)))
    area_errs = (float(m.group(4)), float(m.group(5)), float(m.group(6)))
    return tau_errs, area_errs


def get_barberis_errors(
    deact_row: Deactivation | None = None,
) -> tuple[tuple[float, float, float], tuple[float, float, float]]:
    """Retrieve published component errors from structured field or prose fallback."""
    if deact_row is None:
        deact_row = DEACTIVATIONS[1]
    if getattr(deact_row, "component_errors", ()):
        t_errs, a_errs = deact_row.component_errors
        return tuple(float(x) for x in t_errs), tuple(float(x) for x in a_errs)
    return parse_barberis_component_errors(deact_row.note)


def fit_triexponential_decay(
    t_ms: np.ndarray,
    y: np.ndarray,
    *,
    n_starts: int = 16,
    seed: int = 42,
    spread_tol: float = 0.05,
    extra_starts: Sequence[np.ndarray] | None = None,
) -> DeactFitResult:
    """Fit a normalized triexponential decay Σ aᵢ exp(−t/τᵢ) to a time series.

    The model parameterisation enforces Σ aᵢ = 1.0 (with a₃ = 1 − a₁ − a₂) and
    searches in log10(τ) to cover timescales spanning from 0.01 ms to 10,000 ms.

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

    # Parameter bounds:
    # log10_tau: [-2.0, 4.0] (0.01 ms to 10,000 ms)
    # areas: a1, a2 in [0.0, 1.0]
    lo_b = np.array([-2.0, -1.0, 0.5, 0.0, 0.0])
    hi_b = np.array([1.5, 2.8, 4.2, 1.0, 1.0])
    bounds = (lo_b, hi_b)

    # Multi-start generation
    sampler = qmc.LatinHypercube(d=5, seed=seed)
    sampled = lo_b + sampler.random(n_starts) * (hi_b - lo_b)
    # Ensure ordered tau starts
    for row in sampled:
        row[:3] = np.sort(row[:3])
        # ensure a1 + a2 <= 1.0
        if row[3] + row[4] > 0.95:
            scale = 0.85 / (row[3] + row[4])
            row[3] *= scale
            row[4] *= scale

    # Warm starts: published Barberis parameters, plus equal-area geometry
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

        # Order components fastest-first
        tau_raw = np.array([t1, t2, t3])
        area_raw = np.array([a1, a2, a3])
        perm = np.argsort(tau_raw)
        taus_ordered = tuple(float(tau_raw[i]) for i in perm)
        areas_ordered = tuple(float(area_raw[i]) for i in perm)

        # Singular values of Jacobian at solution
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

    # Sort by residual cost
    optima.sort(key=lambda item: item["cost"])
    best = optima[0]

    # Evaluate spread among top optima (up to 5)
    top = optima[: min(5, len(optima))]
    tau_mat = np.array([item["taus"] for item in top])
    area_mat = np.array([item["areas"] for item in top])

    spread_t = tuple(float(np.ptp(tau_mat[:, i])) for i in range(3))
    spread_a = tuple(float(np.ptp(area_mat[:, i])) for i in range(3))

    # Log10 spread for taus to assess convergence in decades
    log_tau_mat = np.log10(np.maximum(tau_mat, 1e-12))
    log_spread_t = np.ptp(log_tau_mat, axis=0)
    is_unique = bool(np.all(log_spread_t <= spread_tol) and np.all(np.array(spread_a) <= spread_tol))

    # Conditioning analysis:
    # 1. Do any two time constants merge within 5%?
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

    # Compute R-squared
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
    # Grid selection
    if grid is not None:
        t_eval = np.asarray(grid, dtype=float)
    elif uniform_grid:
        t_eval = np.linspace(0.0, float(fit_window_ms), len(_T_BARBERIS))
    else:
        t_eval = _T_BARBERIS

    # Support passing an array or dataset directly (e.g. for ground-truth verification)
    if isinstance(model, np.ndarray):
        y_eval = np.asarray(model, dtype=float)
        return fit_triexponential_decay(t_eval, y_eval, n_starts=n_starts, seed=seed)

    if isinstance(model, DeactivationDataset) or hasattr(model, "normalized_current"):
        ds_t = getattr(model, "time_ms", t_eval)
        ds_y = getattr(model, "normalized_current")
        # If user requested a different evaluation grid, interpolate
        if np.array_equal(ds_t, t_eval):
            y_eval = np.asarray(ds_y, dtype=float)
        else:
            y_eval = np.interp(t_eval, ds_t, ds_y)
        return fit_triexponential_decay(t_eval, y_eval, n_starts=n_starts, seed=seed)

    # Validate that model has state transitions (Markov generator)
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

    # Waveform simulation via model's existing machinery
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
    return fit_triexponential_decay(t_eval, y_grid, n_starts=n_starts, seed=seed)


def deact_components_log_likelihood(
    taus: Sequence[float],
    areas: Sequence[float],
    *,
    deact_row: Deactivation | None = None,
) -> DeactComponentsLikelihoodResult:
    """Score predicted deactivation components against published literature values and errors.

    Weighting (iv):
    Evaluates an exact Gaussian log-likelihood for each predicted component against
    the published values and errors in `circuitpharm.literature.DEACTIVATIONS[1]`
    (Barberis 2007 Table 1):
      τ₁ = 2.8 ± 0.3 ms, area₁ = 0.57 ± 0.04
      τ₂ = 33.4 ± 4.6 ms, area₂ = 0.23 ± 0.02
      τ₃ = 221.35 ± 14.9 ms, area₃ = 0.20 ± 0.03
      τ_w = 52.5 ± 2.9 ms (derived)

    Errors are retrieved directly from `deact_row.component_errors` or parsed from
    `deact_row.note`; never retyped as literals.
    """
    if deact_row is None:
        deact_row = DEACTIVATIONS[1]

    taus_pred_arr = np.asarray(taus, dtype=float)
    areas_pred_arr = np.asarray(areas, dtype=float)
    if taus_pred_arr.shape != (3,) or areas_pred_arr.shape != (3,):
        raise ValueError(f"Expected 3 taus and 3 areas, got {taus_pred_arr.shape} and {areas_pred_arr.shape}")

    # Read published values from literature record
    pub_taus_tup, pub_areas_tup = deact_row.components
    pub_taus = np.asarray(pub_taus_tup, dtype=float)
    pub_areas = np.asarray(pub_areas_tup, dtype=float)

    # Read published errors from structured field or note prose
    err_taus_tup, err_areas_tup = get_barberis_errors(deact_row)
    err_taus = np.asarray(err_taus_tup, dtype=float)
    err_areas = np.asarray(err_areas_tup, dtype=float)

    # Gaussian log-likelihood on taus
    z_tau = (taus_pred_arr - pub_taus) / err_taus
    ll_taus_per_point = -0.5 * (z_tau ** 2 + np.log(2.0 * np.pi * err_taus ** 2))
    ll_taus = float(np.sum(ll_taus_per_point))

    # Gaussian log-likelihood on areas
    z_area = (areas_pred_arr - pub_areas) / err_areas
    ll_areas_per_point = -0.5 * (z_area ** 2 + np.log(2.0 * np.pi * err_areas ** 2))
    ll_areas = float(np.sum(ll_areas_per_point))

    total_ll = ll_taus + ll_areas

    # Derived weighted tau
    tau_w_pred = float(np.sum(areas_pred_arr * taus_pred_arr))
    tau_w_pub = float(deact_row.tau_w_ms)
    tau_w_err = float(deact_row.err_ms)
    z_tw = float((tau_w_pred - tau_w_pub) / tau_w_err)

    return DeactComponentsLikelihoodResult(
        total_log_likelihood=total_ll,
        ll_taus=ll_taus,
        ll_areas=ll_areas,
        z_taus=tuple(float(x) for x in z_tau),
        z_areas=tuple(float(x) for x in z_area),
        taus_pred=tuple(float(x) for x in taus_pred_arr),
        taus_pub=pub_taus_tup,
        taus_err=err_taus_tup,
        areas_pred=tuple(float(x) for x in areas_pred_arr),
        areas_pub=pub_areas_tup,
        areas_err=err_areas_tup,
        tau_w_pred=tau_w_pred,
        tau_w_pub=tau_w_pub,
        tau_w_err=tau_w_err,
        z_tau_w=z_tw,
    )
