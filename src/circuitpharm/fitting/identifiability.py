"""Profile Likelihood and Fisher Information Matrix (FIM) identifiability analysis."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Sequence
import numpy as np
from scipy.optimize import minimize

from ..models.kinetic_jw95 import KineticAllosteryModel
from .data import DOSE_RESPONSE_BENCHMARK


@dataclass(frozen=True)
class ProfileLikelihoodResult:
    """Profile likelihood curve and identifiability classification for a parameter."""
    param_name: str
    optimal_value: float
    grid_values: np.ndarray
    profile_costs: np.ndarray
    ci_95_bounds: tuple[float, float]
    is_identifiable: bool
    diagnostic_reason: str


def compute_kinetic_objective(
    params_vec: np.ndarray,
    param_keys: Sequence[str] = ("kon", "koff", "beta", "alpha", "d", "r"),
) -> float:
    """Compute joint negative log-likelihood / weighted chi-squared against empirical benchmarks."""
    p_dict = {k: max(float(v), 1e-6) for k, v in zip(param_keys, params_vec)}
    model = KineticAllosteryModel(**p_dict)

    # 1. Dose-response residuals
    concs = DOSE_RESPONSE_BENCHMARK.concs_um
    obs = DOSE_RESPONSE_BENCHMARK.mean_response
    sem = DOSE_RESPONSE_BENCHMARK.sem
    pred = model.dose_response(concs, pam_factor=1.0)
    chi2_dr = np.sum(((pred - obs) / sem) ** 2)

    # 2. Deactivation constraint: IPSC decay should be in 10-30 ms range
    # Rough analytical decay time constant approximation: tau ~ (1 + beta/alpha) / (koff * (1 + beta/alpha)) -> 15 ms
    # Or exact: IPSC tau ~ 15 ms
    # We add a mild penalty if po_max is far from 0.75
    chi2_pomax = ((model.po_max - 0.75) / 0.05) ** 2

    return float(chi2_dr + chi2_pomax)


def compute_profile_likelihood(
    param_name: str,
    base_model: KineticAllosteryModel,
    scan_factors: Sequence[float] = (0.3, 0.5, 0.8, 1.0, 1.5, 2.0, 3.0),
    threshold_delta_chi2: float = 3.84,  # chi^2(1) at 95% confidence
) -> ProfileLikelihoodResult:
    """Compute profile likelihood curve for a specific kinetic rate parameter."""
    keys = list(base_model.param_names)
    assert param_name in keys, f"Unknown parameter {param_name}"
    target_idx = keys.index(param_name)

    base_val = getattr(base_model, param_name)
    opt_val = float(base_val)

    # Baseline cost
    x0 = np.array([getattr(base_model, k) for k in keys])
    base_cost = compute_kinetic_objective(x0, keys)

    grid_vals = np.array([opt_val * f for f in scan_factors], dtype=float)
    profile_costs = []

    free_indices = [i for i in range(len(keys)) if i != target_idx]

    for val in grid_vals:
        # Constrain param_name = val, optimize free parameters
        def sub_obj(free_vec: np.ndarray) -> float:
            full = np.zeros(len(keys))
            full[target_idx] = val
            for idx, f_val in zip(free_indices, free_vec):
                full[idx] = f_val
            return compute_kinetic_objective(full, keys)

        init_free = x0[free_indices]
        bounds = [(1e-5, 100.0) for _ in free_indices]
        res = minimize(sub_obj, init_free, method="L-BFGS-B", bounds=bounds)
        cost = float(res.fun) if res.success else float(base_cost + 10.0)
        profile_costs.append(cost)

    prof_costs_arr = np.array(profile_costs)
    delta_costs = prof_costs_arr - base_cost

    # Find where delta_costs crosses 3.84
    ci_mask = delta_costs <= threshold_delta_chi2
    if np.all(ci_mask):
        is_id = False
        reason = "Structurally or practically non-identifiable: cost remains below 95% threshold across entire sweep."
        ci_bounds = (float(grid_vals[0]), float(grid_vals[-1]))
    else:
        is_id = True
        reason = "Identifiable: profile cost exceeds 95% confidence threshold within scanned range."
        valid_vals = grid_vals[ci_mask]
        ci_bounds = (float(valid_vals[0]), float(valid_vals[-1]))

    return ProfileLikelihoodResult(
        param_name=param_name,
        optimal_value=opt_val,
        grid_values=grid_vals,
        profile_costs=prof_costs_arr,
        ci_95_bounds=ci_bounds,
        is_identifiable=is_id,
        diagnostic_reason=reason,
    )


def compute_fisher_information_matrix(
    model: KineticAllosteryModel,
    epsilon: float = 1e-4,
) -> tuple[np.ndarray, float]:
    """Compute numerical Hessian (FIM) and condition number around current parameters.

    Returns
    -------
    fim : np.ndarray
        Hessian matrix (6x6).
    condition_number : float
        Ratio of largest to smallest eigenvalue.
    """
    keys = list(model.param_names)
    theta = np.array([getattr(model, k) for k in keys], dtype=float)
    n = len(theta)
    hessian = np.zeros((n, n), dtype=float)

    f0 = compute_kinetic_objective(theta, keys)

    # Finite difference second derivatives
    for i in range(n):
        for j in range(i, n):
            h_i = epsilon * max(abs(theta[i]), 1e-4)
            h_j = epsilon * max(abs(theta[j]), 1e-4)

            theta_pp = theta.copy()
            theta_pp[i] += h_i
            theta_pp[j] += h_j

            theta_pm = theta.copy()
            theta_pm[i] += h_i
            theta_pm[j] -= h_j

            theta_mp = theta.copy()
            theta_mp[i] -= h_i
            theta_mp[j] += h_j

            theta_mm = theta.copy()
            theta_mm[i] -= h_i
            theta_mm[j] -= h_j

            d2 = (
                compute_kinetic_objective(theta_pp, keys)
                - compute_kinetic_objective(theta_pm, keys)
                - compute_kinetic_objective(theta_mp, keys)
                + compute_kinetic_objective(theta_mm, keys)
            ) / (4.0 * h_i * h_j)

            hessian[i, j] = d2
            hessian[j, i] = d2

    eigvals = np.linalg.eigvalsh(hessian)
    pos_eigvals = eigvals[eigvals > 0]
    cond_num = float(np.max(pos_eigvals) / np.min(pos_eigvals)) if len(pos_eigvals) > 0 else float("inf")

    return hessian, cond_num
