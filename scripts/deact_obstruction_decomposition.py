#!/usr/bin/env python3
"""Deactivation obstruction decomposition across objectives and routes (Agent Brief 9).

Re-derives the amplitude obstruction by decomposing the published deactivation
objective into:
1. `all_six`: all 6 published quantities (tau1, tau2, tau3, area1, area2, area3)
2. `taus_only`: 3 time constants only (tau1, tau2, tau3)
3. `amplitudes_only`: 3 fractional amplitudes only (area1, area2, area3)

Each objective is evaluated twice:
- Route 1: `eigen`: analytic eigen-decomposition (`predict_components_eigen_all`)
  selecting the 3 largest |amplitude| modes.
- Route 2: `fit`: triexponential decay waveform fitting (`predict_components_fit`),
  the canonical experimental protocol route.

Where the fit route degenerates (is_well_conditioned is False, or repeated tau),
it is reported as DEGENERATE with its diagnostic reason, never backfilled from
the eigen route.

All literature target values and standard errors are read directly from
`circuitpharm.literature.DEACTIVATIONS[1]` without any numeric literals.

Usage:
    PYTHONPATH=src /Users/whinz/Biochemistry/.venv/bin/python scripts/deact_obstruction_decomposition.py
"""
from __future__ import annotations

import argparse
import pathlib
import sys
import time
from typing import Any, Sequence

import numpy as np
from scipy.optimize import minimize
from scipy.stats import qmc

# Ensure src is on sys.path
_ROOT = pathlib.Path(__file__).resolve().parent.parent
if str(_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(_ROOT / "src"))

from circuitpharm.fitting.comparison import MODEL_SPECS, ModelSpec
from circuitpharm.fitting.deact_components import (
    DeactFitResult,
    predict_components_eigen_all,
    predict_components_fit,
)
from circuitpharm.literature import DEACTIVATIONS


def load_literature_targets():
    """Read target deactivation components and published errors without literals."""
    record = DEACTIVATIONS[1]
    pub_taus = np.asarray(record.components[0], dtype=float)
    pub_areas = np.asarray(record.components[1], dtype=float)
    err_taus = np.asarray(record.component_errors[0], dtype=float)
    err_areas = np.asarray(record.component_errors[1], dtype=float)
    return pub_taus, pub_areas, err_taus, err_areas, record


def compute_chi2(
    taus_pred: np.ndarray,
    areas_pred: np.ndarray,
    pub_taus: np.ndarray,
    pub_areas: np.ndarray,
    err_taus: np.ndarray,
    err_areas: np.ndarray,
    objective: str,
) -> float:
    """Compute chi^2 for the specified objective."""
    z_t = (taus_pred - pub_taus) / err_taus
    z_a = (areas_pred - pub_areas) / err_areas

    if objective == "taus_only":
        return float(np.sum(z_t ** 2))
    elif objective == "amplitudes_only":
        return float(np.sum(z_a ** 2))
    elif objective == "all_six":
        return float(np.sum(z_t ** 2) + np.sum(z_a ** 2))
    else:
        raise ValueError(f"Unknown objective: {objective!r}")


def get_eigen_3_components(model: Any) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Extract top 3 modes by largest |amplitude|, returned sorted fast-to-slow."""
    taus_all, areas_all = predict_components_eigen_all(model)
    if len(taus_all) <= 3:
        order = np.argsort(taus_all)
        return taus_all[order], areas_all[order], taus_all, areas_all

    top3_idx = np.argsort(-np.abs(areas_all))[:3]
    sort_taus = np.argsort(taus_all[top3_idx])
    chosen = top3_idx[sort_taus]
    return taus_all[chosen], areas_all[chosen], taus_all, areas_all


def check_is_degenerate(fit_res: DeactFitResult) -> tuple[bool, str]:
    """Determine whether a triexponential fit is degenerate."""
    t1, t2, t3 = fit_res.taus
    collapsed_12 = (abs(t2 - t1) / max(t2, 1e-6)) < 0.05
    collapsed_23 = (abs(t3 - t2) / max(t3, 1e-6)) < 0.05

    reasons = []
    if not fit_res.is_well_conditioned:
        reasons.append(fit_res.degeneracy_reason or "ill-conditioned")
    if collapsed_12:
        reasons.append(f"tau1 and tau2 collapsed ({t1:.2f} ms vs {t2:.2f} ms)")
    if collapsed_23:
        reasons.append(f"tau2 and tau3 collapsed ({t2:.2f} ms vs {t3:.2f} ms)")

    is_deg = bool(reasons or not fit_res.is_well_conditioned or collapsed_12 or collapsed_23)
    return is_deg, "; ".join(dict.fromkeys(reasons)) if is_deg else ""


def run_eigen_decomposition(
    spec: ModelSpec,
    objective: str,
    pub_taus: np.ndarray,
    pub_areas: np.ndarray,
    err_taus: np.ndarray,
    err_areas: np.ndarray,
    n_samples: int = 1500,
    n_polish: int = 12,
    seed: int = 42,
    warm_starts: Sequence[np.ndarray] | None = None,
) -> dict[str, Any]:
    """Optimize the model under the eigen route for a specific objective."""
    lo = np.array([b[0] for b in spec.bounds], dtype=float)
    hi = np.array([b[1] for b in spec.bounds], dtype=float)
    d = len(spec.bounds)

    sampler = qmc.LatinHypercube(d=d, seed=seed)
    pts = lo + sampler.random(n_samples) * (hi - lo)

    all_starts = []
    if warm_starts:
        all_starts.extend(warm_starts)
    all_starts.extend(pts)

    evaluated = []
    for pt in all_starts:
        try:
            m = spec.factory(pt)
            t_pred, a_pred, _, _ = get_eigen_3_components(m)
            c2 = compute_chi2(t_pred, a_pred, pub_taus, pub_areas, err_taus, err_areas, objective)
            if np.isfinite(c2):
                evaluated.append((c2, pt))
        except Exception:
            continue

    if not evaluated:
        raise RuntimeError(f"No valid candidates evaluated for {spec.name} under {objective} (eigen)")

    evaluated.sort(key=lambda x: x[0])
    best_chi2, best_vec = evaluated[0]

    # Polish top candidates
    for _, pt in evaluated[:n_polish]:
        def obj_fn(v: np.ndarray) -> float:
            if np.any(v < lo) or np.any(v > hi):
                return 1e12
            try:
                m = spec.factory(v)
                t_pred, a_pred, _, _ = get_eigen_3_components(m)
                val = compute_chi2(t_pred, a_pred, pub_taus, pub_areas, err_taus, err_areas, objective)
                return val if np.isfinite(val) else 1e12
            except Exception:
                return 1e12

        res = minimize(obj_fn, pt, method="Nelder-Mead", options={"maxiter": 2000, "xatol": 1e-5, "fatol": 1e-5})
        if res.fun < best_chi2:
            best_chi2 = float(res.fun)
            best_vec = res.x

    opt_model = spec.factory(best_vec)
    t_pred, a_pred, all_t, all_a = get_eigen_3_components(opt_model)
    full_chi2 = compute_chi2(t_pred, a_pred, pub_taus, pub_areas, err_taus, err_areas, "all_six")
    tau_chi2 = compute_chi2(t_pred, a_pred, pub_taus, pub_areas, err_taus, err_areas, "taus_only")
    amp_chi2 = compute_chi2(t_pred, a_pred, pub_taus, pub_areas, err_taus, err_areas, "amplitudes_only")

    return {
        "spec_name": spec.name,
        "objective": objective,
        "route": "eigen",
        "best_chi2": float(best_chi2),
        "full_chi2": full_chi2,
        "tau_chi2": tau_chi2,
        "amp_chi2": amp_chi2,
        "best_vec": best_vec,
        "predicted_taus": t_pred,
        "predicted_areas": a_pred,
        "all_taus": all_t,
        "all_areas": all_a,
    }


def run_fit_decomposition(
    spec: ModelSpec,
    objective: str,
    pub_taus: np.ndarray,
    pub_areas: np.ndarray,
    err_taus: np.ndarray,
    err_areas: np.ndarray,
    n_samples: int = 150,
    n_polish: int = 6,
    seed: int = 42,
    warm_starts: Sequence[np.ndarray] | None = None,
) -> dict[str, Any]:
    """Optimize the model under the canonical fit route for a specific objective.

    Reports genuine well-conditioned fits or flags as DEGENERATE.
    """
    lo = np.array([b[0] for b in spec.bounds], dtype=float)
    hi = np.array([b[1] for b in spec.bounds], dtype=float)
    d = len(spec.bounds)

    sampler = qmc.LatinHypercube(d=d, seed=seed)
    pts = lo + sampler.random(n_samples) * (hi - lo)

    all_starts = []
    if warm_starts:
        all_starts.extend(warm_starts)
    all_starts.extend(pts)

    evaluated = []
    for pt in all_starts:
        try:
            m = spec.factory(pt)
            fit_res = predict_components_fit(m, n_starts=8)
            t_pred = np.asarray(fit_res.taus, dtype=float)
            a_pred = np.asarray(fit_res.areas, dtype=float)
            c2 = compute_chi2(t_pred, a_pred, pub_taus, pub_areas, err_taus, err_areas, objective)
            is_deg, deg_reason = check_is_degenerate(fit_res)
            if np.isfinite(c2):
                evaluated.append({
                    "chi2": c2,
                    "vec": pt,
                    "fit_res": fit_res,
                    "is_deg": is_deg,
                    "deg_reason": deg_reason,
                })
        except Exception:
            continue

    if not evaluated:
        return {
            "spec_name": spec.name,
            "objective": objective,
            "route": "fit",
            "is_degenerate": True,
            "degeneracy_reason": "No valid fit evaluations converged",
            "best_chi2": float("nan"),
            "display": "DEGENERATE (no converged fits)",
        }

    # Partition into well-conditioned and degenerate
    well_cond = [x for x in evaluated if not x["is_deg"]]
    evaluated.sort(key=lambda x: x["chi2"])

    # If any well-conditioned candidates exist, polish them with preference
    if well_cond:
        well_cond.sort(key=lambda x: x["chi2"])
        polish_pool = [x["vec"] for x in well_cond[:n_polish]]
    else:
        polish_pool = [x["vec"] for x in evaluated[:n_polish]]

    best_cand = well_cond[0] if well_cond else evaluated[0]
    best_chi2 = best_cand["chi2"]
    best_vec = best_cand["vec"]
    best_fit = best_cand["fit_res"]
    best_is_deg = best_cand["is_deg"]
    best_deg_reason = best_cand["deg_reason"]

    for pt in polish_pool:
        def obj_fn(v: np.ndarray) -> float:
            if np.any(v < lo) or np.any(v > hi):
                return 1e12
            try:
                m = spec.factory(v)
                fit = predict_components_fit(m, n_starts=6)
                t_p = np.asarray(fit.taus, dtype=float)
                a_p = np.asarray(fit.areas, dtype=float)
                val = compute_chi2(t_p, a_p, pub_taus, pub_areas, err_taus, err_areas, objective)
                is_d, _ = check_is_degenerate(fit)
                # Penalize degenerate fits heavily so optimizer seeks well-conditioned regions
                penalty = 1e6 if is_d else 0.0
                return val + penalty if np.isfinite(val) else 1e12
            except Exception:
                return 1e12

        res = minimize(obj_fn, pt, method="Nelder-Mead", options={"maxiter": 300, "xatol": 1e-4, "fatol": 1e-4})
        try:
            m_eval = spec.factory(res.x)
            fit_eval = predict_components_fit(m_eval, n_starts=12)
            c2_eval = compute_chi2(np.asarray(fit_eval.taus), np.asarray(fit_eval.areas),
                                   pub_taus, pub_areas, err_taus, err_areas, objective)
            is_d, d_reas = check_is_degenerate(fit_eval)
            if not is_d:
                if best_is_deg or c2_eval < best_chi2:
                    best_chi2 = c2_eval
                    best_vec = res.x
                    best_fit = fit_eval
                    best_is_deg = False
                    best_deg_reason = ""
            elif best_is_deg and c2_eval < best_chi2:
                best_chi2 = c2_eval
                best_vec = res.x
                best_fit = fit_eval
                best_is_deg = True
                best_deg_reason = d_reas
        except Exception:
            continue

    if best_is_deg:
        disp = f"DEGENERATE ({best_deg_reason or 'rank-deficient / repeated tau'})"
    else:
        disp = f"χ² = {best_chi2:.2f}"

    return {
        "spec_name": spec.name,
        "objective": objective,
        "route": "fit",
        "is_degenerate": best_is_deg,
        "degeneracy_reason": best_deg_reason,
        "best_chi2": float(best_chi2) if not best_is_deg else float("nan"),
        "raw_chi2": float(best_chi2),
        "best_vec": best_vec,
        "fit_res": best_fit,
        "display": disp,
    }


def format_decomposition_table(
    schemes: list[str],
    results_eigen: dict[tuple[str, str], dict[str, Any]],
    results_fit: dict[tuple[str, str], dict[str, Any]],
) -> str:
    """Format comparison table across schemes, objectives, and routes."""
    lines = []
    lines.append("| Scheme | Route | All Six ($\\chi^2$) | Time Constants Only ($\\chi^2$) | Amplitudes Only ($\\chi^2$) |")
    lines.append("| :--- | :---: | :---: | :---: | :---: |")

    for scheme in schemes:
        # Eigen row
        r_all_e = results_eigen[(scheme, "all_six")]
        r_tau_e = results_eigen[(scheme, "taus_only")]
        r_amp_e = results_eigen[(scheme, "amplitudes_only")]
        lines.append(
            f"| `{scheme}` | **Eigen** | **{r_all_e['best_chi2']:.2f}** | **{r_tau_e['best_chi2']:.2f}** | **{r_amp_e['best_chi2']:.2f}** |"
        )

        # Fit row
        r_all_f = results_fit[(scheme, "all_six")]
        r_tau_f = results_fit[(scheme, "taus_only")]
        r_amp_f = results_fit[(scheme, "amplitudes_only")]

        lines.append(
            f"| `{scheme}` | **Fit (canonical)** | {r_all_f['display']} | {r_tau_f['display']} | {r_amp_f['display']} |"
        )

    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(description="Deactivation obstruction decomposition")
    parser.add_argument("--schemes", nargs="+", default=[
        "kinetic_jw95_free_rates", "kinetic_flip", "extended_desens", "kinetic_open_desens",
    ], help="Schemes to evaluate (default: 4 required schemes)")
    parser.add_argument("--n-samples-eigen", type=int, default=1500, help="LHS samples for eigen route (default: 1500)")
    parser.add_argument("--n-polish-eigen", type=int, default=12, help="Polish count for eigen route (default: 12)")
    parser.add_argument("--n-samples-fit", type=int, default=100, help="LHS samples for fit route (default: 100)")
    parser.add_argument("--n-polish-fit", type=int, default=6, help="Polish count for fit route (default: 6)")
    parser.add_argument("--seed", type=int, default=42, help="Random seed (default: 42)")
    args = parser.parse_args()

    pub_taus, pub_areas, err_taus, err_areas, record = load_literature_targets()

    print("=" * 80)
    print("DEACTIVATION OBSTRUCTION DECOMPOSITION (Agent Brief 9 Task 2)")
    print("=" * 80)
    print(f"Literature Target: {record.label} (source: {record.source_key})")
    print(f"Published values: taus = {pub_taus} ms, areas = {pub_areas}")
    print(f"Published errors: err_taus = {err_taus} ms, err_areas = {err_areas}")
    print(f"Schemes: {args.schemes}")
    print(f"Eigen search budget: {args.n_samples_eigen} LHS, top {args.n_polish_eigen} polished")
    print(f"Fit search budget: {args.n_samples_fit} LHS, top {args.n_polish_fit} polished")
    print("-" * 80)

    t_start = time.perf_counter()

    results_eigen = {}
    results_fit = {}
    warm_starts_by_scheme = {}

    # 1. Run Eigen Route
    print("\n>>> PHASE 1: EIGEN ROUTE DECOMPOSITION")
    for scheme in args.schemes:
        spec = MODEL_SPECS[scheme]
        print(f"\n[Eigen] Evaluating `{scheme}` ({len(spec.param_names)} free params)...")
        warm = []
        for obj_name in ["all_six", "taus_only", "amplitudes_only"]:
            t0 = time.perf_counter()
            res_e = run_eigen_decomposition(
                spec, obj_name, pub_taus, pub_areas, err_taus, err_areas,
                n_samples=args.n_samples_eigen, n_polish=args.n_polish_eigen,
                seed=args.seed,
            )
            elapsed = time.perf_counter() - t0
            results_eigen[(scheme, obj_name)] = res_e
            warm.append(res_e["best_vec"])
            print(f"  * {obj_name:<16}: chi^2 = {res_e['best_chi2']:8.2f}  ({elapsed:.2f} s)")
        warm_starts_by_scheme[scheme] = warm

    # 2. Run Fit Route
    print("\n>>> PHASE 2: CANONICAL FIT ROUTE DECOMPOSITION")
    for scheme in args.schemes:
        spec = MODEL_SPECS[scheme]
        print(f"\n[Fit] Evaluating `{scheme}` ({len(spec.param_names)} free params)...")
        warm = warm_starts_by_scheme.get(scheme, [])
        for obj_name in ["all_six", "taus_only", "amplitudes_only"]:
            t0 = time.perf_counter()
            res_f = run_fit_decomposition(
                spec, obj_name, pub_taus, pub_areas, err_taus, err_areas,
                n_samples=args.n_samples_fit, n_polish=args.n_polish_fit,
                seed=args.seed, warm_starts=warm,
            )
            elapsed = time.perf_counter() - t0
            results_fit[(scheme, obj_name)] = res_f
            print(f"  * {obj_name:<16}: {res_f['display']}  ({elapsed:.2f} s)")

    t_total = time.perf_counter() - t_start
    print("\n" + "=" * 80)
    print(f"TOTAL EXECUTION TIME: {t_total:.2f} s")
    print("=" * 80 + "\n")

    md_table = format_decomposition_table(args.schemes, results_eigen, results_fit)
    print("## Obstruction Decomposition Summary Table\n")
    print(md_table)
    print("\n")


if __name__ == "__main__":
    main()
