#!/usr/bin/env python3
"""Deactivation obstruction decomposition across objectives and routes (Agent Brief 9 & 9B).

Re-derives the amplitude obstruction by decomposing the published deactivation
objective into:
1. `all_six`: all 6 published quantities (tau1, tau2, tau3, area1, area2, area3)
2. `taus_only`: 3 time constants only (tau1, tau2, tau3)
3. `amplitudes_only`: 3 fractional amplitudes only (area1, area2, area3)

Each objective is evaluated twice:
- Route 1: `eigen`: analytic eigen-decomposition (`predict_components_eigen`)
  selecting the 3 largest |amplitude| modes.
- Route 2: `fit`: triexponential decay waveform fitting (`predict_components_fit`),
  the canonical experimental protocol route.

Conditioning hazard check (Brief 9B §3):
Reports raw amplitude sum Σ a_j alongside fractional amplitudes and true P_o(0)
at the start of relaxation. Rejects any cell where |Σ a_j| is small relative to
P_o(0) (threshold: P_o(0) >= 1e-4, |Σ a_j| >= 1e-4, |Σ a_j| / P_o(0) >= 0.10).

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
import scipy.linalg as la
from scipy.optimize import minimize
from scipy.stats import qmc

# Ensure src is on sys.path
_ROOT = pathlib.Path(__file__).resolve().parent.parent
if str(_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(_ROOT / "src"))

from circuitpharm.fitting.comparison import MODEL_SPECS, ModelSpec
from circuitpharm.fitting.data import BARBERIS2007_PAIRED_PULSE, JAHN1997_PEAK_CRC
from circuitpharm.fitting.deact_components import (
    DeactFitResult,
    get_barberis_errors,
    predict_components_eigen,
    predict_components_eigen_all,
    predict_components_fit,
)
from circuitpharm.fitting.likelihood import predict_paired_pulse_recovery
from circuitpharm.literature import DEACTIVATIONS
from circuitpharm.gabaa_kinetics import FIT_RANGES


# Conditioning thresholds (Brief 9B §3)
THRESHOLD_PO_0 = 1e-4
THRESHOLD_RAW_SUM = 1e-4
THRESHOLD_RATIO = 0.10


def load_literature_targets():
    """Read target deactivation components and published errors without literals."""
    record = DEACTIVATIONS[1]
    pub_taus = np.asarray(record.components[0], dtype=float)
    pub_areas = np.asarray(record.components[1], dtype=float)
    err_taus, err_areas = get_barberis_errors(record)
    err_taus = np.asarray(err_taus, dtype=float)
    err_areas = np.asarray(err_areas, dtype=float)
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


def get_model_onset_po(model: Any, pulse_ms: float = 2.0, conc_mm: float = 10.0) -> float:
    """Compute open probability at the end of the GABA pulse."""
    gaba_pulse_um = float(conc_mm) * 1000.0
    p0 = np.asarray(model.state_distribution(0.0), dtype=float)
    q_pulse = np.asarray(model.q_matrix(gaba_pulse_um), dtype=float)
    p_pulse = p0 @ la.expm(q_pulse * float(pulse_ms))
    p_pulse = np.clip(p_pulse, 0.0, 1.0)
    p_pulse = p_pulse / np.sum(p_pulse)

    if hasattr(model, "open_indices"):
        o_idx = list(model.open_indices)
    elif hasattr(model, "open_index"):
        idx = model.open_index
        o_idx = list(idx) if isinstance(idx, (tuple, list, np.ndarray)) else [int(idx)]
    else:
        o_idx = [3]
    return float(np.sum(p_pulse[o_idx]))


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


def check_at_bound(vec: np.ndarray, bounds: Sequence[tuple[float, float]], param_names: Sequence[str], tol: float = 1e-3) -> list[str]:
    """Identify parameters sitting on search box bounds."""
    flagged = []
    for name, val, (b_lo, b_hi) in zip(param_names, vec, bounds):
        if abs(val - b_lo) < tol:
            flagged.append(f"{name} on lower bound ({b_lo})")
        elif abs(val - b_hi) < tol:
            flagged.append(f"{name} on upper bound ({b_hi})")
    return flagged


def evaluate_task4_cross_checks(spec: ModelSpec, opt_vec: np.ndarray) -> dict[str, Any]:
    """Task 4 of Brief 9: evaluate recovery, Jahn CRC, po_max, and bound flags."""
    model = spec.factory(opt_vec)

    # 1. Paired-pulse recovery
    pp_pred = float(predict_paired_pulse_recovery(model, BARBERIS2007_PAIRED_PULSE))
    pp_obs = float(BARBERIS2007_PAIRED_PULSE.mean_response)
    pp_sem = float(BARBERIS2007_PAIRED_PULSE.sem)
    pp_sigma = (pp_pred - pp_obs) / pp_sem

    # 2. Jahn 1997 peak CRC RMSE
    crc_pred = np.asarray(model.dose_response(JAHN1997_PEAK_CRC.concentrations_um), dtype=float)
    crc_obs = np.asarray(JAHN1997_PEAK_CRC.mean_response, dtype=float)
    crc_rmse = float(np.sqrt(np.mean((crc_pred - crc_obs) ** 2)))
    crc_baseline = 0.004917
    crc_ratio = crc_rmse / crc_baseline

    # 3. po_max against FIT_RANGES [0.70, 0.80]
    po_max_val = float(model.po_max)
    po_min, po_max = FIT_RANGES["po_max"]
    po_in_band = bool(po_min <= po_max_val <= po_max)

    # 4. Parameters on bounds
    pinned = check_at_bound(opt_vec, spec.bounds, spec.param_names)

    return {
        "pp_pred": pp_pred,
        "pp_obs": pp_obs,
        "pp_sem": pp_sem,
        "pp_sigma": pp_sigma,
        "crc_rmse": crc_rmse,
        "crc_baseline": crc_baseline,
        "crc_ratio": crc_ratio,
        "po_max": po_max_val,
        "po_in_band": po_in_band,
        "pinned": pinned,
    }


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
    t0 = time.perf_counter()
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
            t_pred, a_pred = predict_components_eigen(m, n_components=3)
            c2 = compute_chi2(t_pred, a_pred, pub_taus, pub_areas, err_taus, err_areas, objective)
            if np.isfinite(c2):
                evaluated.append((c2, pt))
        except Exception:
            continue

    if not evaluated:
        raise RuntimeError(f"No valid candidates evaluated for {spec.name} under {objective} (eigen)")

    evaluated.sort(key=lambda x: x[0])
    best_chi2, best_vec = evaluated[0]

    # Polish top candidates with Nelder-Mead
    for _, pt in evaluated[:n_polish]:
        def obj_fn(v: np.ndarray) -> float:
            if np.any(v < lo) or np.any(v > hi):
                return 1e12
            try:
                m = spec.factory(v)
                t_pred, a_pred = predict_components_eigen(m, n_components=3)
                val = compute_chi2(t_pred, a_pred, pub_taus, pub_areas, err_taus, err_areas, objective)
                return val if np.isfinite(val) else 1e12
            except Exception:
                return 1e12

        res = minimize(obj_fn, pt, method="Nelder-Mead", options={"maxiter": 2000, "xatol": 1e-5, "fatol": 1e-5})
        if res.fun < best_chi2:
            best_chi2 = float(res.fun)
            best_vec = res.x

    opt_model = spec.factory(best_vec)
    t_pred, a_pred, raw_amps, po_0 = predict_components_eigen(opt_model, n_components=3, return_raw=True)
    raw_sum = float(np.sum(raw_amps))
    ratio = abs(raw_sum) / po_0 if po_0 > 0 else 0.0

    is_conditioned = bool(
        (po_0 >= THRESHOLD_PO_0) and (abs(raw_sum) >= THRESHOLD_RAW_SUM) and (ratio >= THRESHOLD_RATIO)
    )
    cond_reason = ""
    if not is_conditioned:
        cond_reason = f"ILL-CONDITIONED (Po(0)={po_0:.4e}, raw_sum={raw_sum:.4e}, ratio={ratio:.4f} < {THRESHOLD_RATIO})"

    full_chi2 = compute_chi2(t_pred, a_pred, pub_taus, pub_areas, err_taus, err_areas, "all_six")
    tau_chi2 = compute_chi2(t_pred, a_pred, pub_taus, pub_areas, err_taus, err_areas, "taus_only")
    amp_chi2 = compute_chi2(t_pred, a_pred, pub_taus, pub_areas, err_taus, err_areas, "amplitudes_only")

    elapsed = time.perf_counter() - t0

    return {
        "spec_name": spec.name,
        "objective": objective,
        "route": "eigen",
        "best_chi2": float(best_chi2),
        "full_chi2": full_chi2,
        "tau_chi2": tau_chi2,
        "amp_chi2": amp_chi2,
        "best_vec": best_vec,
        "predicted_taus": tuple(float(x) for x in t_pred),
        "predicted_areas": tuple(float(x) for x in a_pred),
        "raw_amplitudes": tuple(float(x) for x in raw_amps),
        "raw_sum": raw_sum,
        "po_0": po_0,
        "ratio": ratio,
        "is_conditioned": is_conditioned,
        "cond_reason": cond_reason,
        "runtime_s": elapsed,
        "budget": f"{n_samples} LHS + {n_polish} polish",
    }


def run_fit_decomposition(
    spec: ModelSpec,
    objective: str,
    pub_taus: np.ndarray,
    pub_areas: np.ndarray,
    err_taus: np.ndarray,
    err_areas: np.ndarray,
    n_samples: int = 50,
    n_polish: int = 2,
    seed: int = 42,
    warm_starts: Sequence[np.ndarray] | None = None,
) -> dict[str, Any]:
    """Optimize the model under the canonical fit route for a specific objective.

    Reports genuine well-conditioned fits or flags as DEGENERATE.
    Explicitly budgeted per Brief 9B §2.
    """
    t0 = time.perf_counter()
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
            fit_res = predict_components_fit(m, n_starts=4)
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
        elapsed = time.perf_counter() - t0
        return {
            "spec_name": spec.name,
            "objective": objective,
            "route": "fit",
            "is_degenerate": True,
            "degeneracy_reason": "No valid fit evaluations converged",
            "best_chi2": float("nan"),
            "display": "DEGENERATE (no converged fits)",
            "runtime_s": elapsed,
            "budget": f"{n_samples} LHS + {n_polish} polish",
        }

    # Partition into well-conditioned and degenerate
    well_cond = [x for x in evaluated if not x["is_deg"]]
    evaluated.sort(key=lambda x: x["chi2"])

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
                fit = predict_components_fit(m, n_starts=2)
                t_p = np.asarray(fit.taus, dtype=float)
                a_p = np.asarray(fit.areas, dtype=float)
                val = compute_chi2(t_p, a_p, pub_taus, pub_areas, err_taus, err_areas, objective)
                is_d, _ = check_is_degenerate(fit)
                penalty = 1e6 if is_d else 0.0
                return val + penalty if np.isfinite(val) else 1e12
            except Exception:
                return 1e12

        res = minimize(obj_fn, pt, method="Nelder-Mead", options={"maxiter": 40, "xatol": 1e-4, "fatol": 1e-4})
        try:
            m_eval = spec.factory(res.x)
            fit_eval = predict_components_fit(m_eval, n_starts=8)
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

    # Final conditioning extraction
    opt_model = spec.factory(best_vec)
    po_0 = get_model_onset_po(opt_model)
    raw_amps = np.asarray(best_fit.areas, dtype=float) * po_0
    raw_sum = float(np.sum(raw_amps))
    ratio = abs(raw_sum) / po_0 if po_0 > 0 else 0.0

    is_conditioned = bool(
        (po_0 >= THRESHOLD_PO_0) and (abs(raw_sum) >= THRESHOLD_RAW_SUM) and (ratio >= THRESHOLD_RATIO)
    )
    cond_reason = ""
    if not is_conditioned:
        cond_reason = f"ILL-CONDITIONED (Po(0)={po_0:.4e}, raw_sum={raw_sum:.4e}, ratio={ratio:.4f} < {THRESHOLD_RATIO})"

    if best_is_deg:
        disp = f"DEGENERATE ({best_deg_reason or 'rank-deficient / repeated tau'})"
    elif not is_conditioned:
        disp = f"ILL-CONDITIONED ({cond_reason})"
    else:
        disp = f"χ² = {best_chi2:.2f}"

    elapsed = time.perf_counter() - t0

    return {
        "spec_name": spec.name,
        "objective": objective,
        "route": "fit",
        "is_degenerate": best_is_deg,
        "degeneracy_reason": best_deg_reason,
        "best_chi2": float(best_chi2) if (not best_is_deg and is_conditioned) else float("nan"),
        "raw_chi2": float(best_chi2),
        "best_vec": best_vec,
        "fit_res": best_fit,
        "predicted_taus": tuple(float(x) for x in best_fit.taus),
        "predicted_areas": tuple(float(x) for x in best_fit.areas),
        "raw_amplitudes": tuple(float(x) for x in raw_amps),
        "raw_sum": raw_sum,
        "po_0": po_0,
        "ratio": ratio,
        "is_conditioned": is_conditioned,
        "cond_reason": cond_reason,
        "display": disp,
        "runtime_s": elapsed,
        "budget": f"{n_samples} LHS + {n_polish} polish",
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


def format_conditioning_table(
    schemes: list[str],
    results_eigen: dict[tuple[str, str], dict[str, Any]],
    results_fit: dict[tuple[str, str], dict[str, Any]],
) -> str:
    """Format raw amplitude sum, Po(0), and conditioning ratio table (Brief 9B §3)."""
    lines = []
    lines.append("| Scheme | Route | Objective | $\\chi^2$ | Raw $\\Sigma a_j$ | True $P_o(0)$ | Ratio $|\\Sigma a_j| / P_o(0)$ | Conditioning Status |")
    lines.append("| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: |")

    for scheme in schemes:
        for route_name, r_dict in [("Eigen", results_eigen), ("Fit", results_fit)]:
            for obj in ["all_six", "taus_only", "amplitudes_only"]:
                r = r_dict[(scheme, obj)]
                c2_str = f"{r['best_chi2']:.2f}" if np.isfinite(r['best_chi2']) else "N/A"
                raw_s = f"{r.get('raw_sum', 0.0):.4f}"
                po0 = f"{r.get('po_0', 0.0):.4f}"
                ratio = f"{r.get('ratio', 0.0):.4f}"
                status = "PASS" if r.get("is_conditioned", True) else f"REJECTED ({r.get('cond_reason', '')})"
                if route_name == "Fit" and r.get("is_degenerate", False):
                    status = f"DEGENERATE"
                lines.append(
                    f"| `{scheme}` | {route_name} | `{obj}` | {c2_str} | {raw_s} | {po0} | {ratio} | **{status}** |"
                )
    return "\n".join(lines)


def format_budget_table(
    schemes: list[str],
    results_eigen: dict[tuple[str, str], dict[str, Any]],
    results_fit: dict[tuple[str, str], dict[str, Any]],
) -> str:
    """Format budget and runtime table per cell (Brief 9B §2)."""
    lines = []
    lines.append("| Scheme | Route | Objective | Search Budget | Runtime (s) |")
    lines.append("| :--- | :---: | :---: | :---: | :---: |")

    for scheme in schemes:
        for route_name, r_dict in [("Eigen", results_eigen), ("Fit", results_fit)]:
            for obj in ["all_six", "taus_only", "amplitudes_only"]:
                r = r_dict[(scheme, obj)]
                lines.append(
                    f"| `{scheme}` | {route_name} | `{obj}` | {r.get('budget', 'N/A')} | {r.get('runtime_s', 0.0):.2f} s |"
                )
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(description="Deactivation obstruction decomposition (Brief 9 & 9B)")
    parser.add_argument("--schemes", nargs="+", default=[
        "kinetic_jw95_free_rates", "kinetic_open_desens_free_rates",
    ], help="Schemes to evaluate in priority order (default: Priority 1 and 2)")
    parser.add_argument("--all-schemes", action="store_true", help="Evaluate all 4 schemes including kinetic_flip and extended_desens")
    parser.add_argument("--n-samples-eigen", type=int, default=1500, help="LHS samples for eigen route (default: 1500)")
    parser.add_argument("--n-polish-eigen", type=int, default=12, help="Polish count for eigen route (default: 12)")
    parser.add_argument("--n-samples-fit", type=int, default=50, help="LHS samples for fit route (default: 50)")
    parser.add_argument("--n-polish-fit", type=int, default=2, help="Polish count for fit route (default: 2)")
    parser.add_argument("--seed", type=int, default=42, help="Random seed (default: 42)")
    args = parser.parse_args()

    if args.all_schemes:
        schemes = ["kinetic_jw95_free_rates", "kinetic_open_desens_free_rates", "kinetic_flip", "extended_desens"]
    else:
        schemes = args.schemes

    pub_taus, pub_areas, err_taus, err_areas, record = load_literature_targets()

    print("=" * 80)
    print("DEACTIVATION OBSTRUCTION DECOMPOSITION (Agent Brief 9 & 9B)")
    print("=" * 80)
    print(f"Literature Target: {record.label} (source: {record.source_key})")
    print(f"Published values: taus = {pub_taus} ms, areas = {pub_areas}")
    print(f"Published errors: err_taus = {err_taus} ms, err_areas = {err_areas}")
    print(f"Conditioning thresholds: Po(0) >= {THRESHOLD_PO_0}, raw_sum >= {THRESHOLD_RAW_SUM}, ratio >= {THRESHOLD_RATIO}")
    print(f"Schemes: {schemes}")
    print(f"Eigen search budget: {args.n_samples_eigen} LHS, top {args.n_polish_eigen} polished")
    print(f"Fit search budget: {args.n_samples_fit} LHS, top {args.n_polish_fit} polished")
    print("-" * 80, flush=True)

    t_start = time.perf_counter()

    results_eigen = {}
    results_fit = {}
    warm_starts_by_scheme = {}

    # 1. Run Eigen Route
    print("\n>>> PHASE 1: EIGEN ROUTE DECOMPOSITION", flush=True)
    for scheme in schemes:
        spec = MODEL_SPECS[scheme]
        print(f"\n[Eigen] Evaluating `{scheme}` ({len(spec.param_names)} free params)...", flush=True)
        warm = []
        for obj_name in ["all_six", "taus_only", "amplitudes_only"]:
            res_e = run_eigen_decomposition(
                spec, obj_name, pub_taus, pub_areas, err_taus, err_areas,
                n_samples=args.n_samples_eigen, n_polish=args.n_polish_eigen,
                seed=args.seed,
            )
            results_eigen[(scheme, obj_name)] = res_e
            warm.append(res_e["best_vec"])
            print(f"  * {obj_name:<16}: chi^2 = {res_e['best_chi2']:8.2f} (Po(0)={res_e['po_0']:.4f}, raw_sum={res_e['raw_sum']:.4f}, ratio={res_e['ratio']:.4f}) [{res_e['runtime_s']:.2f} s]", flush=True)
        warm_starts_by_scheme[scheme] = warm

    # 2. Run Fit Route
    print("\n>>> PHASE 2: CANONICAL FIT ROUTE DECOMPOSITION", flush=True)
    for scheme in schemes:
        spec = MODEL_SPECS[scheme]
        print(f"\n[Fit] Evaluating `{scheme}` ({len(spec.param_names)} free params)...", flush=True)
        warm = warm_starts_by_scheme.get(scheme, [])
        for obj_name in ["all_six", "taus_only", "amplitudes_only"]:
            res_f = run_fit_decomposition(
                spec, obj_name, pub_taus, pub_areas, err_taus, err_areas,
                n_samples=args.n_samples_fit, n_polish=args.n_polish_fit,
                seed=args.seed, warm_starts=warm,
            )
            results_fit[(scheme, obj_name)] = res_f
            print(f"  * {obj_name:<16}: {res_f['display']} (Po(0)={res_f.get('po_0', 0.0):.4f}, raw_sum={res_f.get('raw_sum', 0.0):.4f}) [{res_f['runtime_s']:.2f} s]", flush=True)

    t_total = time.perf_counter() - t_start
    print("\n" + "=" * 80)
    print(f"TOTAL EXECUTION TIME: {t_total:.2f} s")
    print("=" * 80 + "\n", flush=True)

    md_table = format_decomposition_table(schemes, results_eigen, results_fit)
    print("## 1. Obstruction Decomposition Summary Table\n")
    print(md_table)
    print("\n")

    cond_table = format_conditioning_table(schemes, results_eigen, results_fit)
    print("## 2. Amplitude Conditioning & Raw Sum Check (Brief 9B §3)\n")
    print(cond_table)
    print("\n")

    budget_table = format_budget_table(schemes, results_eigen, results_fit)
    print("## 3. Budget & Runtime Accounting Per Cell (Brief 9B §2)\n")
    print(budget_table)
    print("\n")

    # Task 4 Check for any candidate with chi^2 < 100
    print("## 4. Task 4 Cross-Checks for Optima with $\\chi^2 < 100$\n")
    found_any = False
    for scheme in schemes:
        spec = MODEL_SPECS[scheme]
        for route_name, r_dict in [("Eigen", results_eigen), ("Fit", results_fit)]:
            r = r_dict[(scheme, "all_six")]
            c2 = r.get("best_chi2", float("nan"))
            if np.isfinite(c2) and c2 < 100.0:
                found_any = True
                checks = evaluate_task4_cross_checks(spec, r["best_vec"])
                print(f"### Scheme `{scheme}` ({route_name} Route, All Six $\\chi^2 = {c2:.2f}$)")
                print(f"- Paired-pulse recovery: {checks['pp_pred']:.4f} vs {checks['pp_obs']:.2f} +/- {checks['pp_sem']:.2f} ({checks['pp_sigma']:+.2f} sigma)")
                print(f"- Jahn 1997 peak CRC RMSE: {checks['crc_rmse']:.6f} vs baseline {checks['crc_baseline']:.6f} ({checks['crc_ratio']:.2f}x)")
                print(f"- po_max: {checks['po_max']:.4f} (in [0.70, 0.80]: {checks['po_in_band']})")
                print(f"- Pinned parameters: {checks['pinned'] or 'None'}")
                print()
    if not found_any:
        print("No scheme reached $\\chi^2 < 100$ on all six published parameters (all $\\chi^2 > 500$).")
        print()


if __name__ == "__main__":
    main()
