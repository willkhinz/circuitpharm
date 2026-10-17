#!/usr/bin/env python3
"""Symmetric Joint Four-Way Fit Evaluation (Agent Brief 14 Task 1).

Evaluates whether any parameter vector of kinetic_jw95_free_rates, kinetic_two_open,
or kinetic_coop can simultaneously satisfy all four observables within 3 sigma,
with EVERY parametric dataset scored against its published parameters (Weighting iv extended):

1. Deactivation: 6 published components at published errors (tau: +/-0.3, +/-4.6, +/-14.9 ms;
   areas: rank-2 covariance from +/-0.04, +/-0.02, +/-0.03) under L1_MAX <= 1.5.
2. Peak CRC: EC50 = 11.6 +/- 0.9 uM, nH = 2.2 +/- 0.4 over 1-10 uM (Jahn et al. 1997).
   Retains 9-point reconstruction RMSE and ratio to baseline (0.004917) as a reported view.
3. Paired-pulse recovery: 0.33 +/- 0.03 (Barberis et al. 2007).
4. P_o,max gating anchor: band [0.70, 0.80] (target 0.75 +/- 0.05).

NO INVENTED WEIGHTING:
Every quantity has a published error. Residuals are evaluated and reported in sigma.

Objectives:
- obj_natural_joint: sum of squared normalized residuals across all four observables.
- obj_fourway_minimax: maximum normalized residual in sigma across all four observable families.
- obj_two_observable_pareto: deactivation vs CRC minimax trade-off frontier.

Nesting property:
kinetic_coop strictly nests kinetic_two_open at c = 1.0. Warm-start is enforced from
kinetic_two_open's optimum via nesting_embedding, and a strict assertion prevents
reporting any nesting violation.

Usage:
    PYTHONPATH=src /Users/whinz/Biochemistry/.venv/bin/python scripts/joint_fourway_fit.py
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

from circuitpharm.fitting.comparison import MODEL_SPECS, ModelSpec, nesting_embedding
from circuitpharm.fitting.crc_components import (
    CRC_BASELINE_RMSE,
    crc_components_log_likelihood,
    extract_crc_parameters,
    get_jahn_published_records,
    surface_ec50_discrepancy,
)
from circuitpharm.fitting.data import (
    BARBERIS2007_PAIRED_PULSE,
    JAHN1997_PEAK_CRC,
    JAHN_EC50_SEM,
    JAHN_EC50_UM,
    JAHN_HILL,
    JAHN_HILL_SEM,
    JAHN_SLOPE_WINDOW_UM,
)
from circuitpharm.fitting.deact_components import (
    DEFAULT_L1_MAX,
    check_amplitude_admissibility,
    deact_components_log_likelihood,
    predict_components_eigen,
)
from circuitpharm.fitting.likelihood import predict_paired_pulse_recovery
from circuitpharm.gabaa_kinetics import FIT_RANGES, FIT_TARGETS

TARGET_PP_RECOVERY = float(BARBERIS2007_PAIRED_PULSE.mean_response)
SEM_PP_RECOVERY = float(BARBERIS2007_PAIRED_PULSE.sem)
TARGET_PO_MAX = float(FIT_TARGETS["po_max"])
SIGMA_PO_MAX = 0.05  # conventional calibration SEM (half-width of [0.70, 0.80])


def check_at_bound(
    vec: np.ndarray,
    bounds: Sequence[tuple[float, float]],
    param_names: Sequence[str],
    tol: float = 1e-3,
) -> list[str]:
    """Identify parameters sitting on search box bounds."""
    flagged = []
    for name, val, (b_lo, b_hi) in zip(param_names, vec, bounds):
        if abs(val - b_lo) < tol:
            flagged.append(f"{name} on lower bound ({b_lo})")
        elif abs(val - b_hi) < tol:
            flagged.append(f"{name} on upper bound ({b_hi})")
    return flagged


def evaluate_candidate(
    spec: ModelSpec,
    vec: np.ndarray,
    l1_max: float = DEFAULT_L1_MAX,
    compute_recovery: bool = True,
) -> dict[str, Any]:
    """Score one candidate vector across all four observables symmetrically without invented weighting."""
    v = np.asarray(vec, dtype=float)
    m = spec.factory(v)

    # 1. Deactivation components (Weighting iv: 6 published parameters at published errors)
    try:
        taus, areas = predict_components_eigen(m, n_components=3, l1_max=None)
        is_adm, c_ratio = check_amplitude_admissibility(areas, l1_max=l1_max)
        ll_deact = deact_components_log_likelihood(taus, areas, l1_max=None)
        chi2_taus = float(ll_deact.chi2_taus)
        chi2_areas = float(ll_deact.chi2_areas)
        z_taus = tuple(float(x) for x in ll_deact.z_taus)
        z_areas = tuple(float(x) for x in ll_deact.z_areas)
        deact_ok = True
    except Exception:
        deact_ok = False
        taus, areas = (0.0, 0.0, 0.0), (0.0, 0.0, 0.0)
        is_adm, c_ratio = False, 999.0
        chi2_taus, chi2_areas = 1e8, 1e8
        z_taus = (999.0, 999.0, 999.0)
        z_areas = (999.0, 999.0, 999.0)

    # 2. Peak CRC (Weighting iv: EC50 = 11.6 +/- 0.9 uM, nH = 2.2 +/- 0.4 over 1-10 uM)
    try:
        crc_ll = crc_components_log_likelihood(m, compute_9point_view=True)
        ec50_um = crc_ll.ec50_um
        hill_meas = crc_ll.hill_slope
        hill_full = crc_ll.hill_slope_full
        hill_rel = crc_ll.hill_slope_rel
        z_ec50 = crc_ll.z_ec50
        z_hill = crc_ll.z_hill
        z_hill_full = crc_ll.z_hill_full
        z_hill_rel = crc_ll.z_hill_rel
        chi2_ec50 = crc_ll.chi2_ec50
        chi2_hill = crc_ll.chi2_hill
        chi2_crc = crc_ll.chi2_total
        crc_rmse = crc_ll.crc_rmse
        crc_ratio = crc_ll.crc_ratio
    except Exception:
        ec50_um = float("nan")
        hill_meas = float("nan")
        hill_full = float("nan")
        hill_rel = float("nan")
        z_ec50 = 999.0
        z_hill = 999.0
        z_hill_full = 999.0
        z_hill_rel = 999.0
        chi2_ec50 = 1e8
        chi2_hill = 1e8
        chi2_crc = 1e8
        crc_rmse = 999.0
        crc_ratio = 999.0

    # 3. Paired-pulse recovery (0.33 +/- 0.03 at 100 ms gap)
    if compute_recovery:
        try:
            pp_pred = float(predict_paired_pulse_recovery(m, BARBERIS2007_PAIRED_PULSE))
            z_pp = float((pp_pred - TARGET_PP_RECOVERY) / SEM_PP_RECOVERY)
        except Exception:
            pp_pred = -1.0
            z_pp = 999.0
    else:
        pp_pred = float("nan")
        z_pp = float("nan")

    # 4. P_o,max gating anchor (band [0.70, 0.80], target 0.75 +/- 0.05)
    try:
        po_max = float(m.po_max)
        z_po = float((po_max - TARGET_PO_MAX) / SIGMA_PO_MAX)
        po_in_band = bool(FIT_RANGES["po_max"][0] <= po_max <= FIT_RANGES["po_max"][1])
    except Exception:
        po_max = -1.0
        z_po = 999.0
        po_in_band = False

    pinned = check_at_bound(v, spec.bounds, spec.param_names)
    k_eff = len(spec.param_names) - len(pinned)

    # Normalized 3-sigma violation across all four observables
    # Satisfying all 4 within 3 sigma requires max_violation <= 1.000
    terms = [
        abs(z_taus[0]) / 3.0,
        abs(z_taus[1]) / 3.0,
        abs(z_taus[2]) / 3.0,
        np.sqrt(max(chi2_areas, 0.0)) / 3.0,
        abs(z_ec50) / 3.0,
        abs(z_hill) / 3.0,
        (abs(z_pp) / 3.0) if np.isfinite(z_pp) else 0.0,
        abs(z_po) / 3.0,
    ]
    max_violation = float(max(terms)) if (is_adm and deact_ok) else 999.0

    return {
        "vector": v,
        "deact_ok": deact_ok,
        "is_adm": is_adm,
        "c_ratio": c_ratio,
        "taus": taus,
        "areas": areas,
        "z_taus": z_taus,
        "z_areas": z_areas,
        "chi2_taus": chi2_taus,
        "chi2_areas": chi2_areas,
        "chi2_deact": chi2_taus + chi2_areas,
        "ec50_um": ec50_um,
        "hill_meas": hill_meas,
        "hill_full": hill_full,
        "hill_rel": hill_rel,
        "z_ec50": z_ec50,
        "z_hill": z_hill,
        "z_hill_full": z_hill_full,
        "z_hill_rel": z_hill_rel,
        "chi2_ec50": chi2_ec50,
        "chi2_hill": chi2_hill,
        "chi2_crc": chi2_crc,
        "crc_rmse": crc_rmse,
        "crc_ratio": crc_ratio,
        "pp_pred": pp_pred,
        "z_pp": z_pp,
        "po_max": po_max,
        "z_po": z_po,
        "po_in_band": po_in_band,
        "pinned": pinned,
        "k_eff": k_eff,
        "max_violation": max_violation,
    }


def fit_joint_fourway(
    spec: ModelSpec,
    n_samples: int = 2000,
    n_polish: int = 12,
    seed: int = 42,
    l1_max: float = DEFAULT_L1_MAX,
    warm_starts: dict[str, np.ndarray] | None = None,
) -> dict[str, Any]:
    """Execute multi-objective and feasibility evaluation for a model scheme across all four observables."""
    t0 = time.perf_counter()
    lo = np.array([b[0] for b in spec.bounds], dtype=float)
    hi = np.array([b[1] for b in spec.bounds], dtype=float)
    d = len(spec.bounds)

    sampler = qmc.LatinHypercube(d=d, seed=seed)
    pts = lo + sampler.random(n_samples) * (hi - lo)

    # 1. Single-objective: Deactivation alone
    def obj_deact(v):
        try:
            m = spec.factory(v)
            taus, areas = predict_components_eigen(m, n_components=3, l1_max=None)
            is_adm, c = check_amplitude_admissibility(areas, l1_max=l1_max)
            if not is_adm:
                return 1e6 + c * 100
            ll = deact_components_log_likelihood(taus, areas, l1_max=None)
            return float(ll.chi2_total)
        except Exception:
            return 1e8

    # 2. Single-objective: Sourced Peak CRC alone (EC50 = 11.6 +/- 0.9 uM, nH = 2.2 +/- 0.4)
    def obj_crc(v):
        try:
            m = spec.factory(v)
            p = extract_crc_parameters(m)
            z_e = (p.ec50_um - JAHN_EC50_UM) / JAHN_EC50_SEM
            z_h = (p.hill_measured - JAHN_HILL) / JAHN_HILL_SEM
            return float(z_e ** 2 + z_h ** 2)
        except Exception:
            return 1e6

    # 3. Two-observable Pareto trade-off: Deact vs CRC minimax trade-off (relabelled from Brief 12/13)
    def obj_two_observable_pareto(v):
        """Two-observable Pareto trade-off between Deactivation and CRC (Brief 12/13)."""
        try:
            m = spec.factory(v)
            taus, areas = predict_components_eigen(m, n_components=3, l1_max=None)
            is_adm, c = check_amplitude_admissibility(areas, l1_max=l1_max)
            if not is_adm:
                return 1e4 + c * 10
            ll = deact_components_log_likelihood(taus, areas, l1_max=None)
            z_t = np.asarray(ll.z_taus)
            c2_a = float(ll.chi2_areas)
            max_z_deact = max(float(np.max(np.abs(z_t))), float(np.sqrt(max(c2_a, 0.0))))

            p = extract_crc_parameters(m)
            z_e = abs(p.ec50_um - JAHN_EC50_UM) / JAHN_EC50_SEM
            z_h = abs(p.hill_measured - JAHN_HILL) / JAHN_HILL_SEM
            max_z_crc = max(float(z_e), float(z_h))

            return float(max(max_z_deact / 3.0, max_z_crc / 3.0))
        except Exception:
            return 1e6

    # 4. Natural published-error sum of squares (including recovery)
    # Deactivation: chi2_taus + chi2_areas (published errors, 5 terms)
    # CRC: z_ec50^2 + z_nh^2 (published errors, 2 terms)
    # Recovery: z_pp^2 (Barberis 2007 paired-pulse recovery at 100 ms gap, 1 term)
    # Po_max: z_po^2 (gating anchor 0.75 +/- 0.05, 1 term)
    # Total n_data = 9
    def obj_natural_joint(v):
        try:
            m = spec.factory(v)
            taus, areas = predict_components_eigen(m, n_components=3, l1_max=None)
            is_adm, c = check_amplitude_admissibility(areas, l1_max=l1_max)
            if not is_adm:
                return 1e8 + c * 1e4
            ll = deact_components_log_likelihood(taus, areas, l1_max=None)
            c2_deact = float(ll.chi2_total)

            p = extract_crc_parameters(m)
            z_e = (p.ec50_um - JAHN_EC50_UM) / JAHN_EC50_SEM
            z_h = (p.hill_measured - JAHN_HILL) / JAHN_HILL_SEM
            c2_crc = float(z_e ** 2 + z_h ** 2)

            pp_pred = float(predict_paired_pulse_recovery(m, BARBERIS2007_PAIRED_PULSE))
            z_pp = (pp_pred - TARGET_PP_RECOVERY) / SEM_PP_RECOVERY
            c2_pp = float(z_pp ** 2)

            po_max = float(m.po_max)
            z_po = (po_max - TARGET_PO_MAX) / SIGMA_PO_MAX
            c2_po = float(z_po ** 2)

            return float(c2_deact + c2_crc + c2_pp + c2_po)
        except Exception:
            return 1e8

    # 5. True four-way minimax objective (Brief 14 Task 1)
    # Computes maximum normalized residual in sigma across all four observable families:
    # 1. Deactivation: max(|z_tau1|, |z_tau2|, |z_tau3|, sqrt(chi2_areas))
    # 2. Jahn CRC: max(|z_EC50|, |z_nH|)
    # 3. Paired-pulse recovery: |z_pp|
    # 4. Peak open probability: |z_po|
    def obj_fourway_minimax(v):
        try:
            m = spec.factory(v)
            taus, areas = predict_components_eigen(m, n_components=3, l1_max=None)
            is_adm, c = check_amplitude_admissibility(areas, l1_max=l1_max)
            if not is_adm:
                return 1e4 + c * 10
            ll = deact_components_log_likelihood(taus, areas, l1_max=None)
            z_t = np.asarray(ll.z_taus)
            c2_a = float(ll.chi2_areas)
            max_z_deact = max(float(np.max(np.abs(z_t))), float(np.sqrt(max(c2_a, 0.0))))

            p = extract_crc_parameters(m)
            z_e = abs(p.ec50_um - JAHN_EC50_UM) / JAHN_EC50_SEM
            z_h = abs(p.hill_measured - JAHN_HILL) / JAHN_HILL_SEM
            max_z_crc = max(float(z_e), float(z_h))

            pp_pred = float(predict_paired_pulse_recovery(m, BARBERIS2007_PAIRED_PULSE))
            z_pp = abs(pp_pred - TARGET_PP_RECOVERY) / SEM_PP_RECOVERY

            po_max = float(m.po_max)
            z_po = abs(po_max - TARGET_PO_MAX) / SIGMA_PO_MAX

            # Minimax residual across all four observable families (in units of sigma)
            return float(max(max_z_deact, max_z_crc, z_pp, z_po))
        except Exception:
            return 1e6

    print(f"[{spec.name}] Optimizing Deactivation alone...", flush=True)
    best_d_start = min(pts, key=obj_deact)
    if warm_starts and "deact" in warm_starts:
        if obj_deact(warm_starts["deact"]) < obj_deact(best_d_start):
            best_d_start = warm_starts["deact"]
    opt_deact = minimize(obj_deact, best_d_start, bounds=spec.bounds, method="L-BFGS-B")

    print(f"[{spec.name}] Optimizing Sourced Peak CRC alone...", flush=True)
    best_c_start = min(pts, key=obj_crc)
    if warm_starts and "crc" in warm_starts:
        if obj_crc(warm_starts["crc"]) < obj_crc(best_c_start):
            best_c_start = warm_starts["crc"]
    opt_crc = minimize(obj_crc, best_c_start, bounds=spec.bounds, method="L-BFGS-B")

    print(f"[{spec.name}] Optimizing Two-Observable Pareto Trade-off...", flush=True)
    pair_starts = [opt_deact.x, opt_crc.x] + sorted(pts, key=obj_two_observable_pareto)[:n_polish]
    if warm_starts and "pair" in warm_starts:
        pair_starts = [warm_starts["pair"]] + pair_starts
    pair_optima = []
    for s in pair_starts:
        res_p = minimize(obj_two_observable_pareto, s, bounds=spec.bounds, method="Nelder-Mead", options={"maxiter": 1200})
        pair_optima.append((obj_two_observable_pareto(res_p.x), res_p.x))
    pair_optima.sort(key=lambda x: x[0])
    best_pair_vec = pair_optima[0][1]
    best_pair_val = pair_optima[0][0]

    print(f"[{spec.name}] Optimizing Natural Published-Error Joint Objective (with recovery)...", flush=True)
    joint_starts = [best_pair_vec, opt_deact.x, opt_crc.x] + sorted(pts, key=obj_natural_joint)[:n_polish]
    if warm_starts and "joint" in warm_starts:
        joint_starts = [warm_starts["joint"]] + joint_starts
    joint_optima = []
    for s in joint_starts:
        res_j = minimize(obj_natural_joint, s, bounds=spec.bounds, method="Nelder-Mead", options={"maxiter": 1200})
        joint_optima.append((obj_natural_joint(res_j.x), res_j.x))
    joint_optima.sort(key=lambda x: x[0])
    best_joint_vec = joint_optima[0][1]
    best_joint_val = joint_optima[0][0]

    print(f"[{spec.name}] Optimizing True Four-Way Minimax Objective...", flush=True)
    minimax_starts = [best_pair_vec, best_joint_vec, opt_deact.x, opt_crc.x] + sorted(pts, key=obj_fourway_minimax)[:n_polish]
    if warm_starts and "minimax" in warm_starts:
        minimax_starts = [warm_starts["minimax"]] + minimax_starts
    minimax_optima = []
    for s in minimax_starts:
        res_m = minimize(obj_fourway_minimax, s, bounds=spec.bounds, method="Nelder-Mead", options={"maxiter": 1200})
        minimax_optima.append((obj_fourway_minimax(res_m.x), res_m.x))
    minimax_optima.sort(key=lambda x: x[0])
    best_minimax_vec = minimax_optima[0][1]
    best_minimax_val = minimax_optima[0][0]

    # Evaluate complete cross-checks for all profiles
    eval_deact_opt = evaluate_candidate(spec, opt_deact.x, l1_max=l1_max, compute_recovery=True)
    eval_crc_opt = evaluate_candidate(spec, opt_crc.x, l1_max=l1_max, compute_recovery=True)
    eval_pair_tradeoff = evaluate_candidate(spec, best_pair_vec, l1_max=l1_max, compute_recovery=True)
    eval_joint_opt = evaluate_candidate(spec, best_joint_vec, l1_max=l1_max, compute_recovery=True)
    eval_minimax_opt = evaluate_candidate(spec, best_minimax_vec, l1_max=l1_max, compute_recovery=True)

    runtime_s = time.perf_counter() - t0

    return {
        "scheme": spec.name,
        "runtime_s": runtime_s,
        "deact_opt": eval_deact_opt,
        "crc_opt": eval_crc_opt,
        "pair_tradeoff": eval_pair_tradeoff,
        "joint_opt": eval_joint_opt,
        "minimax_opt": eval_minimax_opt,
        "pair_opt_val": best_pair_val,
        "joint_opt_val": best_joint_val,
        "minimax_opt_val": best_minimax_val,
    }


def format_joint_results_table(results: dict[str, dict[str, Any]]) -> str:
    """Format markdown summary table of per-dataset residuals in sigma."""
    lines = [
        "| Scheme | Evaluation Profile | $z_{\\tau_1}$ | $z_{\\tau_2}$ | $z_{\\tau_3}$ | $z_{a_1}$ | $z_{a_2}$ | $z_{a_3}$ | $\\sqrt{\\chi^2_a}$ | $z_{\\text{EC50}}$ | $z_{n_H}$ | $z_{pp}$ (Recovery) | $P_{o,max}$ ($z_{po}$) | CRC RMSE Ratio (Rep.) | Max Viol / $3\\sigma$ | Satisfies All 4 $\\le 3\\sigma$? |",
        "| :--- | :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |",
    ]

    for scheme, r in results.items():
        profiles = [
            ("Deact-Alone Optimum", r["deact_opt"]),
            ("CRC-Alone Optimum", r["crc_opt"]),
            ("Two-Observable Pareto Trade-off", r["pair_tradeoff"]),
            ("Natural Joint Optimum", r["joint_opt"]),
            ("Four-Way Minimax Optimum", r["minimax_opt"]),
        ]
        for prof_name, ev in profiles:
            zt = ev["z_taus"]
            za = ev["z_areas"]
            sq_c2a = np.sqrt(max(ev["chi2_areas"], 0.0))
            zec50 = ev["z_ec50"]
            znh = ev["z_hill"]
            zpp = ev["z_pp"]
            zpp_str = f"{zpp:+.2f}$\\sigma$" if np.isfinite(zpp) else "N/A"
            crc_str = f"{ev['crc_ratio']:.2f}$\\times$" if np.isfinite(ev["crc_ratio"]) else "N/A"
            po_str = f"{ev['po_max']:.3f} ({ev['z_po']:+.2f}$\\sigma$)"
            viol = ev["max_violation"]
            sat_str = "**YES**" if viol <= 1.0 else "**NO**"

            lines.append(
                f"| `{scheme}` | {prof_name} | {zt[0]:+.2f}$\\sigma$ | {zt[1]:+.2f}$\\sigma$ | {zt[2]:+.2f}$\\sigma$ | "
                f"{za[0]:+.2f}$\\sigma$ | {za[1]:+.2f}$\\sigma$ | {za[2]:+.2f}$\\sigma$ | {sq_c2a:.2f}$\\sigma$ | "
                f"{zec50:+.2f}$\\sigma$ | {znh:+.2f}$\\sigma$ | "
                f"{zpp_str} | {po_str} | {crc_str} | **{viol:.2f}** | {sat_str} |"
            )

    return "\n".join(lines)


def format_parameters_table(results: dict[str, dict[str, Any]]) -> str:
    """Format parameter vectors, bound flags, and k_effective for joint optima."""
    lines = [
        "| Scheme | Profile | Free Parameters Vector ($\\log_{10}$) | Bounds Status / Pinned Parameters | $k_{fitted}$ | $k_{eff}$ |",
        "| :--- | :--- | :--- | :--- | :---: | :---: |",
    ]
    for scheme, r in results.items():
        spec = MODEL_SPECS[scheme]
        for prof_name, key in [
            ("Two-Observable Pareto Trade-off", "pair_tradeoff"),
            ("Natural Joint Optimum", "joint_opt"),
            ("Four-Way Minimax Optimum", "minimax_opt"),
        ]:
            ev = r[key]
            v_str = ", ".join(f"{name}={val:.3f}" for name, val in zip(spec.param_names, ev["vector"]))
            pinned_str = "; ".join(ev["pinned"]) if ev["pinned"] else "None (all interior)"
            lines.append(
                f"| `{scheme}` | {prof_name} | `{v_str}` | {pinned_str} | {len(spec.param_names)} | {ev['k_eff']} |"
            )
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(description="Symmetric Joint four-way fit evaluation (Brief 14 Task 1).")
    parser.add_argument("--samples", type=int, default=2000, help="Number of LHS samples (default: 2000).")
    parser.add_argument("--polish", type=int, default=12, help="Number of polished starts (default: 12).")
    parser.add_argument("--seed", type=int, default=42, help="Random seed (default: 42).")
    parser.add_argument(
        "--schemes",
        nargs="+",
        default=["kinetic_jw95_free_rates", "kinetic_two_open", "kinetic_coop"],
        help="Schemes to evaluate (default: kinetic_jw95_free_rates kinetic_two_open kinetic_coop).",
    )
    args = parser.parse_args()

    print("=" * 80)
    print("SYMMETRIC JOINT FOUR-WAY FIT EVALUATION (Agent Brief 14)")
    print("=" * 80)
    print(f"Schemes: {args.schemes}")
    print(f"Budget per scheme: {args.samples} LHS samples, top {args.polish} polished starts")
    print(f"Target 1: Deactivation components (6 published parameters at published errors, L1 <= {DEFAULT_L1_MAX})")
    print(f"Target 2: Sourced Peak CRC (EC50 = {JAHN_EC50_UM} +/- {JAHN_EC50_SEM} uM, nH = {JAHN_HILL} +/- {JAHN_HILL_SEM} over 1-10 uM)")
    print(f"Target 3: Paired-pulse recovery ({TARGET_PP_RECOVERY} +/- {SEM_PP_RECOVERY} at 100 ms gap)")
    print(f"Target 4: P_o,max gating anchor (band {FIT_RANGES['po_max']}, target {TARGET_PO_MAX} +/- {SIGMA_PO_MAX})")
    print("NO INVENTED WEIGHTING: reporting per-dataset residuals in sigma.")
    print("-" * 80, flush=True)

    # Discrepancy report
    disc = surface_ec50_discrepancy()
    print("EC50 DISCREPANCY REGISTER:")
    print(f"  * FIT_TARGETS['ec50_um']: {disc['fit_targets_ec50_um']:.1f} uM ({disc['fit_targets_source']})")
    print(f"  * Jahn 1997 sourced:      {disc['jahn_sourced_ec50_um']:.1f} +/- {disc['jahn_sourced_ec50_sem']:.1f} uM ({disc['jahn_source']})")
    print(f"  * Discrepancy:            {disc['discrepancy_um']:.1f} uM ({disc['discrepancy_sigma']:+.2f} sigma)")
    print(f"  * Resolution:             {disc['resolution']}")
    print("-" * 80, flush=True)

    results = {}
    for scheme in args.schemes:
        spec = MODEL_SPECS[scheme]
        warm_starts = None
        if scheme == "kinetic_coop" and "kinetic_two_open" in results:
            spec_sub = MODEL_SPECS["kinetic_two_open"]
            res_sub = results["kinetic_two_open"]
            warm_starts = {
                "deact": nesting_embedding(spec_sub, spec, res_sub["deact_opt"]["vector"]),
                "crc": nesting_embedding(spec_sub, spec, res_sub["crc_opt"]["vector"]),
                "pair": nesting_embedding(spec_sub, spec, res_sub["pair_tradeoff"]["vector"]),
                "joint": nesting_embedding(spec_sub, spec, res_sub["joint_opt"]["vector"]),
                "minimax": nesting_embedding(spec_sub, spec, res_sub["minimax_opt"]["vector"]),
            }
            print(f"[{scheme}] Warm-starting from kinetic_two_open optima via nesting_embedding (c = 1.0)...", flush=True)

        res = fit_joint_fourway(
            spec,
            n_samples=args.samples,
            n_polish=args.polish,
            seed=args.seed,
            warm_starts=warm_starts,
        )
        results[scheme] = res

    # Enforce Nesting Assertion for kinetic_coop from kinetic_two_open
    if "kinetic_coop" in results and "kinetic_two_open" in results:
        res_coop = results["kinetic_coop"]
        res_two = results["kinetic_two_open"]

        obj_j_coop = res_coop["joint_opt_val"]
        obj_j_two = res_two["joint_opt_val"]
        assert obj_j_coop <= obj_j_two + 1e-4, (
            f"Nesting violation (obj_natural_joint): coop ({obj_j_coop:.6f}) > two_open ({obj_j_two:.6f})"
        )

        obj_m_coop = res_coop["minimax_opt_val"]
        obj_m_two = res_two["minimax_opt_val"]
        assert obj_m_coop <= obj_m_two + 1e-4, (
            f"Nesting violation (obj_fourway_minimax): coop ({obj_m_coop:.6f}) > two_open ({obj_m_two:.6f})"
        )
        print("\n" + "=" * 80)
        print("[NESTING ASSERTION PASSED] kinetic_coop objective values strictly bounded by kinetic_two_open optimum.")
        print(f"  * obj_natural_joint:   coop {obj_j_coop:.6f} <= two_open {obj_j_two:.6f}")
        print(f"  * obj_fourway_minimax: coop {obj_m_coop:.6f} <= two_open {obj_m_two:.6f}")
        print("=" * 80 + "\n", flush=True)

    print("\n" + "=" * 80)
    print("PER-DATASET RESIDUALS SUMMARY TABLE (IN SIGMA)")
    print("=" * 80 + "\n")
    tbl = format_joint_results_table(results)
    print(tbl)
    print("\n")

    print("=" * 80)
    print("OPTIMAL PARAMETER VECTORS & BOUND FLAGS")
    print("=" * 80 + "\n")
    tbl_params = format_parameters_table(results)
    print(tbl_params)
    print("\n")

    # Definitive Verdict
    print("=" * 80)
    print("DEFINITIVE VERDICT")
    print("=" * 80)
    for scheme, r in results.items():
        spec = MODEL_SPECS[scheme]
        ev_joint = r["joint_opt"]
        ev_pair = r["pair_tradeoff"]
        ev_mini = r["minimax_opt"]

        min_viol = min(ev_joint["max_violation"], ev_pair["max_violation"], ev_mini["max_violation"])
        k_fit = len(spec.param_names)
        n_data = 9  # 3 taus + 2 independent amplitudes + EC50 + nH + recovery + po_max
        res_dof = n_data - k_fit
        print(f"\nScheme `{scheme}`:")
        print(f"  * n_data = {n_data}, k_fitted = {k_fit}, residual dof = {res_dof}")
        if res_dof == 0:
            print("  * NOTE: Fit is EXACTLY SATURATED (residual dof = 0). This is an existence result, not a test.")
        elif res_dof < 0:
            print(f"  * NOTE: Fit is OVERPARAMETERIZED (residual dof = {res_dof} < 0).")

        print(f"  * Four-Way Minimax Optimum residual: {ev_mini['max_violation'] * 3.0:.2f} sigma")
        print(f"  * Natural Joint Optimum max residual: {ev_joint['max_violation'] * 3.0:.2f} sigma")
        print(f"  * Two-Observable Pareto max residual: {ev_pair['max_violation'] * 3.0:.2f} sigma")

        if min_viol <= 1.0:
            print("  * Satisfies all four observables simultaneously within 3-sigma: YES")
            print(f"  * Optimal vector (minimax): {ev_mini['vector']}")
            print(f"  * Pinned parameters: {ev_mini['pinned']}")
            print(f"  * Effective free parameters (k_eff): {ev_mini['k_eff']}")
            print(f"  * Best minimax residual: {min_viol * 3.0:.2f} sigma")
        else:
            print("  * Satisfies all four observables simultaneously within 3-sigma: NO")
            print(f"  * Minimum minimax violation factor: {min_viol:.2f}x beyond 3-sigma threshold")
            print(f"  * At best compromise: max residual = {min_viol * 3.0:.2f} sigma")

        if scheme == "kinetic_coop":
            m_opt = spec.factory(ev_mini["vector"])
            print(f"  * Cooperativity factor c (minimax opt): {m_opt.c:.4f} (log10_c = {np.log10(m_opt.c):.4f})")
            print(f"  * koff: {m_opt.koff:.6f} ms^-1, koff2: {m_opt.koff2:.6f} ms^-1")
    print("=" * 80 + "\n")


if __name__ == "__main__":
    main()
