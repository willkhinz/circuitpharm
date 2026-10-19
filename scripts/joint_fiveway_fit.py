#!/usr/bin/env python3
"""Symmetric Joint Five-Way Fit Evaluation (Agent Brief 15 Task 2).

Evaluates whether any parameter vector of kinetic_jw95_free_rates, kinetic_two_open,
or kinetic_coop can simultaneously satisfy all FIVE observables within 3 sigma,
with EVERY parametric dataset scored against its published parameters (Weighting iv extended):

1. Deactivation: 6 published components at published errors (tau: +/-0.3, +/-4.6, +/-14.9 ms;
   areas: rank-2 covariance from +/-0.04, +/-0.02, +/-0.03) under L1_MAX <= 1.5. (5 data terms)
2. Peak CRC: EC50 = 11.6 +/- 0.9 uM, nH = 2.2 +/- 0.4 over 1-10 uM (Jahn et al. 1997). (2 data terms)
3. Paired-pulse recovery: 0.33 +/- 0.03 (Barberis et al. 2007). (1 data term)
4. P_o,max gating anchor: band [0.70, 0.80] (target 0.75 +/- 0.05). (1 data term)
5. Mean open time: 2.68 ms, Table III band [2.54, 2.83] ms (Keramidas & Harrison 2008). (1 data term)

Total data degrees of freedom: n_data = 10.
Model parameters k_fitted:
- kinetic_jw95_free_rates: k = 6  => res_dof = 10 - 6 = +4
- kinetic_two_open:        k = 8  => res_dof = 10 - 8 = +2
- kinetic_coop:            k = 9  => res_dof = 10 - 9 = +1 (TESTED FIT)

Nesting property:
kinetic_coop strictly nests kinetic_two_open at c = 1.0. Warm-start is enforced from
kinetic_two_open's optimum via nesting_embedding, and a strict assertion prevents
reporting any nesting violation (coop <= two_open + 1e-4).

Usage:
    PYTHONPATH=src /Users/whinz/Biochemistry/.venv/bin/python scripts/joint_fiveway_fit.py
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
    KERAMIDAS2008_MEAN_OPEN_TIME_M,
    KERAMIDAS_M_MODE_BAND,
)
from circuitpharm.fitting.deact_components import (
    DEFAULT_L1_MAX,
    check_amplitude_admissibility,
    deact_components_log_likelihood,
    predict_components_eigen,
)
from circuitpharm.fitting.likelihood import (
    predict_mean_open_time,
    predict_paired_pulse_recovery,
)
from circuitpharm.gabaa_kinetics import FIT_RANGES, FIT_TARGETS

TARGET_PP_RECOVERY = float(BARBERIS2007_PAIRED_PULSE.mean_response)
SEM_PP_RECOVERY = float(BARBERIS2007_PAIRED_PULSE.sem)
TARGET_PO_MAX = float(FIT_TARGETS["po_max"])
SIGMA_PO_MAX = 0.05  # conventional calibration SEM (half-width of [0.70, 0.80])

TARGET_MOT_MS = float(KERAMIDAS2008_MEAN_OPEN_TIME_M.mean_open_time_ms)  # 2.68 ms
MOT_BAND_LO, MOT_BAND_HI = KERAMIDAS_M_MODE_BAND                        # (2.54, 2.83) ms
SIGMA_MOT_MS = (MOT_BAND_HI - MOT_BAND_LO) / 2.0                        # 0.145 ms


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
    mot_mode: str = "band",
) -> dict[str, Any]:
    """Score one candidate vector across all five observables symmetrically without invented weighting."""
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

    # 5. Mean open time (Keramidas 2008 Table III, M-mode target 2.68 ms, band [2.54, 2.83] ms)
    try:
        mot_ms = float(predict_mean_open_time(m, KERAMIDAS2008_MEAN_OPEN_TIME_M))
        mot_in_band = bool(MOT_BAND_LO <= mot_ms <= MOT_BAND_HI)
        # Band z-score: 0 if in band; distance to nearest boundary / sigma if outside
        if mot_ms < MOT_BAND_LO:
            z_mot_band = float((mot_ms - MOT_BAND_LO) / SIGMA_MOT_MS)
        elif mot_ms > MOT_BAND_HI:
            z_mot_band = float((mot_ms - MOT_BAND_HI) / SIGMA_MOT_MS)
        else:
            z_mot_band = 0.0
        z_mot_midpoint = float((mot_ms - TARGET_MOT_MS) / SIGMA_MOT_MS)
    except Exception:
        mot_ms = float("nan")
        mot_in_band = False
        z_mot_band = 999.0
        z_mot_midpoint = 999.0

    z_mot = z_mot_band if mot_mode == "band" else z_mot_midpoint

    pinned = check_at_bound(v, spec.bounds, spec.param_names)
    k_eff = len(spec.param_names) - len(pinned)

    # Normalized 3-sigma violation across all five observables
    # Satisfying all 5 within 3 sigma requires max_violation <= 1.000
    terms = [
        abs(z_taus[0]) / 3.0,
        abs(z_taus[1]) / 3.0,
        abs(z_taus[2]) / 3.0,
        np.sqrt(max(chi2_areas, 0.0)) / 3.0,
        abs(z_ec50) / 3.0,
        abs(z_hill) / 3.0,
        (abs(z_pp) / 3.0) if np.isfinite(z_pp) else 0.0,
        abs(z_po) / 3.0,
        abs(z_mot) / 3.0,
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
        "mot_ms": mot_ms,
        "mot_in_band": mot_in_band,
        "z_mot_band": z_mot_band,
        "z_mot_midpoint": z_mot_midpoint,
        "z_mot": z_mot,
        "pinned": pinned,
        "k_eff": k_eff,
        "max_violation": max_violation,
    }


def fit_joint_fiveway(
    spec: ModelSpec,
    n_samples: int = 2000,
    n_polish: int = 12,
    seed: int = 42,
    l1_max: float = DEFAULT_L1_MAX,
    warm_starts: dict[str, np.ndarray] | None = None,
    mot_mode: str = "band",
) -> dict[str, Any]:
    """Execute multi-objective and feasibility evaluation for a model scheme across all five observables."""
    t0 = time.perf_counter()
    lo = np.array([b[0] for b in spec.bounds], dtype=float)
    hi = np.array([b[1] for b in spec.bounds], dtype=float)
    d = len(spec.bounds)

    sampler = qmc.LatinHypercube(d=d, seed=seed)
    pts = lo + sampler.random(n_samples) * (hi - lo)

    def mot_penalty(m):
        try:
            mot = float(predict_mean_open_time(m, KERAMIDAS2008_MEAN_OPEN_TIME_M))
            if mot_mode == "band":
                if mot < MOT_BAND_LO:
                    return float(((mot - MOT_BAND_LO) / SIGMA_MOT_MS) ** 2)
                elif mot > MOT_BAND_HI:
                    return float(((mot - MOT_BAND_HI) / SIGMA_MOT_MS) ** 2)
                else:
                    return 0.0
            else:
                return float(((mot - TARGET_MOT_MS) / SIGMA_MOT_MS) ** 2)
        except Exception:
            return 1e6

    def mot_z(m):
        try:
            mot = float(predict_mean_open_time(m, KERAMIDAS2008_MEAN_OPEN_TIME_M))
            if mot_mode == "band":
                if mot < MOT_BAND_LO:
                    return abs(mot - MOT_BAND_LO) / SIGMA_MOT_MS
                elif mot > MOT_BAND_HI:
                    return abs(mot - MOT_BAND_HI) / SIGMA_MOT_MS
                else:
                    return 0.0
            else:
                return abs(mot - TARGET_MOT_MS) / SIGMA_MOT_MS
        except Exception:
            return 999.0

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

    # 2. Single-objective: Sourced Peak CRC alone
    def obj_crc(v):
        try:
            m = spec.factory(v)
            p = extract_crc_parameters(m)
            z_e = (p.ec50_um - JAHN_EC50_UM) / JAHN_EC50_SEM
            z_h = (p.hill_measured - JAHN_HILL) / JAHN_HILL_SEM
            return float(z_e ** 2 + z_h ** 2)
        except Exception:
            return 1e6

    # 3. Natural Published-Error Joint Objective (all 5 observables, n_data = 10)
    # Deact (5) + CRC (2) + Recovery (1) + Po_max (1) + MOT (1) = 10
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

            c2_mot = mot_penalty(m)

            return float(c2_deact + c2_crc + c2_pp + c2_po + c2_mot)
        except Exception:
            return 1e8

    # 4. True Five-Way Minimax Objective
    def obj_fiveway_minimax(v):
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

            z_m = mot_z(m)

            return float(max(max_z_deact, max_z_crc, z_pp, z_po, z_m))
        except Exception:
            return 1e6

    print(f"[{spec.name}] Optimizing Deactivation alone...", flush=True)
    d_candidates = [min(pts, key=obj_deact)]
    if warm_starts and "deact" in warm_starts:
        starts = warm_starts["deact"] if isinstance(warm_starts["deact"], list) else [warm_starts["deact"]]
        d_candidates.extend(starts)
    best_d_start = min(d_candidates, key=obj_deact)
    opt_deact = minimize(obj_deact, best_d_start, bounds=spec.bounds, method="L-BFGS-B")

    print(f"[{spec.name}] Optimizing Sourced Peak CRC alone...", flush=True)
    c_candidates = [min(pts, key=obj_crc)]
    if warm_starts and "crc" in warm_starts:
        starts = warm_starts["crc"] if isinstance(warm_starts["crc"], list) else [warm_starts["crc"]]
        c_candidates.extend(starts)
    best_c_start = min(c_candidates, key=obj_crc)
    opt_crc = minimize(obj_crc, best_c_start, bounds=spec.bounds, method="L-BFGS-B")

    print(f"[{spec.name}] Optimizing Natural Published-Error Joint Objective (5-way)...", flush=True)
    joint_starts = [opt_deact.x, opt_crc.x] + sorted(pts, key=obj_natural_joint)[:n_polish]
    if warm_starts and "joint" in warm_starts:
        starts = warm_starts["joint"] if isinstance(warm_starts["joint"], list) else [warm_starts["joint"]]
        joint_starts = starts + joint_starts
    joint_optima = []
    for s in joint_starts:
        res_j = minimize(obj_natural_joint, s, bounds=spec.bounds, method="Nelder-Mead", options={"maxiter": 1500})
        joint_optima.append((obj_natural_joint(res_j.x), res_j.x))
    joint_optima.sort(key=lambda x: x[0])
    best_joint_vec = joint_optima[0][1]
    best_joint_val = joint_optima[0][0]

    print(f"[{spec.name}] Optimizing True Five-Way Minimax Objective...", flush=True)
    minimax_starts = [best_joint_vec, opt_deact.x, opt_crc.x] + sorted(pts, key=obj_fiveway_minimax)[:n_polish]
    if warm_starts and "minimax" in warm_starts:
        starts = warm_starts["minimax"] if isinstance(warm_starts["minimax"], list) else [warm_starts["minimax"]]
        minimax_starts = starts + minimax_starts
    minimax_optima = []
    for s in minimax_starts:
        res_m = minimize(obj_fiveway_minimax, s, bounds=spec.bounds, method="Nelder-Mead", options={"maxiter": 1500})
        minimax_optima.append((obj_fiveway_minimax(res_m.x), res_m.x))
    minimax_optima.sort(key=lambda x: x[0])
    best_minimax_vec = minimax_optima[0][1]
    best_minimax_val = minimax_optima[0][0]

    eval_deact_opt = evaluate_candidate(spec, opt_deact.x, l1_max=l1_max, compute_recovery=True, mot_mode=mot_mode)
    eval_crc_opt = evaluate_candidate(spec, opt_crc.x, l1_max=l1_max, compute_recovery=True, mot_mode=mot_mode)
    eval_joint_opt = evaluate_candidate(spec, best_joint_vec, l1_max=l1_max, compute_recovery=True, mot_mode=mot_mode)
    eval_minimax_opt = evaluate_candidate(spec, best_minimax_vec, l1_max=l1_max, compute_recovery=True, mot_mode=mot_mode)

    runtime_s = time.perf_counter() - t0

    return {
        "scheme": spec.name,
        "runtime_s": runtime_s,
        "deact_opt": eval_deact_opt,
        "crc_opt": eval_crc_opt,
        "joint_opt": eval_joint_opt,
        "minimax_opt": eval_minimax_opt,
        "joint_opt_val": best_joint_val,
        "minimax_opt_val": best_minimax_val,
    }


def format_joint_results_table(results: dict[str, dict[str, Any]]) -> str:
    """Format markdown summary table of per-dataset residuals in sigma."""
    lines = [
        "| Scheme | Evaluation Profile | $z_{\\tau_1}$ | $z_{\\tau_2}$ | $z_{\\tau_3}$ | $z_{a_1}$ | $z_{a_2}$ | $z_{a_3}$ | $\\sqrt{\\chi^2_a}$ | $z_{\\text{EC50}}$ | $z_{n_H}$ | $z_{pp}$ | $P_{o,max}$ | $\\tau_O$ (ms) | $z_{\\text{MOT}}$ (band) | $z_{\\text{MOT}}$ (mid) | Max Viol / $3\\sigma$ | Satisfies All 5 $\\le 3\\sigma$? |",
        "| :--- | :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |",
    ]

    for scheme, r in results.items():
        profiles = [
            ("Deact-Alone Optimum", r["deact_opt"]),
            ("CRC-Alone Optimum", r["crc_opt"]),
            ("Natural Joint Optimum", r["joint_opt"]),
            ("Five-Way Minimax Optimum", r["minimax_opt"]),
        ]
        for prof_name, ev in profiles:
            zt = ev["z_taus"]
            za = ev["z_areas"]
            sq_c2a = np.sqrt(max(ev["chi2_areas"], 0.0))
            zec50 = ev["z_ec50"]
            znh = ev["z_hill"]
            zpp = ev["z_pp"]
            zpp_str = f"{zpp:+.2f}$\\sigma$" if np.isfinite(zpp) else "N/A"
            po_str = f"{ev['po_max']:.3f} ({ev['z_po']:+.2f}$\\sigma$)"
            mot_val_str = f"{ev['mot_ms']:.2f}" if np.isfinite(ev["mot_ms"]) else "N/A"
            mot_band_str = f"{ev['z_mot_band']:+.2f}$\\sigma$" if np.isfinite(ev["z_mot_band"]) else "N/A"
            mot_mid_str = f"{ev['z_mot_midpoint']:+.2f}$\\sigma$" if np.isfinite(ev["z_mot_midpoint"]) else "N/A"
            viol = ev["max_violation"]
            sat_str = "**YES**" if viol <= 1.0 else "**NO**"

            lines.append(
                f"| `{scheme}` | {prof_name} | {zt[0]:+.2f}$\\sigma$ | {zt[1]:+.2f}$\\sigma$ | {zt[2]:+.2f}$\\sigma$ | "
                f"{za[0]:+.2f}$\\sigma$ | {za[1]:+.2f}$\\sigma$ | {za[2]:+.2f}$\\sigma$ | {sq_c2a:.2f}$\\sigma$ | "
                f"{zec50:+.2f}$\\sigma$ | {znh:+.2f}$\\sigma$ | "
                f"{zpp_str} | {po_str} | {mot_val_str} | {mot_band_str} | {mot_mid_str} | **{viol:.2f}** | {sat_str} |"
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
            ("Natural Joint Optimum", "joint_opt"),
            ("Five-Way Minimax Optimum", "minimax_opt"),
        ]:
            ev = r[key]
            v_str = ", ".join(f"{name}={val:.3f}" for name, val in zip(spec.param_names, ev["vector"]))
            pinned_str = "; ".join(ev["pinned"]) if ev["pinned"] else "None (all interior)"
            lines.append(
                f"| `{scheme}` | {prof_name} | `{v_str}` | {pinned_str} | {len(spec.param_names)} | {ev['k_eff']} |"
            )
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(description="Symmetric Joint Five-Way Fit Evaluation (Brief 15 Task 2).")
    parser.add_argument("--samples", type=int, default=2000, help="Number of LHS samples (default: 2000).")
    parser.add_argument("--polish", type=int, default=12, help="Number of polished starts (default: 12).")
    parser.add_argument("--seed", type=int, default=42, help="Random seed (default: 42).")
    parser.add_argument(
        "--schemes",
        nargs="+",
        default=["kinetic_two_open", "kinetic_coop", "kinetic_three_open"],
        help="Schemes to evaluate (default: kinetic_two_open kinetic_coop kinetic_three_open).",
    )
    parser.add_argument(
        "--mot-mode",
        choices=["band", "midpoint"],
        default="band",
        help="How to score Keramidas 2008: 'band' (zero penalty within [2.54, 2.83]) or 'midpoint' (scored vs 2.68).",
    )
    args = parser.parse_args()

    print("=" * 80)
    print("SYMMETRIC JOINT FIVE-WAY FIT EVALUATION (Agent Brief 15)")
    print("=" * 80)
    print(f"Schemes: {args.schemes}")
    print(f"Budget per scheme: {args.samples} LHS samples, top {args.polish} polished starts")
    print(f"Target 1: Deactivation components (6 published parameters at published errors, L1 <= {DEFAULT_L1_MAX})")
    print(f"Target 2: Sourced Peak CRC (EC50 = {JAHN_EC50_UM} +/- {JAHN_EC50_SEM} uM, nH = {JAHN_HILL} +/- {JAHN_HILL_SEM} over 1-10 uM)")
    print(f"Target 3: Paired-pulse recovery ({TARGET_PP_RECOVERY} +/- {SEM_PP_RECOVERY} at 100 ms gap)")
    print(f"Target 4: P_o,max gating anchor (band {FIT_RANGES['po_max']}, target {TARGET_PO_MAX} +/- {SIGMA_PO_MAX})")
    print(f"Target 5: Mean open time (target {TARGET_MOT_MS} ms, Table III band [{MOT_BAND_LO}, {MOT_BAND_HI}] ms, mode={args.mot_mode})")
    print("NO INVENTED WEIGHTING: reporting per-dataset residuals in sigma.")
    print("-" * 80, flush=True)

    results = {}
    for scheme in args.schemes:
        spec = MODEL_SPECS[scheme]
        warm_starts = None

        if scheme == "kinetic_two_open":
            v_two_b14_m = np.array([0.459, 0.257, 1.773, -1.795, -0.718, 0.228, -1.200, -1.235])
            v_two_b14_j = np.array([0.454, -0.042, 1.423, -1.832, -1.015, 1.531, -1.100, -1.476])
            v_two_5way = np.array([0.4680, 0.2782, 1.9409, -1.8517, -0.6366, 0.2463, -0.9550, -1.3674])
            warm_starts = {
                "deact": [v_two_b14_m],
                "joint": [v_two_b14_j, v_two_5way],
                "minimax": [v_two_b14_m, v_two_5way],
            }
        elif scheme == "kinetic_coop" and "kinetic_two_open" in results:
            spec_sub = MODEL_SPECS["kinetic_two_open"]
            res_sub = results["kinetic_two_open"]
            v_two_opt_j = res_sub["joint_opt"]["vector"]
            v_two_opt_m = res_sub["minimax_opt"]["vector"]
            v_two_opt_d = res_sub["deact_opt"]["vector"]
            v_two_opt_c = res_sub["crc_opt"]["vector"]

            v_coop_b14_m = np.array([2.530, -0.073, 1.271, -1.858, 1.878, 1.413, -1.146, -1.460, 2.912])
            v_coop_b14_j = np.array([2.438, -0.071, 1.381, -1.825, 1.937, 1.333, -1.082, -1.483, 2.981])
            v_coop_5way = np.array([2.7420, -0.0828, 1.3951, -1.8619, 2.2637, 0.2393, -1.1410, -1.5737, 3.3895])

            def _clip(v_in, b_in):
                lo = np.array([b[0] for b in b_in], dtype=float)
                hi = np.array([b[1] for b in b_in], dtype=float)
                return np.clip(v_in, lo, hi)

            warm_starts = {
                "deact": [nesting_embedding(spec_sub, spec, _clip(v_two_opt_d, spec_sub.bounds)), v_coop_b14_m],
                "crc": [nesting_embedding(spec_sub, spec, _clip(v_two_opt_c, spec_sub.bounds)), v_coop_b14_m],
                "joint": [nesting_embedding(spec_sub, spec, _clip(v_two_opt_j, spec_sub.bounds)), v_coop_b14_j, v_coop_5way],
                "minimax": [nesting_embedding(spec_sub, spec, _clip(v_two_opt_m, spec_sub.bounds)), v_coop_b14_m, v_coop_5way],
            }
            print(f"[{scheme}] Warm-starting from kinetic_two_open optima via nesting_embedding (c = 1.0) and cooperativity islands...", flush=True)
        elif scheme == "kinetic_three_open" and "kinetic_coop" in results:
            spec_sub = MODEL_SPECS["kinetic_coop"]
            res_sub = results["kinetic_coop"]
            v_coop_opt_j = res_sub["joint_opt"]["vector"]
            v_coop_opt_m = res_sub["minimax_opt"]["vector"]
            v_coop_opt_d = res_sub["deact_opt"]["vector"]
            v_coop_opt_c = res_sub["crc_opt"]["vector"]

            def _clip(v_in, b_in):
                lo = np.array([b[0] for b in b_in], dtype=float)
                hi = np.array([b[1] for b in b_in], dtype=float)
                return np.clip(v_in, lo, hi)

            warm_starts = {
                "deact": [nesting_embedding(spec_sub, spec, _clip(v_coop_opt_d, spec_sub.bounds))],
                "crc": [nesting_embedding(spec_sub, spec, _clip(v_coop_opt_c, spec_sub.bounds))],
                "joint": [nesting_embedding(spec_sub, spec, _clip(v_coop_opt_j, spec_sub.bounds))],
                "minimax": [nesting_embedding(spec_sub, spec, _clip(v_coop_opt_m, spec_sub.bounds))],
            }
            print(f"[{scheme}] Warm-starting from kinetic_coop optima via nesting_embedding (delta2 -> 0)...", flush=True)

        res = fit_joint_fiveway(
            spec,
            n_samples=args.samples,
            n_polish=args.polish,
            seed=args.seed,
            warm_starts=warm_starts,
            mot_mode=args.mot_mode,
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
            f"Nesting violation (obj_fiveway_minimax): coop ({obj_m_coop:.6f}) > two_open ({obj_m_two:.6f})"
        )
        print("\n" + "=" * 80)
        print("[NESTING ASSERTION PASSED] kinetic_coop objective values strictly bounded by kinetic_two_open optimum.")
        print(f"  * obj_natural_joint:   coop {obj_j_coop:.6f} <= two_open {obj_j_two:.6f}")
        print(f"  * obj_fiveway_minimax: coop {obj_m_coop:.6f} <= two_open {obj_m_two:.6f}")
        print("=" * 80 + "\n", flush=True)

    # Enforce Nesting Assertion for kinetic_three_open from kinetic_coop
    if "kinetic_three_open" in results and "kinetic_coop" in results:
        res_three = results["kinetic_three_open"]
        res_coop = results["kinetic_coop"]

        obj_j_three = res_three["joint_opt_val"]
        obj_j_coop = res_coop["joint_opt_val"]
        assert obj_j_three <= obj_j_coop + 1e-4, (
            f"Nesting violation (obj_natural_joint): three_open ({obj_j_three:.6f}) > coop ({obj_j_coop:.6f})"
        )

        obj_m_three = res_three["minimax_opt_val"]
        obj_m_coop = res_coop["minimax_opt_val"]
        assert obj_m_three <= obj_m_coop + 1e-4, (
            f"Nesting violation (obj_fiveway_minimax): three_open ({obj_m_three:.6f}) > coop ({obj_m_coop:.6f})"
        )
        print("\n" + "=" * 80)
        print("[NESTING ASSERTION PASSED] kinetic_three_open objective values strictly bounded by kinetic_coop optimum.")
        print(f"  * obj_natural_joint:   three_open {obj_j_three:.6f} <= coop {obj_j_coop:.6f}")
        print(f"  * obj_fiveway_minimax: three_open {obj_m_three:.6f} <= coop {obj_m_coop:.6f}")
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
        ev_mini = r["minimax_opt"]

        min_viol = min(ev_joint["max_violation"], ev_mini["max_violation"])
        k_fit = len(spec.param_names)
        n_data = 10  # 3 taus + 2 independent amplitudes + EC50 + nH + recovery + po_max + mean_open_time
        res_dof = n_data - k_fit
        print(f"\nScheme `{scheme}`:")
        print(f"  * n_data = {n_data}, k_fitted = {k_fit}, residual dof = {res_dof}")
        if res_dof >= 1:
            print(f"  * res_dof = +{res_dof} >= +1: Genuine tested hypothesis fit!")
        elif res_dof == 0:
            print("  * NOTE: Fit is EXACTLY SATURATED (residual dof = 0). This is an existence result, not a test.")
        elif res_dof < 0:
            print(f"  * NOTE: Fit is OVERPARAMETERIZED (residual dof = {res_dof} < 0). Fits with residual dof <= 0 cannot serve as statistical tests; they evaluate parameter saturation.")

        print(f"  * Five-Way Minimax Optimum worst residual: {ev_mini['max_violation'] * 3.0:.2f} sigma")
        print(f"  * Natural Joint Optimum worst residual: {ev_joint['max_violation'] * 3.0:.2f} sigma")

        if min_viol <= 1.0:
            print("  * Satisfies all five observables simultaneously within 3-sigma: YES")
            print(f"  * Optimal vector (minimax): {ev_mini['vector']}")
            print(f"  * Pinned parameters: {ev_mini['pinned']}")
            print(f"  * Effective free parameters (k_eff): {ev_mini['k_eff']}")
            print(f"  * Worst residual: {min_viol * 3.0:.2f} sigma")
        else:
            print("  * Satisfies all five observables simultaneously within 3-sigma: NO")
            print(f"  * Minimum minimax violation factor: {min_viol:.2f}x beyond 3-sigma threshold")
            print(f"  * At best compromise: max residual = {min_viol * 3.0:.2f} sigma")

        if scheme in ("kinetic_coop", "kinetic_three_open"):
            m_opt = spec.factory(ev_mini["vector"])
            print(f"  * Cooperativity factor c (minimax opt): {m_opt.c:.4f} (log10_c = {np.log10(m_opt.c):.4f})")
            print(f"  * koff: {m_opt.koff:.6f} ms^-1, koff2: {m_opt.koff2:.6f} ms^-1")
            print(f"  * Mean open time: {ev_mini['mot_ms']:.3f} ms (Keramidas Table III M-mode band: [{MOT_BAND_LO}, {MOT_BAND_HI}] ms)")
            print(f"  * Hill slope n_H: {ev_mini['hill_meas']:.3f} (z_nH = {ev_mini['z_hill']:+.2f} sigma vs target 2.2 +/- 0.4)")
            if scheme == "kinetic_three_open":
                print(f"  * delta2: {m_opt.delta2:.6f} ms^-1, gamma2: {m_opt.gamma2:.6f} ms^-1")
    print("=" * 80 + "\n")


if __name__ == "__main__":
    main()
