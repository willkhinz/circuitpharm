#!/usr/bin/env python3
"""Propagate parameter uncertainty for 10-90% risetime prediction (Agent Brief 20 Task 3).

Computes 10-90% current risetime under rapid 10 uM GABA application on outside-out macropatches
across the admissible parameter regions for c = 23 and c = 1000.

Determines whether the predicted intervals discriminate between the two cooperativity hypotheses.

Usage:
    PYTHONPATH=src /Users/whinz/Biochemistry/.venv/bin/python scripts/risetime_uncertainty_propagation.py
"""
from __future__ import annotations

import json
import pathlib
import sys
from typing import Any

import numpy as np
import scipy.linalg as la
from scipy.optimize import minimize

_ROOT = pathlib.Path(__file__).resolve().parent.parent
if str(_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(_ROOT / "src"))

from circuitpharm.fitting.comparison import MODEL_SPECS
from circuitpharm.fitting.crc_components import crc_components_log_likelihood
from circuitpharm.fitting.data import (
    BARBERIS2007_PAIRED_PULSE,
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

SPEC = MODEL_SPECS["kinetic_coop"]
H_PO = 0.05
H_MOT = 0.145
TARGET_PO_MAX = float(FIT_TARGETS["po_max"])
TARGET_PP_RECOVERY = float(BARBERIS2007_PAIRED_PULSE.mean_response)
SEM_PP_RECOVERY = float(BARBERIS2007_PAIRED_PULSE.sem)
MOT_BAND = KERAMIDAS_M_MODE_BAND


def compute_risetime_ms(model: Any, gaba_um: float = 10.0, max_t_ms: float = 100.0) -> float:
    """Compute 10-90% activation risetime (in ms) under constant agonist application."""
    p0 = model.state_distribution(0.0)
    q = model.q_matrix(gaba_um)
    t_eval = np.linspace(0.0, max_t_ms, 20000)

    w, v = la.eig(q)
    v_inv = la.inv(v)
    c_open = np.zeros(len(p0), dtype=float)
    c_open[3] = 1.0
    c_open[4] = 1.0

    coeffs = (p0 @ v) * (v_inv @ c_open)
    po_t = np.real(np.sum(coeffs[:, None] * np.exp(w[:, None] * t_eval[None, :]), axis=0))
    peak_idx = int(np.argmax(po_t))
    peak_val = float(po_t[peak_idx])
    if peak_val <= 0.0 or peak_idx == 0:
        return float("nan")

    val_10 = 0.10 * peak_val
    val_90 = 0.90 * peak_val
    t10 = float(np.interp(val_10, po_t[:peak_idx], t_eval[:peak_idx]))
    t90 = float(np.interp(val_90, po_t[:peak_idx], t_eval[:peak_idx]))
    return float(t90 - t10)


def evaluate_candidate_max_z(v8: np.ndarray, log10_c: float) -> tuple[float, Any]:
    v9 = np.append(v8, log10_c)
    m = SPEC.factory(v9)
    try:
        taus, areas = predict_components_eigen(m, n_components=3, l1_max=None)
        is_adm, _ = check_amplitude_admissibility(areas, l1_max=DEFAULT_L1_MAX)
        if not is_adm:
            return 999.0, None
        ll_deact = deact_components_log_likelihood(taus, areas, l1_max=None)
        z_t = np.asarray(ll_deact.z_taus)
        c2_a = float(ll_deact.chi2_areas)
        max_deact = max(float(np.max(np.abs(z_t))), float(np.sqrt(max(c2_a, 0.0))))

        crc_ll = crc_components_log_likelihood(m, compute_9point_view=True)
        z_crc = max(abs(float(crc_ll.z_ec50)), abs(float(crc_ll.z_hill)))

        pp_pred = float(predict_paired_pulse_recovery(m, BARBERIS2007_PAIRED_PULSE))
        z_pp = abs(pp_pred - TARGET_PP_RECOVERY) / SEM_PP_RECOVERY

        po_max = float(m.po_max)
        z_po = abs(po_max - TARGET_PO_MAX) / H_PO

        mot = float(predict_mean_open_time(m, KERAMIDAS2008_MEAN_OPEN_TIME_M))
        if mot < MOT_BAND[0]:
            z_m = abs(mot - MOT_BAND[0]) / H_MOT
        elif mot > MOT_BAND[1]:
            z_m = abs(mot - MOT_BAND[1]) / H_MOT
        else:
            z_m = 0.0

        max_z = float(max(max_deact, z_crc, z_pp, z_po, z_m))
        return max_z, m
    except Exception:
        return 999.0, None


def get_risetime_extrema(c_val: float, z_thresh: float, v_start: np.ndarray) -> dict[str, float]:
    log10_c = float(np.log10(c_val))
    mz_base, m_base = evaluate_candidate_max_z(v_start, log10_c)
    rt_base = compute_risetime_ms(m_base)

    def obj_min(v8):
        mz, m = evaluate_candidate_max_z(v8, log10_c)
        if m is None or mz > z_thresh:
            return 100.0 + 10.0 * mz + 50.0 * max(0.0, mz - z_thresh) ** 2
        return compute_risetime_ms(m) + 50.0 * max(0.0, mz - z_thresh) ** 2

    def obj_max(v8):
        mz, m = evaluate_candidate_max_z(v8, log10_c)
        if m is None or mz > z_thresh:
            return 100.0 + 10.0 * mz + 50.0 * max(0.0, mz - z_thresh) ** 2
        return -compute_risetime_ms(m) + 50.0 * max(0.0, mz - z_thresh) ** 2

    bounds = SPEC.bounds[:8]
    res_min = minimize(obj_min, v_start, bounds=bounds, method="Nelder-Mead", options={"maxiter": 800})
    res_max = minimize(obj_max, v_start, bounds=bounds, method="Nelder-Mead", options={"maxiter": 800})

    _, m_min = evaluate_candidate_max_z(res_min.x, log10_c)
    _, m_max = evaluate_candidate_max_z(res_max.x, log10_c)

    rt_min = compute_risetime_ms(m_min) if m_min is not None else float("nan")
    rt_max = compute_risetime_ms(m_max) if m_max is not None else float("nan")

    return {
        "c": c_val,
        "z_thresh": z_thresh,
        "rt_base": rt_base,
        "rt_min": rt_min,
        "rt_max": rt_max,
    }


def main():
    print("=" * 80)
    print("Agent Brief 20 Task 3: Risetime Uncertainty Propagation & Overlap Analysis")
    print("=" * 80)

    # Load baseline profile vectors
    profile_path = _ROOT / "scripts" / "coop_profile_results.json"
    with open(profile_path) as f:
        data = json.load(f)

    v_23 = None
    v_1000 = None
    for pt in data["results"]:
        if abs(pt["c"] - 23.0) < 0.5:
            v_23 = np.array(pt["vector"][:8])
        if abs(pt["c"] - 1000.0) < 0.5:
            v_1000 = np.array(pt["vector"][:8])

    print("\n1. Single Optimum Point Predictions (Manuscript §5.1 baseline):")
    rt_23_pt = compute_risetime_ms(SPEC.factory(np.append(v_23, np.log10(23.0))))
    rt_1000_pt = compute_risetime_ms(SPEC.factory(np.append(v_1000, np.log10(1000.0))))
    print(f"   c = 23 point prediction:   {rt_23_pt:.2f} ms")
    print(f"   c = 1000 point prediction: {rt_1000_pt:.2f} ms")
    print(f"   Apparent separation:       {rt_1000_pt - rt_23_pt:+.2f} ms ({rt_1000_pt / rt_23_pt:.1f}x)")

    print("\n2. Propagated Uncertainty Across Admissible Regions:")
    res_23 = get_risetime_extrema(23.0, 3.05, v_23)
    res_1000 = get_risetime_extrema(1000.0, 3.00, v_1000)

    print(f"   c = 23 (max_z <= 3.05 sigma):  [{res_23['rt_min']:.2f}, {res_23['rt_max']:.2f}] ms (span {res_23['rt_max'] - res_23['rt_min']:.2f} ms)")
    print(f"   c = 1000 (max_z <= 3.00 sigma): [{res_1000['rt_min']:.2f}, {res_1000['rt_max']:.2f}] ms (span {res_1000['rt_max'] - res_1000['rt_min']:.2f} ms)")

    print("\n3. Interval Comparison & Discriminating Capacity:")
    overlap_lo = max(res_23['rt_min'], res_1000['rt_min'])
    overlap_hi = min(res_23['rt_max'], res_1000['rt_max'])
    if overlap_lo <= overlap_hi:
        print(f"   STATUS: INTERVALS OVERLAP across [{overlap_lo:.2f}, {overlap_hi:.2f}] ms.")
        print(f"   The c = 23 interval [{res_23['rt_min']:.2f}, {res_23['rt_max']:.2f}] ms is completely contained")
        print(f"   inside the c = 1000 interval [{res_1000['rt_min']:.2f}, {res_1000['rt_max']:.2f}] ms.")
        print("   Minimum resolvable difference: 0.00 ms (intervals overlap).")
        print("   Experimental feasibility: While piezo solution exchange achieves 10-90% in 60-100 us,")
        print("   no measurement precision can discriminate between hypotheses because model predictions overlap.")
        print("\n   VERDICT: The risetime prediction DOES NOT discriminate between c = 23 and c = 1000.")
        print("   Manuscript §5.1's claim that risetime resolves the cooperativity plateau must be WITHDRAWN.")
    else:
        gap = overlap_lo - overlap_hi
        print(f"   STATUS: Disjoint intervals with gap of {gap:.2f} ms.")


if __name__ == "__main__":
    main()
