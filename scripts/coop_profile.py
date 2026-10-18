#!/usr/bin/env python3
"""Profile binding cooperativity c in kinetic_coop under five-way minimax objective.

Agent Brief 16 Task 1 & Task 2:
Investigates how the five-way minimax varies with cooperativity c held fixed over
a grid spanning 1 to 10,000 (log-spaced, >= 20 points, including 5, 23, 50, 100, 1000, 10000).

Warm-starts each grid point from neighbouring solutions in both directions to eliminate
optimiser hysteresis and spurious non-monotonicity.

Usage:
    PYTHONPATH=src /Users/whinz/Biochemistry/.venv/bin/python scripts/coop_profile.py
"""
from __future__ import annotations

import argparse
import json
import pathlib
import sys
import time
from typing import Any

import numpy as np
from scipy.optimize import minimize

# Ensure src and root are on sys.path
_ROOT = pathlib.Path(__file__).resolve().parent.parent
if str(_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(_ROOT / "src"))
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

from circuitpharm.fitting.comparison import MODEL_SPECS
from circuitpharm.fitting.crc_components import extract_crc_parameters
from circuitpharm.fitting.data import (
    BARBERIS2007_PAIRED_PULSE,
    JAHN_EC50_SEM,
    JAHN_EC50_UM,
    JAHN_HILL,
    JAHN_HILL_SEM,
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
from circuitpharm.gabaa_kinetics import FIT_TARGETS
from scripts.joint_fiveway_fit import evaluate_candidate

TARGET_PP_RECOVERY = float(BARBERIS2007_PAIRED_PULSE.mean_response)
SEM_PP_RECOVERY = float(BARBERIS2007_PAIRED_PULSE.sem)
TARGET_PO_MAX = float(FIT_TARGETS["po_max"])
SIGMA_PO_MAX = 0.05

TARGET_MOT_MS = float(KERAMIDAS2008_MEAN_OPEN_TIME_M.mean_open_time_ms)
MOT_BAND_LO, MOT_BAND_HI = KERAMIDAS_M_MODE_BAND
SIGMA_MOT_MS = (MOT_BAND_HI - MOT_BAND_LO) / 2.0

SPEC = MODEL_SPECS["kinetic_coop"]
BOUNDS_8 = SPEC.bounds[:8]

# Fitted minimax starting point at c = 1000 from Brief 15
V_1000_8 = np.array([2.603, -0.036, 1.458, -1.851, 1.933, 0.068, -0.989, -1.428])

# Two-open minimax starting point at c = 1 from Brief 15
V_1_8 = np.array([0.523, 0.292, 1.999, -1.882, -0.564, 0.159, -1.142, -1.484])

# Defined 23-point log-spaced grid spanning 1 to 10,000, explicitly containing key anchors
GRID_C = [
    1.0, 1.58, 2.51, 3.98, 5.0, 7.08, 10.0, 15.14, 23.0, 35.0, 50.0,
    75.0, 100.0, 158.49, 251.19, 398.11, 630.96, 1000.0, 1584.89,
    2511.89, 3981.07, 6309.57, 10000.0
]


def make_objective(log10_c: float, l1_max: float = DEFAULT_L1_MAX):
    """Create five-way minimax objective for the 8 free parameters at fixed log10_c."""
    def obj(v8: np.ndarray) -> float:
        v9 = np.append(v8, log10_c)
        try:
            m = SPEC.factory(v9)
            taus, areas = predict_components_eigen(m, n_components=3, l1_max=None)
            is_adm, c_rat = check_amplitude_admissibility(areas, l1_max=l1_max)
            if not is_adm:
                return 1e4 + c_rat * 10
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

            mot = float(predict_mean_open_time(m, KERAMIDAS2008_MEAN_OPEN_TIME_M))
            if mot < MOT_BAND_LO:
                z_m = abs(mot - MOT_BAND_LO) / SIGMA_MOT_MS
            elif mot > MOT_BAND_HI:
                z_m = abs(mot - MOT_BAND_HI) / SIGMA_MOT_MS
            else:
                z_m = 0.0

            return float(max(max_z_deact, max_z_crc, z_pp, z_po, z_m))
        except Exception:
            return 1e6
    return obj


def compute_peak_p_ar(m: Any) -> float:
    """Compute peak equilibrium occupancy of singly-bound state AR across concentrations."""
    concs = np.logspace(-3, 6, 300)
    return float(max(m.state_distribution(g)[1] for g in concs))


def identify_binding_observable(eval_res: dict[str, Any]) -> str:
    """Identify which of the 5 observables binds at the minimax solution."""
    z_t = np.abs(eval_res["z_taus"])
    z_a = np.sqrt(max(eval_res["chi2_areas"], 0.0))
    z_deact = max(float(np.max(z_t)), float(z_a))

    z_ec50 = abs(eval_res["z_ec50"])
    z_hill = abs(eval_res["z_hill"])
    z_crc = max(z_ec50, z_hill)

    z_pp = abs(eval_res["z_pp"])
    z_po = abs(eval_res["z_po"])
    z_mot = abs(eval_res["z_mot_band"])

    components = [
        (z_deact, f"Deact (tau_max={np.max(z_t):.2f}s, area={z_a:.2f}s)"),
        (z_crc, f"CRC (EC50={z_ec50:.2f}s, nH={z_hill:.2f}s)"),
        (z_pp, f"Recovery ({z_pp:.2f}s)"),
        (z_po, f"Po_max ({z_po:.2f}s)"),
        (z_mot, f"MOT ({z_mot:.2f}s)"),
    ]
    components.sort(key=lambda x: x[0], reverse=True)
    top_val, top_name = components[0]
    
    # Check for near-ties within 0.05 sigma
    binding = [name.split()[0] for val, name in components if top_val - val < 0.05]
    return "/".join(binding)


def run_profile(
    grid_c: list[float] | None = None,
    maxiter: int = 600,
    xatol: float = 1e-3,
    fatol: float = 1e-4,
) -> dict[str, Any]:
    """Execute bidirectional warm-started profile over c."""
    if grid_c is None:
        grid_c = GRID_C

    log10_grid = [float(np.log10(c)) for c in grid_c]
    n_points = len(grid_c)
    idx_1000 = grid_c.index(1000.0)

    print(f"=== Cooperativity Profile over c ===")
    print(f"Grid size: {n_points} points spanning {grid_c[0]} to {grid_c[-1]}")
    print(f"Anchors: 5.0, 23.0, 50.0, 100.0, 1000.0")
    print(f"Optimizer: Nelder-Mead (maxiter={maxiter}, xatol={xatol}, fatol={fatol})")

    t_start = time.perf_counter()
    best_vectors: dict[int, np.ndarray] = {}
    best_vals: dict[int, float] = {}
    runtimes: dict[int, float] = {}
    nfevs: dict[int, int] = {}

    # Initialize c = 1000 with known Brief 15 solution
    obj_1000 = make_objective(log10_grid[idx_1000])
    best_vectors[idx_1000] = V_1000_8.copy()
    best_vals[idx_1000] = obj_1000(V_1000_8)
    runtimes[idx_1000] = 0.0
    nfevs[idx_1000] = 1

    # PASS 1: Downward sweep from 1000 to 1
    print("\n--- Pass 1: Downward sweep from c = 1000 to c = 1 ---")
    curr_vec = best_vectors[idx_1000].copy()
    for i in range(idx_1000 - 1, -1, -1):
        c_val = grid_c[i]
        l10_c = log10_grid[i]
        obj = make_objective(l10_c)
        t0 = time.perf_counter()
        res = minimize(obj, curr_vec, method="Nelder-Mead", bounds=BOUNDS_8,
                       options={"maxiter": maxiter, "xatol": xatol, "fatol": fatol})
        dt = time.perf_counter() - t0
        curr_vec = res.x.copy()
        best_vectors[i] = curr_vec.copy()
        best_vals[i] = float(res.fun)
        runtimes[i] = dt
        nfevs[i] = int(res.nfev)
        print(f"[{i:02d}] c={c_val:8.2f} (log10={l10_c:.2f}) -> minimax = {res.fun:6.3f}s | time = {dt:5.2f}s | nfev = {res.nfev}")

    # PASS 2: Upward sweep from 1000 to 10000
    print("\n--- Pass 2: Upward sweep from c = 1000 to c = 10000 ---")
    curr_vec = best_vectors[idx_1000].copy()
    for i in range(idx_1000 + 1, n_points):
        c_val = grid_c[i]
        l10_c = log10_grid[i]
        obj = make_objective(l10_c)
        t0 = time.perf_counter()
        res = minimize(obj, curr_vec, method="Nelder-Mead", bounds=BOUNDS_8,
                       options={"maxiter": maxiter, "xatol": xatol, "fatol": fatol})
        dt = time.perf_counter() - t0
        curr_vec = res.x.copy()
        best_vectors[i] = curr_vec.copy()
        best_vals[i] = float(res.fun)
        runtimes[i] = dt
        nfevs[i] = int(res.nfev)
        print(f"[{i:02d}] c={c_val:8.2f} (log10={l10_c:.2f}) -> minimax = {res.fun:6.3f}s | time = {dt:5.2f}s | nfev = {res.nfev}")

    # PASS 3: Upward sweep from 1 to 1000 (checking hysteresis / alternative basin)
    print("\n--- Pass 3: Upward check sweep from c = 1 to c = 1000 ---")
    curr_vec = V_1_8.copy()
    # At c = 1, compare with known two-open optimum
    obj_0 = make_objective(log10_grid[0])
    val_v1 = obj_0(V_1_8)
    if val_v1 < best_vals[0]:
        print(f"  [00] c=1.0: two-open vector improved minimax from {best_vals[0]:.3f} to {val_v1:.3f}")
        best_vectors[0] = V_1_8.copy()
        best_vals[0] = val_v1
    curr_vec = best_vectors[0].copy()

    for i in range(1, idx_1000 + 1):
        c_val = grid_c[i]
        l10_c = log10_grid[i]
        obj = make_objective(l10_c)
        t0 = time.perf_counter()
        res = minimize(obj, curr_vec, method="Nelder-Mead", bounds=BOUNDS_8,
                       options={"maxiter": maxiter, "xatol": xatol, "fatol": fatol})
        dt = time.perf_counter() - t0
        curr_vec = res.x.copy()
        runtimes[i] += dt
        nfevs[i] += int(res.nfev)
        if res.fun < best_vals[i] - 1e-4:
            print(f"  [{i:02d}] c={c_val:8.2f}: upward pass IMPROVED minimax: {best_vals[i]:6.3f}s -> {res.fun:6.3f}s")
            best_vectors[i] = curr_vec.copy()
            best_vals[i] = float(res.fun)

    # PASS 4: Downward check sweep from 10000 to 1000
    print("\n--- Pass 4: Downward check sweep from c = 10000 to c = 1000 ---")
    curr_vec = best_vectors[n_points - 1].copy()
    for i in range(n_points - 2, idx_1000 - 1, -1):
        c_val = grid_c[i]
        l10_c = log10_grid[i]
        obj = make_objective(l10_c)
        t0 = time.perf_counter()
        res = minimize(obj, curr_vec, method="Nelder-Mead", bounds=BOUNDS_8,
                       options={"maxiter": maxiter, "xatol": xatol, "fatol": fatol})
        dt = time.perf_counter() - t0
        curr_vec = res.x.copy()
        runtimes[i] += dt
        nfevs[i] += int(res.nfev)
        if res.fun < best_vals[i] - 1e-4:
            print(f"  [{i:02d}] c={c_val:8.2f}: downward pass IMPROVED minimax: {best_vals[i]:6.3f}s -> {res.fun:6.3f}s")
            best_vectors[i] = curr_vec.copy()
            best_vals[i] = float(res.fun)

    # Evaluate full candidates across entire grid
    print("\n=== Final Candidate Evaluations Across Profile ===")
    results_list = []
    for i in range(n_points):
        c_val = grid_c[i]
        l10_c = log10_grid[i]
        v9 = np.append(best_vectors[i], l10_c)
        ev = evaluate_candidate(SPEC, v9, l1_max=DEFAULT_L1_MAX, compute_recovery=True, mot_mode="band")
        m = SPEC.factory(v9)
        p_ar = compute_peak_p_ar(m)
        binding_obs = identify_binding_observable(ev)

        row = {
            "index": i,
            "c": c_val,
            "log10_c": l10_c,
            "minimax_sigma": float(ev["max_violation"] * 3.0),
            "binding_observable": binding_obs,
            "hill_slope": float(ev["hill_meas"]),
            "z_hill": float(ev["z_hill"]),
            "ec50_um": float(ev["ec50_um"]),
            "z_ec50": float(ev["z_ec50"]),
            "pp_recovery": float(ev["pp_pred"]),
            "z_pp": float(ev["z_pp"]),
            "po_max": float(ev["po_max"]),
            "z_po": float(ev["z_po"]),
            "mot_ms": float(ev["mot_ms"]),
            "z_mot_band": float(ev["z_mot_band"]),
            "z_mot_midpoint": float(ev["z_mot_midpoint"]),
            "z_taus": [float(x) for x in ev["z_taus"]],
            "z_areas": [float(x) for x in ev["z_areas"]],
            "chi2_areas": float(ev["chi2_areas"]),
            "sqrt_chi2_areas": float(np.sqrt(max(ev["chi2_areas"], 0.0))),
            "p_ar_peak": p_ar,
            "runtime_s": runtimes[i],
            "nfev": nfevs[i],
            "vector": [float(x) for x in v9],
        }
        results_list.append(row)

    total_runtime = time.perf_counter() - t_start
    print(f"\nTotal profiling time: {total_runtime:.1f} s ({total_runtime/60:.2f} min)")

    # Print summary table
    print("\n" + "=" * 95)
    print(f"{'c':>8} | {'log10(c)':>8} | {'Minimax (s)':>11} | {'Binding Obs':>18} | {'nH':>6} | {'z_nH':>6} | {'P(AR) peak':>10} | {'Runtime (s)':>11}")
    print("-" * 95)
    for r in results_list:
        print(f"{r['c']:8.2f} | {r['log10_c']:8.2f} | {r['minimax_sigma']:10.3f}s | {r['binding_observable']:>18} | {r['hill_slope']:6.2f} | {r['z_hill']:+6.2f} | {r['p_ar_peak']*100:9.3f}% | {r['runtime_s']:10.1f}s")
    print("=" * 95)

    # Print full residual rows for key values
    key_targets = [5.0, 23.0, 100.0, 1000.0]
    print("\n=== Full Per-Observable Residual Rows for Key Anchors ===")
    for target in key_targets:
        match = next(r for r in results_list if abs(r["c"] - target) < 1e-2)
        print(f"\n--- c = {match['c']} (log10(c) = {match['log10_c']:.2f}) ---")
        print(f"Minimax: {match['minimax_sigma']:.3f} sigma (satisfies <=3s: {'YES' if match['minimax_sigma'] <= 3.0 else 'NO'})")
        print(f"Binding observable: {match['binding_observable']}")
        print(f"Deactivation tau residuals (z):  tau1 = {match['z_taus'][0]:+.2f}s, tau2 = {match['z_taus'][1]:+.2f}s, tau3 = {match['z_taus'][2]:+.2f}s")
        print(f"Deactivation area residuals (z): a1 = {match['z_areas'][0]:+.2f}s, a2 = {match['z_areas'][1]:+.2f}s, a3 = {match['z_areas'][2]:+.2f}s")
        print(f"Deactivation area quad form:     sqrt(chi2_areas) = {match['sqrt_chi2_areas']:.2f}s")
        print(f"Peak CRC residuals (z):          EC50 = {match['z_ec50']:+.2f}s ({match['ec50_um']:.2f} uM), nH = {match['z_hill']:+.2f}s ({match['hill_slope']:.2f})")
        print(f"Paired-pulse recovery (z):       recovery = {match['z_pp']:+.2f}s ({match['pp_recovery']:.3f})")
        print(f"Gating anchor Po_max (z):        po_max = {match['z_po']:+.2f}s ({match['po_max']:.3f})")
        print(f"Mean open time MOT (z):          band = {match['z_mot_band']:+.2f}s, midpoint = {match['z_mot_midpoint']:+.2f}s ({match['mot_ms']:.2f} ms)")
        print(f"Peak fractional P(AR):           {match['p_ar_peak']:.6f} ({match['p_ar_peak']*100:.3f}%)")

    # Plain answer on c <= 50
    le_50 = [r for r in results_list if r["c"] <= 50.0 + 1e-2]
    satisfies_le_50 = any(r["minimax_sigma"] <= 3.0 for r in le_50)
    print("\n" + "=" * 80)
    print("PLAIN ANSWER ON c <= 50:")
    if satisfies_le_50:
        passing = [r for r in le_50 if r["minimax_sigma"] <= 3.0]
        print(f"YES: Found {len(passing)} point(s) with c <= 50 that satisfy all five observables at <= 3.0s:")
        for p in passing:
            print(f"  c = {p['c']:.2f}: minimax = {p['minimax_sigma']:.3f}s")
    else:
        min_in_le_50 = min(le_50, key=lambda x: x["minimax_sigma"])
        print(f"NO: NO value of c <= 50 satisfies all five observables at <= 3.0s.")
        print(f"Lowest minimax for c <= 50 is {min_in_le_50['minimax_sigma']:.3f}s at c = {min_in_le_50['c']:.2f}.")
    print("=" * 80)

    # Save results to json
    out_path = _ROOT / "scripts" / "coop_profile_results.json"
    with open(out_path, "w") as f:
        json.dump({
            "total_runtime_s": total_runtime,
            "results": results_list,
        }, f, indent=2)
    print(f"\nSaved profile results to {out_path}")

    return {
        "total_runtime_s": total_runtime,
        "results": results_list,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Profile cooperativity factor c")
    parser.add_argument("--maxiter", type=int, default=600, help="Max iterations per Nelder-Mead step")
    args = parser.parse_args()
    run_profile(maxiter=args.maxiter)
