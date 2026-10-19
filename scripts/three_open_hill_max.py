#!/usr/bin/env python3
"""Maximum Achievable Hill Slope Analysis (Agent Brief 18 Task 5 / §3 Question 1).

Measures the maximum achievable Hill coefficient (n_H) across:
1. kinetic_two_open (8 params)
2. kinetic_coop (9 params)
3. kinetic_three_open (series, 11 params)
4. kinetic_three_open_parallel (parallel, 11 params)

Evaluations performed:
(a) Unconstrained over the parameter search box (with EC50 in physiological window 1-100 uM).
(b) Sourced Peak CRC alone (EC50 = 11.6 +/- 0.9 uM, nH = 2.2 +/- 0.4).
(c) At Five-Way Minimax Optima.
(d) Mathematical analysis of independent vs sequential open transitions on Hill slope.
"""
from __future__ import annotations

import pathlib
import sys
import numpy as np
from scipy.optimize import minimize, differential_evolution

_ROOT = pathlib.Path(__file__).resolve().parent.parent
if str(_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(_ROOT / "src"))

from circuitpharm.fitting.comparison import MODEL_SPECS
from circuitpharm.fitting.crc_components import extract_crc_parameters
from circuitpharm.fitting.data import JAHN_EC50_SEM, JAHN_EC50_UM, JAHN_HILL, JAHN_HILL_SEM
from circuitpharm.models.kinetic_three_open import KineticThreeOpenModel, KineticThreeOpenParallelModel


def parallel_factory(v):
    param_names = MODEL_SPECS["kinetic_three_open"].param_names
    p = dict(zip(param_names, v))
    c_val = 10.0 ** p["log10_c"]
    kd_val = 10.0 ** p["log10_kd"]
    eff_val = 10.0 ** p["log10_E"]
    d_val = 10.0 ** p["log10_D"]
    r_val = 10.0 ** p["log10_r"]
    koff_val = 10.0 ** p["log10_koff"]
    alpha_val = 10.0 ** p["log10_alpha"]
    delta_val = 10.0 ** p["log10_delta"]
    gamma_val = 10.0 ** p["log10_gamma"]
    delta2_val = 10.0 ** p["log10_delta2"]
    gamma2_val = 10.0 ** p["log10_gamma2"]

    kon_val = koff_val / kd_val
    koff2_val = koff_val / c_val
    beta_val = eff_val * alpha_val
    d_rate = d_val * r_val
    return KineticThreeOpenParallelModel(
        kon=kon_val,
        koff=koff_val,
        koff2=koff2_val,
        beta=beta_val,
        alpha=alpha_val,
        delta=delta_val,
        gamma=gamma_val,
        beta_p=delta2_val,
        alpha_p=gamma2_val,
        d=d_rate,
        r=r_val,
    )


def evaluate_hill_profiles(scheme_name: str, factory_fn, bounds, seed: int = 42):
    # 1. Unconstrained maximum Hill slope over box subject to EC50 matching Jahn window (5 - 20 uM)
    def neg_hill_jahn_window(v):
        try:
            m = factory_fn(v)
            p = extract_crc_parameters(m)
            if not np.isfinite(p.hill_measured) or p.hill_measured <= 0:
                return 1e6
            if not np.isfinite(p.ec50_um) or p.ec50_um < 5.0 or p.ec50_um > 25.0:
                # Penalty for straying outside Jahn concentration regime
                penalty = 100.0 * (min(abs(p.ec50_um - 11.6), 50.0)) ** 2
                return 1e4 + penalty
            return -float(p.hill_measured)
        except Exception:
            return 1e6

    # 2. Pure CRC-alone objective (z_ec50^2 + z_hill^2)
    def obj_crc(v):
        try:
            m = factory_fn(v)
            p = extract_crc_parameters(m)
            ze = (p.ec50_um - JAHN_EC50_UM) / JAHN_EC50_SEM
            zh = (p.hill_measured - JAHN_HILL) / JAHN_HILL_SEM
            return float(ze ** 2 + zh ** 2)
        except Exception:
            return 1e6

    print(f"[{scheme_name}] Maximizing Hill slope with EC50 in Jahn window (5-25 uM)...", flush=True)
    de_jahn = differential_evolution(
        neg_hill_jahn_window,
        bounds=bounds,
        seed=seed,
        maxiter=20,
        popsize=10,
        tol=1e-3,
        workers=1,
    )
    m_jahn = factory_fn(de_jahn.x)
    p_jahn = extract_crc_parameters(m_jahn)

    print(f"[{scheme_name}] Optimizing Sourced Peak CRC alone...", flush=True)
    de_crc = differential_evolution(
        obj_crc,
        bounds=bounds,
        seed=seed,
        maxiter=20,
        popsize=10,
        tol=1e-3,
        workers=1,
    )
    opt_crc = minimize(obj_crc, de_crc.x, bounds=bounds, method="L-BFGS-B")
    m_crc = factory_fn(opt_crc.x)
    p_crc = extract_crc_parameters(m_crc)

    return {
        "scheme": scheme_name,
        "max_hill_jahn": p_jahn.hill_measured,
        "ec50_jahn": p_jahn.ec50_um,
        "po_max_jahn": m_jahn.po_max,
        "crc_opt_hill": p_crc.hill_measured,
        "crc_opt_ec50": p_crc.ec50_um,
        "crc_opt_z_hill": (p_crc.hill_measured - JAHN_HILL) / JAHN_HILL_SEM,
        "crc_opt_z_ec50": (p_crc.ec50_um - JAHN_EC50_UM) / JAHN_EC50_SEM,
        "crc_opt_po_max": m_crc.po_max,
    }


def main():
    print("=" * 80)
    print("MAXIMUM ACHIEVABLE HILL SLOPE TOURNAMENT (Brief 18 Task 5 / §3 Question 1)")
    print(f"Published target: nH = {JAHN_HILL} +/- {JAHN_HILL_SEM} (Jahn et al. 1997)")
    print("=" * 80 + "\n")

    schemes = [
        ("kinetic_two_open", MODEL_SPECS["kinetic_two_open"].factory, MODEL_SPECS["kinetic_two_open"].bounds),
        ("kinetic_coop", MODEL_SPECS["kinetic_coop"].factory, MODEL_SPECS["kinetic_coop"].bounds),
        ("kinetic_three_open", MODEL_SPECS["kinetic_three_open"].factory, MODEL_SPECS["kinetic_three_open"].bounds),
        ("kinetic_three_open_parallel", parallel_factory, MODEL_SPECS["kinetic_three_open"].bounds),
    ]

    for name, factory, bnds in schemes:
        res = evaluate_hill_profiles(name, factory, bnds)
        print(f"\n--- {name} ---")
        print(f"  * Max nH (with EC50 in Jahn window 5-25 uM): {res['max_hill_jahn']:.3f} (EC50 = {res['ec50_jahn']:.2f} uM, Po,max = {res['po_max_jahn']:.3f})")
        print(f"  * CRC-Alone Optimum: nH = {res['crc_opt_hill']:.3f} (z_nH = {res['crc_opt_z_hill']:+.2f} sigma), EC50 = {res['crc_opt_ec50']:.2f} uM (z_EC50 = {res['crc_opt_z_ec50']:+.2f} sigma)")
        print(f"  * Gap from 2.2 target at CRC-alone: {res['crc_opt_hill'] - 2.2:+.3f}")

    print("\n" + "=" * 80)
    print("HILL SLOPE AT FIVE-WAY MINIMAX OPTIMA (From scripts/joint_fiveway_fit.py)")
    print("=" * 80)
    v_two_mini = np.array([0.523, 0.292, 1.999, -1.882, -0.564, 0.159, -1.142, -1.484])
    m_two = MODEL_SPECS["kinetic_two_open"].factory(v_two_mini)
    p_two = extract_crc_parameters(m_two)
    print(f"kinetic_two_open (minimax opt):")
    print(f"  * nH: {p_two.hill_measured:.3f} (z_nH = {(p_two.hill_measured - JAHN_HILL) / JAHN_HILL_SEM:+.2f} sigma)")
    print(f"  * EC50: {p_two.ec50_um:.2f} uM, Po,max: {m_two.po_max:.3f}")

    v_coop_mini = np.array([2.60289828, -0.03629142, 1.4577392, -1.8505241, 1.93281986, 0.06792833, -0.98895444, -1.42765121, 2.9999705])
    m_coop = MODEL_SPECS["kinetic_coop"].factory(v_coop_mini)
    p_coop = extract_crc_parameters(m_coop)
    print(f"\nkinetic_coop (minimax opt):")
    print(f"  * nH: {p_coop.hill_measured:.3f} (z_nH = {(p_coop.hill_measured - JAHN_HILL) / JAHN_HILL_SEM:+.2f} sigma)")
    print(f"  * EC50: {p_coop.ec50_um:.2f} uM, Po,max: {m_coop.po_max:.3f}")

    v_three_mini = np.array([2.60971152, -0.03575962, 1.45888640, -1.84936750, 1.93211494, 0.07095051, -0.98912109, -1.42705089, 2.99986950, -5.0, 0.000359])
    m_three = MODEL_SPECS["kinetic_three_open"].factory(v_three_mini)
    p_three = extract_crc_parameters(m_three)
    print(f"\nkinetic_three_open (minimax opt):")
    print(f"  * nH: {p_three.hill_measured:.3f} (z_nH = {(p_three.hill_measured - JAHN_HILL) / JAHN_HILL_SEM:+.2f} sigma)")
    print(f"  * EC50: {p_three.ec50_um:.2f} uM, Po,max: {m_three.po_max:.3f}")
    print("=" * 80 + "\n")


if __name__ == "__main__":
    main()
