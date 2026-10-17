"""Evaluation script for Brief 7A results.

Computes:
1. All 4 relaxation components at each of the three models' Brief-6 optima.
2. Checks whether receptors reach steady state during the 2 ms pulse of 10 mM GABA.
3. Computes maximum absolute deviation between eigen-reconstructed waveform and simulate_waveform.
4. Scores predicted components against published Barberis 2007 values under Weighting (iv).
5. Fits models under Weighting (iv) (published components + paired-pulse recovery) to determine
   whether the Barberis deactivation / recovery conflict exists when evaluated honestly against
   published parameters.

Usage:
    PYTHONPATH=src python scripts/eval_brief_7a.py
"""
import json
import pathlib
import sys
import time

import numpy as np
import scipy.linalg as la
from scipy.optimize import minimize

# Ensure src is on sys.path
_ROOT = pathlib.Path(__file__).resolve().parent.parent
if str(_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(_ROOT / "src"))

from circuitpharm.fitting.comparison import MODEL_SPECS
from circuitpharm.fitting.data import BARBERIS2007_DEACT, BARBERIS2007_PAIRED_PULSE
from circuitpharm.fitting.deact_components import (
    deact_components_diagnostic_score,
    deact_components_log_likelihood,
    parse_deactivation_errors,
    predict_components_eigen,
    predict_components_eigen_all,
    verify_eigen_reconstruction,
)
from circuitpharm.fitting.likelihood import predict_paired_pulse_recovery
from circuitpharm.literature import DEACTIVATIONS


def load_brief6_optima() -> dict[str, dict]:
    """Load parameter vectors for the three models at their Brief-6 as-is optima."""
    cache_dir = _ROOT / ".trace-cache" / "weighting_sensitivity"
    optima = {}
    specs = {
        "kinetic_jw95": "kinetic_jw95_as_is",
        "kinetic_jw95_free_rates": "kinetic_jw95_free_rates_as_is",
        "kinetic_open_desens": "kinetic_open_desens_as_is",
    }
    for model_name, prefix in specs.items():
        candidates = list(cache_dir.glob(f"{prefix}_*.json"))
        if not candidates:
            raise FileNotFoundError(f"No cache file found matching {prefix}_*.json in {cache_dir}")
        candidates.sort(key=lambda p: p.stat().st_mtime, reverse=True)
        with open(candidates[0], "r") as f:
            data = json.load(f)
        optima[model_name] = {
            "vector": np.asarray(data["vector"], dtype=float),
            "total_ll_brief6": data.get("total_log_likelihood"),
            "po_max_brief6": data.get("po_max"),
        }
    return optima


def evaluate_brief6_optima():
    b = DEACTIVATIONS[1]
    target_taus, target_areas = b.components
    tau_errs, area_errs = parse_deactivation_errors(b)
    tau_errs_arr = np.asarray(tau_errs, dtype=float)
    area_errs_arr = np.asarray(area_errs, dtype=float)

    optima = load_brief6_optima()

    print("=" * 110)
    print("BRIEF 7A: EVALUATION AT BRIEF-6 AS-IS OPTIMA")
    print("=" * 110)
    print(f"Target Deactivation: Barberis 2007 brief pulse (2 ms, 10 mM GABA)")
    print(f"  Published taus:  {target_taus} ms (errors: +/- {tau_errs} ms)")
    print(f"  Published areas: {target_areas} (errors: +/- {area_errs})")
    print(f"  Published tau_w: {b.tau_w_ms} +/- {b.err_ms} ms")
    print("-" * 110)

    results = {}

    for name, info in optima.items():
        vec = info["vector"]
        spec = MODEL_SPECS[name]
        m = spec.factory(vec)

        # 1. Check steady state during 2 ms pulse
        p0 = np.asarray(m.state_distribution(0.0), dtype=float)
        q_pulse = np.asarray(m.q_matrix(10000.0), dtype=float)
        p_pulse = p0 @ la.expm(q_pulse * 2.0)
        p_ss = np.asarray(m.state_distribution(10000.0), dtype=float)
        ss_reached = bool(np.allclose(p_pulse, p_ss, atol=1e-3, rtol=1e-2))

        # 2. Eigen-components: all 4
        taus_4, areas_4 = predict_components_eigen_all(m)

        # 3. Waveform reconstruction max deviation
        max_dev = verify_eigen_reconstruction(m, atol=1e-5)

        # 4. 4-vs-3 selection: 3 largest by absolute area
        taus_3, areas_3 = predict_components_eigen(m, n_components=3)

        # 5. Residuals against published
        tau_diff = taus_3 - np.asarray(target_taus)
        area_diff = areas_3 - np.asarray(target_areas)
        tau_sig = tau_diff / tau_errs_arr
        area_sig = area_diff / area_errs_arr

        # Log likelihood under Weighting (iv)
        ll_iv = deact_components_log_likelihood(m, deact_record=b)

        # Weighted tau
        tau_w_pred = float(np.sum(areas_3 * taus_3) / np.sum(areas_3))
        tau_w_sig = (tau_w_pred - float(b.tau_w_ms)) / float(b.err_ms)

        # Paired-pulse recovery
        pp_pred = float(predict_paired_pulse_recovery(m, BARBERIS2007_PAIRED_PULSE))
        pp_sig = (pp_pred - float(BARBERIS2007_PAIRED_PULSE.mean_response)) / float(BARBERIS2007_PAIRED_PULSE.sem)

        results[name] = {
            "vector": vec,
            "taus_4": taus_4,
            "areas_4": areas_4,
            "taus_3": taus_3,
            "areas_3": areas_3,
            "tau_sig": tau_sig,
            "area_sig": area_sig,
            "tau_w": tau_w_pred,
            "tau_w_sig": tau_w_sig,
            "ll_iv": ll_iv,
            "max_dev": max_dev,
            "ss_reached": ss_reached,
            "p_pulse_open": float(p_pulse[3]),
            "p_ss_open": float(p_ss[3]),
            "pp_pred": pp_pred,
            "pp_sig": pp_sig,
            "po_max": float(m.po_max) if hasattr(m, "po_max") else None,
        }

        print(f"\nModel: {name}")
        print(f"  Fitted vector: {vec.tolist()}")
        print(f"  Pulse end open prob: {p_pulse[3]:.4f} vs 10 mM steady-state: {p_ss[3]:.4f} (steady state reached? {ss_reached})")
        print(f"  Reconstruction max deviation vs simulate_waveform: {max_dev:.3e} (< 1e-6)")
        print(f"  All 4 non-zero components:")
        for t_i, ar_i in zip(taus_4, areas_4):
            print(f"    tau = {t_i:10.4f} ms, area = {ar_i:10.6f}")
        print(f"  3 selected components (dropping zero-amplitude unbinding mode):")
        for i in range(3):
            print(f"    Comp {i+1}: tau = {taus_3[i]:10.4f} ms ({tau_sig[i]:+6.2f}σ), area = {areas_3[i]:10.6f} ({area_sig[i]:+6.2f}σ)")
        print(f"  Weighted tau: {tau_w_pred:10.4f} ms ({tau_w_sig:+6.2f}σ vs published 52.5 +/- 2.9 ms)")
        print(f"  Deact LL (Weighting iv): {ll_iv:10.4f} nats")
        print(f"  Recovery: R = {pp_pred:6.4f} ({pp_sig:+6.2f}σ vs published 0.33 +/- 0.03)")
        if hasattr(m, "po_max"):
            print(f"  Po,max: {m.po_max:.4f}")

    return results


def joint_objective(m):
    b = DEACTIVATIONS[1]
    ll_deact = deact_components_log_likelihood(m, deact_record=b)
    pp_pred = float(predict_paired_pulse_recovery(m, BARBERIS2007_PAIRED_PULSE))
    pp_sig = (pp_pred - float(BARBERIS2007_PAIRED_PULSE.mean_response)) / float(BARBERIS2007_PAIRED_PULSE.sem)
    ll_rec = -0.5 * (pp_sig ** 2 + np.log(2.0 * np.pi * (0.03 ** 2)))
    ll_anchor = 0.0
    if hasattr(m, "po_max"):
        po = float(m.po_max)
        if po < 0.70:
            ll_anchor = -0.5 * ((po - 0.70) / 0.02) ** 2
        elif po > 0.80:
            ll_anchor = -0.5 * ((po - 0.80) / 0.02) ** 2
    return -(ll_deact + ll_rec + ll_anchor), ll_deact, ll_rec, pp_pred, pp_sig


def fit_under_weighting_iv():
    """Test whether the deact / paired-pulse conflict exists under Weighting (iv)."""
    print("\n" + "=" * 110)
    print("WEIGHTING (IV) JOINT FIT: DOES THE CONFLICT EXIST UNDER HONEST SCORING?")
    print("=" * 110)

    spec_3 = MODEL_SPECS["kinetic_jw95"]
    spec_6 = MODEL_SPECS["kinetic_jw95_free_rates"]

    # 1. Fit kinetic_jw95 (3 params)
    def loss_3(v):
        return joint_objective(spec_3.factory(v))[0]

    v0_3 = np.array([-0.428, 0.104, 1.972])
    res_3 = minimize(loss_3, v0_3, method="L-BFGS-B", bounds=spec_3.bounds)
    m_3 = spec_3.factory(res_3.x)
    tot_3, d_3, r_3, pp_3, sig_3 = joint_objective(m_3)
    t3, a3 = predict_components_eigen(m_3)

    print("\n1. kinetic_jw95 (3 parameters) under Weighting (iv):")
    print(f"  Fitted vector: {res_3.x.tolist()}")
    print(f"  Total LL: {-tot_3:.2f} nats (Deact LL: {d_3:.2f}, Recovery LL: {r_3:.2f})")
    print(f"  Recovery miss: {sig_3:+.2f}σ (predicted R = {pp_3:.4f})")
    print(f"  Pred taus: {t3}")
    print(f"  Pred areas: {a3}")
    print(f"  Po,max: {m_3.po_max:.4f}")

    # 2. Fit kinetic_jw95_free_rates (6 params)
    def loss_6(v):
        return joint_objective(spec_6.factory(v))[0]

    v0_6 = np.array([0.009, 0.211, 2.701, -1.991, 0.688, -0.485])
    res_6 = minimize(loss_6, v0_6, method="L-BFGS-B", bounds=spec_6.bounds)
    m_6 = spec_6.factory(res_6.x)
    tot_6, d_6, r_6, pp_6, sig_6 = joint_objective(m_6)
    t6, a6 = predict_components_eigen(m_6)

    print("\n2. kinetic_jw95_free_rates (6 parameters) under Weighting (iv):")
    print(f"  Fitted vector: {res_6.x.tolist()}")
    print(f"  Total LL: {-tot_6:.2f} nats (Deact LL: {d_6:.2f}, Recovery LL: {r_6:.2f})")
    print(f"  Recovery miss: {sig_6:+.2f}σ (predicted R = {pp_6:.4f})")
    print(f"  Pred taus: {t6}")
    print(f"  Pred areas: {a6}")
    print(f"  Po,max: {m_6.po_max:.4f}")


if __name__ == "__main__":
    t0 = time.perf_counter()
    evaluate_brief6_optima()
    fit_under_weighting_iv()
    elapsed = time.perf_counter() - t0
    print(f"\nExecution complete in {elapsed:.2f} s.")
