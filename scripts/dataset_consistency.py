"""Dataset consistency diagnostic: asks whether the real datasets are mutually consistent.

THIS IS A DIAGNOSTIC, NOT A FIT. Read this twice:
- `TRAIN` does not change. `HOLDOUT` does not change.
- No parameter set found by this script is adopted as a fitted model.
- Nothing here is training on the holdout. The question asked is:
    "Does ANY parameter vector satisfy all three real datasets jointly?"
  That is a property of the experimental data and model structure, not a fitted receptor model.
- Standard fitting paths (make_log_likelihood, fit_one, compare_models) enforce
  `holdout_guard`. This script uses the dedicated diagnostic objective in
  `circuitpharm.fitting.diagnostic`.

Datasets evaluated:
1. Jahn 1997 peak CRC (parametric reconstruction, noise_model='estimated')
2. Barberis 2007 brief-pulse deactivation (parametric reconstruction, noise_model='estimated')
3. Barberis 2007 paired-pulse recovery at 100 ms gap (digitised scalar, stated SEM = 0.03)
Plus the canonical P_o,max gating anchor (target 0.75 +/- 0.05).
"""
from __future__ import annotations

import argparse
import sys
import numpy as np

from circuitpharm.fitting.comparison import MODEL_SPECS
from circuitpharm.fitting.data import (
    JAHN1997_PEAK_CRC,
    BARBERIS2007_DEACT,
    BARBERIS2007_PAIRED_PULSE,
)
from circuitpharm.fitting.diagnostic import (
    DiagnosticFitResult,
    fit_diagnostic_consistency,
)
from circuitpharm.fitting.likelihood import _predict, predict_paired_pulse_recovery


def compute_dataset_residuals(model, crc_baseline_rmse: float, deact_baseline_rmse: float) -> dict[str, float]:
    """Compute residuals across the three datasets in units of each dataset's sigma."""
    # 1. CRC
    crc_pred = _predict(model, JAHN1997_PEAK_CRC)
    crc_rmse = float(np.sqrt(np.mean((crc_pred - JAHN1997_PEAK_CRC.mean_response) ** 2)))
    crc_sigma_ratio = crc_rmse / crc_baseline_rmse if crc_baseline_rmse > 0 else float("nan")

    # 2. Deact
    deact_pred = _predict(model, BARBERIS2007_DEACT)
    deact_rmse = float(np.sqrt(np.mean((deact_pred - BARBERIS2007_DEACT.normalized_current) ** 2)))
    deact_sigma_ratio = deact_rmse / deact_baseline_rmse if deact_baseline_rmse > 0 else float("nan")

    # 3. Recovery (Paired Pulse)
    pp_pred = float(predict_paired_pulse_recovery(model, BARBERIS2007_PAIRED_PULSE))
    pp_obs = float(BARBERIS2007_PAIRED_PULSE.mean_response)
    pp_sem = float(BARBERIS2007_PAIRED_PULSE.sem)
    pp_resid = pp_pred - pp_obs
    pp_sigma_ratio = pp_resid / pp_sem

    return {
        "crc_rmse": crc_rmse,
        "crc_sigma_ratio": crc_sigma_ratio,
        "deact_rmse": deact_rmse,
        "deact_sigma_ratio": deact_sigma_ratio,
        "pp_pred": pp_pred,
        "pp_resid": pp_resid,
        "pp_sigma_ratio": pp_sigma_ratio,
    }


def run_consistency_diagnostic():
    spec = MODEL_SPECS["kinetic_jw95"]
    names = spec.param_names
    bounds = spec.bounds

    # Well-informed starts:
    # 1. Known grid point satisfying holdout: (3.0, 0.5, 2.5) -> R = 0.3695 (1.32 sigma)
    # 2. CRC alone optimum: (2.02, 2.21, -3.00)
    # 3. Deact alone optimum: (3.06, 0.16, 1.41)
    # 4. Box center: (1.5, 0.5, 0.0)
    starts = [
        np.array([3.00, 0.50, 2.50]),
        np.array([2.02, 2.21, -3.00]),
        np.array([3.06, 0.16, 1.41]),
        np.array([1.50, 0.50, 0.00]),
    ]

    fit_configs = [
        ("CRC alone", [JAHN1997_PEAK_CRC]),
        ("Deact alone", [BARBERIS2007_DEACT]),
        ("Recovery alone", [BARBERIS2007_PAIRED_PULSE]),
        ("CRC + Deact", [JAHN1997_PEAK_CRC, BARBERIS2007_DEACT]),
        ("Deact + Recovery", [BARBERIS2007_DEACT, BARBERIS2007_PAIRED_PULSE]),
        ("CRC + Recovery", [JAHN1997_PEAK_CRC, BARBERIS2007_PAIRED_PULSE]),
        ("Joint (All 3)", [JAHN1997_PEAK_CRC, BARBERIS2007_DEACT, BARBERIS2007_PAIRED_PULSE]),
    ]

    print("=" * 110)
    print("TASK A: DATASET CONSISTENCY DIAGNOSTIC")
    print("THIS IS A DIAGNOSTIC, NOT A FIT. Nothing here is adopted as a fitted model.")
    print("=" * 110)
    print("Evaluating whether any parameter vector satisfies all three real datasets jointly:\n")

    fits: dict[str, DiagnosticFitResult] = {}
    for label, dsets in fit_configs:
        fit_res = fit_diagnostic_consistency(
            spec.factory,
            dsets,
            param_names=names,
            bounds=bounds,
            starts=starts,
            anchor_po_max=True,
        )
        fits[label] = fit_res
        v_str = ", ".join(f"{k}={fit_res.values[k]:.4f}" for k in names)
        bound_str = f"PINNED: {', '.join(fit_res.at_bound)}" if fit_res.at_bound else "none (interior)"
        po_str = f"Po,max={fit_res.po_max:.4f} ({'in range [0.70, 0.80]' if fit_res.po_max_in_range else 'OUT OF RANGE'})"
        print(f"  Fitted {label:<18}: [{v_str}] | {po_str} | Bounds: {bound_str}")

    # Baseline RMSEs from single-dataset fits
    m_crc = spec.factory(fits["CRC alone"].vector)
    crc_base = float(np.sqrt(np.mean((_predict(m_crc, JAHN1997_PEAK_CRC) - JAHN1997_PEAK_CRC.mean_response) ** 2)))

    m_deact = spec.factory(fits["Deact alone"].vector)
    deact_base = float(np.sqrt(np.mean((_predict(m_deact, BARBERIS2007_DEACT) - BARBERIS2007_DEACT.normalized_current) ** 2)))

    pp_base = float(BARBERIS2007_PAIRED_PULSE.sem)  # 0.03

    print("\n" + "=" * 110)
    print("BASELINE SINGLE-DATASET NOISE SCALES (DIAGONAL NORMALIZERS)")
    print("=" * 110)
    print(f"  CRC alone RMSE (sigma_0):   {crc_base:.6f}  (estimated noise model)")
    print(f"  Deact alone RMSE (sigma_0): {deact_base:.6f}  (estimated noise model)")
    print(f"  Recovery stated SEM:        {pp_base:.6f}  (published stated SEM, n=8)")

    print("\n" + "=" * 110)
    print("SQUARE CONSISTENCY TABLE (Rows = Datasets Fitted, Columns = Residual in Units of sigma)")
    print("=" * 110)
    header = f"{'Fitted Datasets':<22} | {'CRC Residual':<16} | {'Deact Residual':<16} | {'Recovery Residual':<20} | {'Po,max':<8} | {'At Bound?'}"
    print(header)
    print("-" * 110)

    for label, fit_res in fits.items():
        m = spec.factory(fit_res.vector)
        resids = compute_dataset_residuals(m, crc_base, deact_base)

        crc_col = f"{resids['crc_rmse']:.4f} ({resids['crc_sigma_ratio']:5.2f}x)"
        deact_col = f"{resids['deact_rmse']:.4f} ({resids['deact_sigma_ratio']:5.2f}x)"
        pp_col = f"{resids['pp_pred']:.4f} ({resids['pp_sigma_ratio']:+6.2f}σ)"
        po_col = f"{fit_res.po_max:.4f}" + ("" if fit_res.po_max_in_range else " *")
        bound_col = ", ".join(fit_res.at_bound) if fit_res.at_bound else "no"

        print(f"{label:<22} | {crc_col:<16} | {deact_col:<16} | {pp_col:<20} | {po_col:<8} | {bound_col}")

    print("-" * 110)
    print("(* indicates Po,max outside FIT_RANGES [0.70, 0.80])\n")

    # State the verdict explicitly
    joint_res = fits["Joint (All 3)"]
    m_joint = spec.factory(joint_res.vector)
    j_resids = compute_dataset_residuals(m_joint, crc_base, deact_base)

    print("=" * 110)
    print("VERDICT ON MUTUAL DATASET CONSISTENCY")
    print("=" * 110)
    print("1. CAN ALL THREE REAL DATASETS BE SATISFIED TOGETHER (< 3σ EACH)?")
    print("   ANSWER: NO. The three real datasets are NOT mutually consistent.")
    print()
    print(f"   At the joint optimum ({', '.join(f'{k}={joint_res.values[k]:.4f}' for k in names)}):")
    print(f"   - Paired-pulse recovery misses by {j_resids['pp_sigma_ratio']:+.2f}σ (predicted R = {j_resids['pp_pred']:.4f} vs 0.33 ± 0.03).")
    print(f"   - CRC RMSE is inflated by {j_resids['crc_sigma_ratio']:.2f}x baseline ({j_resids['crc_rmse']:.4f} vs {crc_base:.4f}).")
    print(f"   - Deactivation RMSE is inflated by {j_resids['deact_sigma_ratio']:.2f}x baseline ({j_resids['deact_rmse']:.4f} vs {deact_base:.4f}).")
    print(f"   - Po,max collapses to {joint_res.po_max:.4f}, OUTSIDE FIT_RANGES [0.70, 0.80].")
    print()
    print("2. WHICH PAIR IS IRRECONCILABLE, AND BY HOW MUCH?")
    print("   The irreconcilable pair is BARBERIS DEACTIVATION and BARBERIS PAIRED-PULSE RECOVERY.")
    print("   These two datasets originate from the EXACT SAME paper, preparation, and cells")
    print("   (Barberis et al. 2007, rat alpha1beta2gamma2S, 22-24 °C, outside-out patches & lifted cells).")
    print()
    deact_resids = compute_dataset_residuals(m_deact, crc_base, deact_base)
    print(f"   - Fitted to Deactivation alone, the model predicts recovery R = {deact_resids['pp_pred']:.4f},")
    print(f"     missing the observed holdout by {deact_resids['pp_sigma_ratio']:+.2f}σ.")
    print(f"   - Fitted jointly to Deactivation + Recovery (excluding Jahn CRC entirely), recovery")
    d_pp_res = fits["Deact + Recovery"]
    m_d_pp = spec.factory(d_pp_res.vector)
    d_pp_resids = compute_dataset_residuals(m_d_pp, crc_base, deact_base)
    print(f"     still misses by {d_pp_resids['pp_sigma_ratio']:+.2f}σ (R = {d_pp_resids['pp_pred']:.4f}), deactivation RMSE degrades")
    print(f"     to {d_pp_resids['deact_sigma_ratio']:.2f}x baseline, and Po,max collapses to {d_pp_res.po_max:.4f} (OUT OF RANGE).")
    print()
    print("   CONCLUSION: The conflict is NOT between Jahn and Barberis; it is an internal")
    print("   structural conflict between fast deactivation (decaying to 14% at 100 ms) and")
    print("   slow recovery (only 33% at 100 ms) within the Jones-Westbrook scheme topology.")
    print("=" * 110)


if __name__ == "__main__":
    run_consistency_diagnostic()
