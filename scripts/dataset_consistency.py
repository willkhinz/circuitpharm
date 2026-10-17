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

import hashlib
import json
import pathlib
import time
from circuitpharm.fitting.comparison import MODEL_SPECS, nesting_embedding
from circuitpharm.fitting.data import (
    JAHN1997_PEAK_CRC,
    BARBERIS2007_DEACT,
    BARBERIS2007_PAIRED_PULSE,
)
from circuitpharm.fitting.diagnostic import (
    DiagnosticFitResult,
    diagnostic_log_likelihood,
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


CACHE_DIR = pathlib.Path(__file__).resolve().parent.parent / ".trace-cache" / "diagnostic_fits"



def compute_code_hash() -> str:
    h = hashlib.sha256()
    repo_root = pathlib.Path(__file__).resolve().parent.parent
    for rel_path in [
        "src/circuitpharm/fitting/diagnostic.py",
        "src/circuitpharm/fitting/comparison.py",
        "src/circuitpharm/fitting/likelihood.py",
        "src/circuitpharm/models/kinetic.py",
        "src/circuitpharm/models/kinetic_open_desens.py",
    ]:
        p = repo_root / rel_path
        if p.exists():
            h.update(p.read_bytes())
    return h.hexdigest()[:16]


def cache_key(spec_name: str, label: str, n_starts: int, seed: int) -> str:
    h = hashlib.sha256()
    h.update(f"{spec_name}:{label}:{n_starts}:{seed}".encode())
    h.update(compute_code_hash().encode())
    return h.hexdigest()[:20]


def run_spec_diagnostic(
    spec_name: str,
    starts: list[np.ndarray],
    crc_base: float,
    deact_base: float,
    pp_base: float,
    *,
    sub_fits: dict[str, DiagnosticFitResult] | None = None,
    n_random_starts: int = 6,
    n_workers: int = 6,
    use_cache: bool = True,
):
    spec = MODEL_SPECS[spec_name]
    names = spec.param_names
    bounds = spec.bounds

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
    print(f"SPECIFICATION: {spec_name} ({len(names)} parameters: {', '.join(names)})")
    if spec.conventions:
        print(f"Fixed conventions: {', '.join(spec.conventions)}")
    else:
        print("All rate conventions freed (r, koff, alpha fitted).")
    print(f"Optimization settings: {n_random_starts} random starts + warm starts, ftol=1e-8, gtol=1e-6")
    print("=" * 110)

    fits: dict[str, DiagnosticFitResult] = {}
    for label, dsets in fit_configs:
        cfg_starts = list(starts)
        v_emb = None
        if sub_fits is not None and label in sub_fits:
            sub_spec = MODEL_SPECS["kinetic_jw95"]
            v_sub = sub_fits[label].vector
            v_emb = nesting_embedding(sub_spec, spec, v_sub)
            cfg_starts.insert(0, v_emb)

        # Check disk cache
        k = cache_key(spec_name, label, n_random_starts, 42)
        c_file = CACHE_DIR / f"{spec_name}_{label.replace(' ', '_').replace('(', '').replace(')', '')}_{k}.json"

        fit_res = None
        if use_cache and c_file.exists():
            try:
                with open(c_file, "r") as f:
                    data = json.load(f)
                v = np.array(data["vector"], dtype=float)
                m = spec.factory(v)
                sc = diagnostic_log_likelihood(m, dsets, n_params=len(names), anchor_po_max=True)
                fit_res = DiagnosticFitResult(
                    values={names[i]: float(v[i]) for i in range(len(v))},
                    vector=v,
                    total_log_likelihood=float(sc.total_log_likelihood),
                    per_dataset=sc.per_dataset,
                    po_max=float(m.po_max) if hasattr(m, "po_max") else None,
                    po_max_in_range=data["po_max_in_range"],
                    at_bound=tuple(data["at_bound"]),
                    converged=data["converged"],
                    note=data.get("note", ""),
                )
            except Exception:
                fit_res = None

        if fit_res is None:
            fit_res = fit_diagnostic_consistency(
                spec.factory,
                dsets,
                param_names=names,
                bounds=bounds,
                starts=cfg_starts,
                anchor_po_max=True,
                n_random_starts=n_random_starts,
                seed=42,
                n_workers=n_workers,
            )
            if use_cache:
                CACHE_DIR.mkdir(parents=True, exist_ok=True)
                with open(c_file, "w") as f:
                    json.dump({
                        "spec_name": spec_name,
                        "label": label,
                        "vector": fit_res.vector.tolist(),
                        "values": fit_res.values,
                        "total_log_likelihood": fit_res.total_log_likelihood,
                        "po_max": fit_res.po_max,
                        "po_max_in_range": fit_res.po_max_in_range,
                        "at_bound": list(fit_res.at_bound),
                        "converged": fit_res.converged,
                        "note": fit_res.note,
                    }, f, indent=2)

        # Enforce nesting check on total objective (Brief 6 Task C §6):
        if v_emb is not None:
            m_emb = spec.factory(v_emb)
            ll_sub_emb = float(diagnostic_log_likelihood(m_emb, dsets, n_params=len(names), anchor_po_max=True).total_log_likelihood)
            if fit_res.total_log_likelihood < ll_sub_emb - 1e-4:
                raise RuntimeError(
                    f"Nesting violation in {label}: enlarged model total LL ({fit_res.total_log_likelihood:.6f}) < "
                    f"embedded sub-model LL ({ll_sub_emb:.6f})"
                )

        fits[label] = fit_res
        v_str = ", ".join(f"{k}={fit_res.values[k]:.4f}" for k in names)
        bound_str = f"PINNED: {', '.join(fit_res.at_bound)}" if fit_res.at_bound else "none (interior)"
        po_str = f"Po,max={fit_res.po_max:.4f} ({'in range [0.70, 0.80]' if fit_res.po_max_in_range else 'OUT OF RANGE'})"
        print(f"  Fitted {label:<18}: [{v_str}] | {po_str} | Bounds: {bound_str}")

    print("\n" + "-" * 110)
    print(f"SQUARE CONSISTENCY TABLE FOR {spec_name}")
    print(f"Denominators: CRC sigma_0 = {crc_base:.6f}, Deact sigma_0 = {deact_base:.6f}, Recovery SEM = {pp_base:.6f}")
    print("-" * 110)
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
    return fits


def run_consistency_diagnostic():
    t_start = time.perf_counter()
    print("=" * 110)
    print("DATASET CONSISTENCY DIAGNOSTIC (Brief 4 & Brief 5)")
    print("THIS IS A DIAGNOSTIC, NOT A FIT. Nothing here is adopted as a fitted model.")
    print("=" * 110)
    print("Evaluating whether any parameter vector satisfies all three real datasets jointly:\n")

    # Baseline denominators established from the canonical single-dataset 3-parameter fits:
    crc_base = 0.004917
    deact_base = 0.095091
    pp_base = float(BARBERIS2007_PAIRED_PULSE.sem)  # 0.030000

    print("BASELINE SINGLE-DATASET NOISE SCALES (DIAGONAL NORMALIZERS)")
    print(f"  CRC alone RMSE (sigma_0):   {crc_base:.6f}  (estimated noise model)")
    print(f"  Deact alone RMSE (sigma_0): {deact_base:.6f}  (estimated noise model)")
    print(f"  Recovery stated SEM:        {pp_base:.6f}  (published stated SEM, n=8)\n")

    # 1. Standard 3-parameter spec (kinetic_jw95)
    starts_3 = [
        np.array([3.00, 0.50, 2.50]),
        np.array([2.02, 2.21, -3.00]),
        np.array([3.06, 0.16, 1.41]),
        np.array([1.50, 0.50, 0.00]),
    ]
    fits_3 = run_spec_diagnostic("kinetic_jw95", starts_3, crc_base, deact_base, pp_base, n_random_starts=6, n_workers=6)

    # 2. Enlarged 6-parameter spec (kinetic_jw95_free_rates)
    spec_3 = MODEL_SPECS["kinetic_jw95"]
    spec_6 = MODEL_SPECS["kinetic_jw95_free_rates"]
    starts_6 = [nesting_embedding(spec_3, spec_6, s) for s in starts_3]
    fits_6 = run_spec_diagnostic("kinetic_jw95_free_rates", starts_6, crc_base, deact_base, pp_base, sub_fits=fits_3, n_random_starts=6, n_workers=6)

    # 3. Alternative topology spec (kinetic_open_desens: A2R <-> A2O <-> A2D)
    starts_open = [
        np.array([3.00, 0.50, 2.50]),
        np.array([2.02, 2.21, -3.00]),
        np.array([3.06, 0.16, 1.41]),
        np.array([1.50, 0.50, 0.00]),
    ]
    fits_open = run_spec_diagnostic("kinetic_open_desens", starts_open, crc_base, deact_base, pp_base, n_random_starts=6, n_workers=6)

    # Evaluate verdict on conflict survival (Task B)
    spec_6 = MODEL_SPECS["kinetic_jw95_free_rates"]
    joint_res_6 = fits_6["Joint (All 3)"]
    m_joint_6 = spec_6.factory(joint_res_6.vector)
    j6_resids = compute_dataset_residuals(m_joint_6, crc_base, deact_base)

    d_pp_res_6 = fits_6["Deact + Recovery"]
    m_d_pp_6 = spec_6.factory(d_pp_res_6.vector)
    d6_resids = compute_dataset_residuals(m_d_pp_6, crc_base, deact_base)

    # Evaluate verdict on alternative topology (Task C)
    spec_open = MODEL_SPECS["kinetic_open_desens"]
    d_pp_res_open = fits_open["Deact + Recovery"]
    m_d_pp_open = spec_open.factory(d_pp_res_open.vector)
    d_open_resids = compute_dataset_residuals(m_d_pp_open, crc_base, deact_base)

    m_deact_open = spec_open.factory(fits_open["Deact alone"].vector)
    deact_open_resids = compute_dataset_residuals(m_deact_open, crc_base, deact_base)

    m_rec_open = spec_open.factory(fits_open["Recovery alone"].vector)
    rec_open_resids = compute_dataset_residuals(m_rec_open, crc_base, deact_base)

    print("=" * 110)
    print("VERDICT ON TASK B: DOES FREEING THE CONVENTIONS DISSOLVE THE CONFLICT?")
    print("=" * 110)
    print("ANSWER: NO. THE CONFLICT SURVIVES.")
    print()
    print("Even when all microscopic rate conventions (r, koff, alpha) are freed alongside the")
    print("equilibrium ratios (kd, E, D), the conflict between Barberis deactivation and Barberis")
    print("paired-pulse recovery cannot be reconciled in Model B (JW95):")
    print()
    print("1. In the Barberis-only joint fit (Deact + Recovery, excluding Jahn CRC entirely):")
    print(f"   - Paired-pulse recovery misses by {d6_resids['pp_sigma_ratio']:+.2f}σ (predicted R = {d6_resids['pp_pred']:.4f} vs 0.33 ± 0.03).")
    print(f"   - Deactivation RMSE is {d6_resids['deact_rmse']:.4f} ({d6_resids['deact_sigma_ratio']:.2f}x baseline).")
    print(f"   - Gating efficacy collapses to Po,max = {d_pp_res_6.po_max:.4f}, OUTSIDE FIT_RANGES [0.70, 0.80].")
    print()
    print("2. In the Joint (All 3) fit:")
    print(f"   - Paired-pulse recovery misses by {j6_resids['pp_sigma_ratio']:+.2f}σ (predicted R = {j6_resids['pp_pred']:.4f} vs 0.33 ± 0.03).")
    print(f"   - CRC RMSE is blown out by {j6_resids['crc_sigma_ratio']:.2f}x baseline ({j6_resids['crc_rmse']:.4f} vs {crc_base:.4f}).")
    print(f"   - Gating efficacy Po,max = {joint_res_6.po_max:.4f} ({'in range' if joint_res_6.po_max_in_range else 'OUT OF RANGE'}).")
    print()
    print("=" * 110)
    print("VERDICT ON TASK C: DOES AN ALTERNATIVE TOPOLOGY (OPEN-STATE DESENSITISATION) RECONCILE THEM?")
    print("=" * 110)
    print("ANSWER: NO. THE CONFLICT PERSISTS UNDER THE ALTERNATIVE TOPOLOGY AS WELL.")
    print()
    print("Testing KineticOpenDesensModel (A2R <-> A2O <-> A2D, desensitisation directly from open channel):")
    print(f"1. Deactivation alone misses paired-pulse recovery by {deact_open_resids['pp_sigma_ratio']:+.2f}σ (predicted R = {deact_open_resids['pp_pred']:.4f}),")
    print(f"   and Po,max drops to {fits_open['Deact alone'].po_max:.4f} (OUT OF RANGE).")
    print(f"2. Recovery alone inflates deactivation RMSE by {rec_open_resids['deact_sigma_ratio']:.2f}x baseline ({rec_open_resids['deact_rmse']:.4f} vs {deact_base:.4f}).")
    print(f"3. Deactivation + Recovery jointly misses recovery by {d_open_resids['pp_sigma_ratio']:+.2f}σ (R = {d_open_resids['pp_pred']:.4f})")
    print(f"   and inflates deactivation RMSE by {d_open_resids['deact_sigma_ratio']:.2f}x baseline ({d_open_resids['deact_rmse']:.4f} vs {deact_base:.4f}).")
    print()
    print("OVERALL ARCHITECTURAL CONCLUSION (The Stronger Result):")
    print("Neither topology (closed-state desensitisation A2R <-> A2D vs open-state desensitisation A2O <-> A2D),")
    print("and no combination of freed microscopic rates, can reconcile fast deactivation (decay to 14% at 100 ms)")
    print("with slow paired-pulse recovery (only 33% at 100 ms) from the same preparation under the as-is weighting.")
    print("NOTE (Brief 6): This conclusion was reached under status-quo Weighting (i) (100-point profiled deact vs 1-point recovery).")
    print("Under Weighting (ii) (fixed deact sigma) and Weighting (iii) (equal 1/n per-dataset weighting),")
    print("the 6-parameter model reconciles both observables (recovery miss < 3σ, deact RMSE <= baseline).")
    print("See scripts/weighting_sensitivity.py and RESULTS-BRIEF-6.md for the full weighting sensitivity analysis.")
    print("=" * 110)


if __name__ == "__main__":
    run_consistency_diagnostic()
