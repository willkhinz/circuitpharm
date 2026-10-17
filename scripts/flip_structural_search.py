"""Structural search script for pre-open flip mechanism vs JW95 free rates (Agent Brief 8B).

Compares:
1. `kinetic_jw95_free_rates` (6 free rates)
2. `kinetic_flip` (8 free rates, pre-open conformational step A2R <-> A2F <-> A2O)

Search budget and method:
- Latin Hypercube sampling over spec bounds (default 1500 starts).
- Best 12 starts polished via bounded Nelder-Mead.
- Minimises chi^2 over all six published deactivation numbers from Barberis et al. 2007
  read directly from circuitpharm.literature.DEACTIVATIONS[1] (no literals):
    * tau_1, tau_2, tau_3
    * area_1, area_2, area_3
- Evaluates relaxation using `predict_components_eigen_all` and selecting the 3 modes
  with largest |fractional amplitude|.
- Reports per-quantity sigma residuals, mode counts (|amp| >= 0.10), positivity of amplitudes,
  and parameter bound flags (`at_bound`).
- Evaluates the three evidence cross-checks:
    1. Paired-pulse recovery at 100 ms (BARBERIS2007_PAIRED_PULSE, 0.33 +/- 0.03)
    2. Jahn 1997 peak CRC RMSE (relative to baseline 0.004917)
    3. po_max against physiological band [0.70, 0.80]

Usage:
    PYTHONPATH=src /Users/whinz/Biochemistry/.venv/bin/python scripts/flip_structural_search.py
"""
from __future__ import annotations

import argparse
import pathlib
import sys
import time
from typing import Any

import numpy as np
from scipy.optimize import minimize
from scipy.stats import qmc

# Ensure src is on sys.path
_ROOT = pathlib.Path(__file__).resolve().parent.parent
if str(_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(_ROOT / "src"))

from circuitpharm.fitting.comparison import MODEL_SPECS, ModelSpec
from circuitpharm.fitting.data import BARBERIS2007_PAIRED_PULSE, JAHN1997_PEAK_CRC
from circuitpharm.fitting.deact_components import (
    parse_deactivation_errors,
    predict_components_eigen_all,
)
from circuitpharm.fitting.likelihood import _predict, predict_paired_pulse_recovery
from circuitpharm.gabaa_kinetics import FIT_RANGES
from circuitpharm.literature import DEACTIVATIONS


def check_at_bound(
    vec: np.ndarray,
    bounds: tuple[tuple[float, float], ...],
    param_names: tuple[str, ...],
    tol_rel: float = 1e-3,
) -> list[str]:
    """Identify parameters sitting within tol_rel of their lower or upper bound."""
    pinned = []
    for i, (val, (lo, hi)) in enumerate(zip(vec, bounds)):
        width = hi - lo
        if (val - lo) <= tol_rel * width:
            pinned.append(f"{param_names[i]} (L)")
        elif (hi - val) <= tol_rel * width:
            pinned.append(f"{param_names[i]} (U)")
    return pinned


def evaluate_spec_vector(
    spec: ModelSpec,
    vec: np.ndarray,
    target_taus: np.ndarray,
    target_areas: np.ndarray,
    tau_errs: np.ndarray,
    area_errs: np.ndarray,
) -> dict[str, Any]:
    """Evaluate chi^2 and full diagnostic metrics for a parameter vector."""
    # Box penalty
    penalty = 0.0
    for i, (lo, hi) in enumerate(spec.bounds):
        if vec[i] < lo:
            penalty += 1e6 * (lo - vec[i])
        elif vec[i] > hi:
            penalty += 1e6 * (vec[i] - hi)
    if penalty > 0.0:
        return {"chi2": 1e9 + penalty, "valid": False}

    try:
        model = spec.factory(vec)
        open_idx = getattr(model, "open_index", 3)
        taus_all, areas_all = predict_components_eigen_all(model, open_index=open_idx)

        if len(taus_all) < 3:
            return {"chi2": 1e9, "valid": False}

        # Select top 3 modes by absolute amplitude/area
        top3 = np.argsort(-np.abs(areas_all))[:3]
        order = np.argsort(taus_all[top3])
        chosen = top3[order]

        p_taus = taus_all[chosen]
        p_areas = areas_all[chosen]

        tau_resids_sigma = (p_taus - target_taus) / tau_errs
        area_resids_sigma = (p_areas - target_areas) / area_errs

        chi2_tau = float(np.sum(tau_resids_sigma ** 2))
        chi2_area = float(np.sum(area_resids_sigma ** 2))
        chi2_total = chi2_tau + chi2_area

        # Mode metrics
        appreciable_modes = int(np.sum(np.abs(areas_all) >= 0.10))
        three_positive = bool(np.all(p_areas > 0.0))

        return {
            "chi2": chi2_total,
            "chi2_tau": chi2_tau,
            "chi2_area": chi2_area,
            "pred_taus": p_taus,
            "pred_areas": p_areas,
            "tau_resids_sigma": tau_resids_sigma,
            "area_resids_sigma": area_resids_sigma,
            "all_taus": taus_all,
            "all_areas": areas_all,
            "appreciable_modes": appreciable_modes,
            "total_modes": len(taus_all),
            "three_positive": three_positive,
            "model": model,
            "valid": True,
        }
    except Exception:
        return {"chi2": 1e9, "valid": False}


def run_structural_search(
    spec_name: str,
    n_samples: int = 1500,
    n_polish: int = 12,
    seed: int = 42,
) -> dict[str, Any]:
    """Run Latin-hypercube sampling followed by bounded Nelder-Mead polishing."""
    spec = MODEL_SPECS[spec_name]
    bounds_arr = np.array(spec.bounds)
    d = len(spec.param_names)

    # Read published targets and errors with zero literals
    b = DEACTIVATIONS[1]
    target_taus = np.asarray(b.components[0], dtype=float)
    target_areas = np.asarray(b.components[1], dtype=float)
    raw_tau_errs, raw_area_errs = parse_deactivation_errors(b)
    tau_errs = np.asarray(raw_tau_errs, dtype=float)
    area_errs = np.asarray(raw_area_errs, dtype=float)

    def obj_fn(x: np.ndarray) -> float:
        res = evaluate_spec_vector(spec, x, target_taus, target_areas, tau_errs, area_errs)
        return float(res["chi2"])

    t0 = time.perf_counter()

    # 1. Latin Hypercube Sampling
    sampler = qmc.LatinHypercube(d=d, seed=seed)
    unit_samples = sampler.random(n=n_samples)
    lhs_points = qmc.scale(unit_samples, bounds_arr[:, 0], bounds_arr[:, 1])

    scores = np.array([obj_fn(pt) for pt in lhs_points])
    best_candidate_idx = np.argsort(scores)[:n_polish]

    # 2. Polishing via bounded Nelder-Mead
    best_chi2 = float("inf")
    best_vec = lhs_points[best_candidate_idx[0]].copy()

    for idx in best_candidate_idx:
        start_pt = lhs_points[idx]
        opt_res = minimize(
            obj_fn,
            start_pt,
            method="Nelder-Mead",
            bounds=spec.bounds,
            options={"maxiter": 1000, "xatol": 1e-4, "fatol": 1e-4},
        )
        if opt_res.fun < best_chi2:
            best_chi2 = float(opt_res.fun)
            best_vec = opt_res.x.copy()

    elapsed = time.perf_counter() - t0

    # Final detailed evaluation at best vector
    eval_res = evaluate_spec_vector(spec, best_vec, target_taus, target_areas, tau_errs, area_errs)
    pinned_params = check_at_bound(best_vec, spec.bounds, spec.param_names)

    # Cross-checks
    model = eval_res["model"]
    # 1. Paired-pulse recovery at 100 ms
    pp_pred = float(predict_paired_pulse_recovery(model, BARBERIS2007_PAIRED_PULSE))
    pp_pub = float(BARBERIS2007_PAIRED_PULSE.mean_response)
    pp_sem = float(BARBERIS2007_PAIRED_PULSE.sem)
    pp_sigma = (pp_pred - pp_pub) / pp_sem

    # 2. Jahn 1997 peak CRC RMSE
    crc_pred = _predict(model, JAHN1997_PEAK_CRC)
    crc_obs = np.asarray(JAHN1997_PEAK_CRC.mean_response, dtype=float)
    crc_rmse = float(np.sqrt(np.mean((crc_pred - crc_obs) ** 2)))
    crc_baseline_rmse = 0.004917
    crc_rmse_ratio = crc_rmse / crc_baseline_rmse

    # 3. po_max against FIT_RANGES
    po_min, po_max = FIT_RANGES["po_max"]
    po_max_val = float(model.po_max)
    po_in_band = bool(po_min <= po_max_val <= po_max)

    return {
        "spec_name": spec_name,
        "n_params": d,
        "param_names": spec.param_names,
        "bounds": spec.bounds,
        "best_vec": best_vec,
        "best_chi2": best_chi2,
        "chi2_tau": eval_res["chi2_tau"],
        "chi2_area": eval_res["chi2_area"],
        "pred_taus": eval_res["pred_taus"],
        "pred_areas": eval_res["pred_areas"],
        "target_taus": target_taus,
        "target_areas": target_areas,
        "tau_errs": tau_errs,
        "area_errs": area_errs,
        "tau_resids_sigma": eval_res["tau_resids_sigma"],
        "area_resids_sigma": eval_res["area_resids_sigma"],
        "all_taus": eval_res["all_taus"],
        "all_areas": eval_res["all_areas"],
        "appreciable_modes": eval_res["appreciable_modes"],
        "total_modes": eval_res["total_modes"],
        "three_positive": eval_res["three_positive"],
        "pinned_params": pinned_params,
        "elapsed_s": elapsed,
        "pp_pred": pp_pred,
        "pp_pub": pp_pub,
        "pp_sem": pp_sem,
        "pp_sigma": pp_sigma,
        "crc_rmse": crc_rmse,
        "crc_baseline_rmse": crc_baseline_rmse,
        "crc_rmse_ratio": crc_rmse_ratio,
        "po_max": po_max_val,
        "po_min": po_min,
        "po_max_bound": po_max,
        "po_in_band": po_in_band,
    }


def format_results_markdown(results: list[dict[str, Any]]) -> str:
    """Format comparison results into clean markdown report tables."""
    lines = []
    lines.append("| Scheme | Free Params | Best $\\chi^2$ (6 items) | Modes $\\ge 0.10$ | 3 Positive Amplitudes? | Pinned at Bound (`at_bound`) |")
    lines.append("|---|---|---|---|---|---|")
    for r in results:
        pinned_str = ", ".join(r["pinned_params"]) if r["pinned_params"] else "none (interior)"
        pos_str = "yes" if r["three_positive"] else "no"
        modes_str = f"{r['appreciable_modes']} of {r['total_modes']}"
        lines.append(
            f"| `{r['spec_name']}` | {r['n_params']} | **{r['best_chi2']:.1f}** | "
            f"{modes_str} | {pos_str} | {pinned_str} |"
        )
    lines.append("")

    lines.append("### Component Breakdown and Residuals")
    lines.append("")
    for r in results:
        lines.append(f"#### Model: `{r['spec_name']}` (Best $\\chi^2 = {r['best_chi2']:.2f}$)")
        lines.append("")
        lines.append("| Quantity | Published Target | Published Error | Predicted Value | Residual ($\\sigma$) |")
        lines.append("|---|---|---|---|---|")
        for i in range(3):
            lines.append(
                f"| $\\tau_{i+1}$ (ms) | {r['target_taus'][i]:.2f} | $\\pm {r['tau_errs'][i]:.2f}$ | "
                f"{r['pred_taus'][i]:.4f} | {r['tau_resids_sigma'][i]:+.2f}$\\sigma$ |"
            )
        for i in range(3):
            lines.append(
                f"| $\\text{{area}}_{i+1}$ | {r['target_areas'][i]:.2f} | $\\pm {r['area_errs'][i]:.2f}$ | "
                f"{r['pred_areas'][i]:.6f} | {r['area_resids_sigma'][i]:+.2f}$\\sigma$ |"
            )
        lines.append("")
        lines.append(f"**All relaxation modes:**")
        lines.append(f"* $\\tau$: `[{', '.join(f'{t:.4g}' for t in r['all_taus'])}]` ms")
        lines.append(f"* $\\text{{areas}}$: `[{', '.join(f'{a:.4e}' for a in r['all_areas'])}]`")
        lines.append("")
        lines.append(f"**Fitted parameter vector:**")
        lines.append(f"* `[{', '.join(f'{v:.6f}' for v in r['best_vec'])}]`")
        lines.append(f"* Rates:")
        for name, val in zip(r["param_names"], r["best_vec"]):
            rate_val = 10.0 ** val
            lines.append(f"  * `{name}` = {val:.4f} (value: {rate_val:.4g})")
        lines.append("")

    lines.append("### Task 3: Cross-Checks on Independent Evidence")
    lines.append("")
    lines.append("| Scheme | Paired-Pulse (100 ms) | PP $\\sigma$ Miss | Jahn Peak CRC RMSE | Jahn RMSE Ratio | $P_{o,\\max}$ | $P_{o,\\max} \\in [0.70, 0.80]$? |")
    lines.append("|---|---|---|---|---|---|---|")
    for r in results:
        in_band_str = "YES" if r["po_in_band"] else "NO"
        lines.append(
            f"| `{r['spec_name']}` | {r['pp_pred']:.4f} (pub {r['pp_pub']} $\\pm$ {r['pp_sem']}) | "
            f"{r['pp_sigma']:+.2f}$\\sigma$ | {r['crc_rmse']:.6f} | {r['crc_rmse_ratio']:.1f}$\\times$ | "
            f"{r['po_max']:.4f} | {in_band_str} |"
        )
    lines.append("")

    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(description="Structural search on deactivation components")
    parser.add_argument("--n-samples", type=int, default=1500, help="Number of Latin hypercube starts (default: 1500)")
    parser.add_argument("--n-polish", type=int, default=12, help="Number of top candidates to polish (default: 12)")
    parser.add_argument("--seed", type=int, default=42, help="Random seed for Latin hypercube sampling")
    args = parser.parse_args()

    print("=" * 80)
    print("AGENT BRIEF 8B: STRUCTURAL SEARCH FOR PRE-OPEN FLIP MECHANISM")
    print("=" * 80)
    print(f"Command: {' '.join(sys.argv)}")
    print(f"Search budget: {args.n_samples} Latin Hypercube samples, top {args.n_polish} polished via Nelder-Mead")
    print(f"Random seed: {args.seed}")
    print("-" * 80)

    t_all = time.perf_counter()

    results = []
    for spec_name in ["kinetic_jw95_free_rates", "kinetic_flip"]:
        print(f"\nRunning search for `{spec_name}` ({len(MODEL_SPECS[spec_name].param_names)} parameters)...")
        r = run_structural_search(
            spec_name,
            n_samples=args.n_samples,
            n_polish=args.n_polish,
            seed=args.seed,
        )
        results.append(r)
        print(f"  Best chi^2: {r['best_chi2']:.2f} (elapsed: {r['elapsed_s']:.2f} s)")
        print(f"  Modes >= 0.10: {r['appreciable_modes']} of {r['total_modes']}")
        print(f"  Three positive amplitudes: {r['three_positive']}")
        print(f"  Pinned at bounds: {r['pinned_params'] or 'none'}")
        print(f"  Paired-pulse recovery: {r['pp_pred']:.4f} ({r['pp_sigma']:+.2f} sigma)")
        print(f"  Jahn CRC RMSE: {r['crc_rmse']:.6f} ({r['crc_rmse_ratio']:.1f}x baseline)")
        print(f"  po_max: {r['po_max']:.4f} (in band: {r['po_in_band']})")

    total_time = time.perf_counter() - t_all
    print("\n" + "=" * 80)
    print(f"TOTAL EXECUTION TIME: {total_time:.2f} s")
    print("=" * 80 + "\n")

    md_table = format_results_markdown(results)
    print(md_table)


if __name__ == "__main__":
    main()
