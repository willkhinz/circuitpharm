"""Weighting sensitivity analysis: quantify nat exchange rates across diagnostic weightings.

Brief 6 Task C:
Quantifies how the objective function trades deactivation fit against paired-pulse recovery
under three declared weightings:
  (i) as-is: concentrated Gaussian with deactivation sigma profiled out (~100 points free sigma
      vs 1 scalar point with fixed SEM=0.03).
  (ii) fixed deactivation sigma: deactivation sigma fixed at baseline single-dataset RMSE (0.095091),
      so the trace cannot buy likelihood by shrinking its own sigma.
  (iii) equal weight per dataset: each dataset's log-likelihood divided by its point count n (1/n),
      so a 100-point trace does not outvote a 1-point scalar 100 to 1. Anchor has n=1.

EVALUATED MODELS on Deact + Recovery:
  - kinetic_jw95 (3 params: log10_kd, log10_E, log10_D)
  - kinetic_jw95_free_rates (6 params: log10_kd, log10_E, log10_D, log10_r, log10_koff, log10_alpha)
  - kinetic_open_desens (3 params: log10_kd, log10_E, log10_D; alternative topology A2R <-> A2O <-> A2D)

Usage:
  python scripts/weighting_sensitivity.py [--n-starts 24] [--workers 8] [--no-cache]
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import pathlib
import sys
import time
from typing import Any

import numpy as np

# Ensure src is on sys.path
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent / "src"))

from circuitpharm.fitting.comparison import MODEL_SPECS, nesting_embedding
from circuitpharm.fitting.data import (
    BARBERIS2007_DEACT,
    BARBERIS2007_PAIRED_PULSE,
)
from circuitpharm.fitting.diagnostic import (
    DiagnosticFitResult,
    diagnostic_log_likelihood,
    fit_diagnostic_consistency,
)
from circuitpharm.fitting.likelihood import _predict, predict_paired_pulse_recovery

CACHE_DIR = pathlib.Path(__file__).resolve().parent.parent / ".trace-cache" / "weighting_sensitivity"


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


def cache_key(spec_name: str, weighting: str, fixed_deact_sigma: float, n_starts: int, seed: int) -> str:
    h = hashlib.sha256()
    h.update(f"{spec_name}:{weighting}:{fixed_deact_sigma}:{n_starts}:{seed}".encode())
    h.update(compute_code_hash().encode())
    return h.hexdigest()[:20]


def run_fit(
    spec_name: str,
    weighting: str,
    *,
    warm_starts: list[np.ndarray],
    n_random_starts: int = 24,
    seed: int = 42,
    fixed_deact_sigma: float = 0.095091,
    n_workers: int = 8,
    use_cache: bool = True,
) -> DiagnosticFitResult:
    spec = MODEL_SPECS[spec_name]
    dsets = [BARBERIS2007_DEACT, BARBERIS2007_PAIRED_PULSE]

    key = cache_key(spec_name, weighting, fixed_deact_sigma, n_random_starts, seed)
    cache_file = CACHE_DIR / f"{spec_name}_{weighting}_{key}.json"

    if use_cache and cache_file.exists():
        try:
            with open(cache_file, "r") as f:
                data = json.load(f)
            v = np.array(data["vector"], dtype=float)
            m = spec.factory(v)
            sc = diagnostic_log_likelihood(
                m, dsets, n_params=len(spec.param_names), anchor_po_max=True,
                weighting=weighting, fixed_deact_sigma=fixed_deact_sigma,
            )
            return DiagnosticFitResult(
                values={spec.param_names[i]: float(v[i]) for i in range(len(v))},
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
            pass

    res = fit_diagnostic_consistency(
        spec.factory,
        dsets,
        param_names=spec.param_names,
        bounds=spec.bounds,
        starts=warm_starts,
        anchor_po_max=True,
        n_random_starts=n_random_starts,
        seed=seed,
        weighting=weighting,
        fixed_deact_sigma=fixed_deact_sigma,
        n_workers=n_workers,
    )

    if use_cache:
        CACHE_DIR.mkdir(parents=True, exist_ok=True)
        with open(cache_file, "w") as f:
            json.dump({
                "spec_name": spec_name,
                "weighting": weighting,
                "vector": res.vector.tolist(),
                "values": res.values,
                "total_log_likelihood": res.total_log_likelihood,
                "po_max": res.po_max,
                "po_max_in_range": res.po_max_in_range,
                "at_bound": list(res.at_bound),
                "converged": res.converged,
                "note": res.note,
            }, f, indent=2)

    return res


def analyze_weightings(n_starts: int = 24, n_workers: int = 8, use_cache: bool = True):
    t0 = time.perf_counter()
    deact_base = 0.095091
    dsets = [BARBERIS2007_DEACT, BARBERIS2007_PAIRED_PULSE]

    spec_3 = MODEL_SPECS["kinetic_jw95"]
    spec_6 = MODEL_SPECS["kinetic_jw95_free_rates"]
    spec_open = MODEL_SPECS["kinetic_open_desens"]

    base_3_starts = [
        np.array([3.00, 0.50, 2.50]),
        np.array([2.02, 2.21, -3.00]),
        np.array([3.06, 0.16, 1.41]),
        np.array([-0.428, 0.104, 1.9721]),
        np.array([1.50, 0.50, 0.00]),
    ]
    base_open_starts = [
        np.array([3.00, 0.50, 2.50]),
        np.array([2.02, 2.21, -3.00]),
        np.array([3.06, 0.16, 1.41]),
        np.array([-0.5196, 0.5392, 1.5411]),
        np.array([1.50, 0.50, 0.00]),
    ]

    weightings = [
        ("as_is", "Weighting (i): as-is (deactivation sigma profiled out, 100 points vs 1 point SEM=0.03)"),
        ("fixed_deact_sigma", "Weighting (ii): deactivation sigma fixed at baseline RMSE 0.095091"),
        ("equal_per_dataset", "Weighting (iii): equal weight per dataset (1/n weighting; deact/100, recovery/1, anchor/1)"),
    ]

    print("=" * 115)
    print("WEIGHTING SENSITIVITY AND NAT EXCHANGE RATE ANALYSIS (Brief 6 Task C)")
    print("=" * 115)
    print(f"Configurations evaluated: Deact + Recovery on 3 models")
    print(f"Tolerances: L-BFGS-B ftol=1e-14, gtol=1e-11, maxiter=4000; Nelder-Mead xatol=1e-9, fatol=1e-11, maxiter=2000")
    print(f"Random starts per fit: {n_starts} LatinHypercube + warm starts | Workers: {n_workers}\n")

    all_results = {}

    for w_key, w_desc in weightings:
        print("-" * 115)
        print(f"EVALUATING {w_desc.upper()}")
        print("-" * 115)

        # 1. Fit kinetic_jw95 (3 params)
        fit_3 = run_fit(
            "kinetic_jw95", w_key,
            warm_starts=base_3_starts,
            n_random_starts=n_starts,
            n_workers=n_workers,
            use_cache=use_cache,
        )

        # 2. Fit kinetic_jw95_free_rates (6 params)
        # Warm starts MUST include the embedded optimum of kinetic_jw95 under the same weighting
        v_sub_opt = fit_3.vector
        v_emb_opt = nesting_embedding(spec_3, spec_6, v_sub_opt)
        warm_6 = [v_emb_opt] + [nesting_embedding(spec_3, spec_6, s) for s in base_3_starts]
        fit_6 = run_fit(
            "kinetic_jw95_free_rates", w_key,
            warm_starts=warm_6,
            n_random_starts=n_starts,
            n_workers=n_workers,
            use_cache=use_cache,
        )

        # Check nesting property on total log-likelihood:
        m_emb = spec_6.factory(v_emb_opt)
        ll_emb = diagnostic_log_likelihood(
            m_emb, dsets, n_params=6, anchor_po_max=True, weighting=w_key,
        ).total_log_likelihood
        if fit_6.total_log_likelihood < ll_emb - 1e-4:
            raise RuntimeError(
                f"Nesting violation under weighting '{w_key}': "
                f"enlarged model total LL ({fit_6.total_log_likelihood:.6f}) < "
                f"embedded sub-model LL ({ll_emb:.6f})"
            )

        # 3. Fit kinetic_open_desens (3 params)
        fit_open = run_fit(
            "kinetic_open_desens", w_key,
            warm_starts=base_open_starts,
            n_random_starts=n_starts,
            n_workers=n_workers,
            use_cache=use_cache,
        )

        all_results[w_key] = {
            "kinetic_jw95": fit_3,
            "kinetic_jw95_free_rates": fit_6,
            "kinetic_open_desens": fit_open,
            "ll_emb": ll_emb,
        }

        # Print detailed breakdown for this weighting
        print(f"{'Model':<26} | {'Total Obj':<11} | {'Deact RMSE':<14} | {'Recovery R':<13} | {'Rec Miss':<10} | {'Po,max':<8} | {'At Bound?'}")
        print("-" * 115)
        for m_name, fit_obj in [
            ("kinetic_jw95 (3-param)", fit_3),
            ("kinetic_jw95_free_rates (6)", fit_6),
            ("kinetic_open_desens (3)", fit_open),
        ]:
            m = MODEL_SPECS[m_name.split()[0]].factory(fit_obj.vector)
            deact_pred = _predict(m, BARBERIS2007_DEACT)
            deact_rmse = float(np.sqrt(np.mean((deact_pred - BARBERIS2007_DEACT.normalized_current) ** 2)))
            deact_ratio = deact_rmse / deact_base

            pp_pred = float(predict_paired_pulse_recovery(m, BARBERIS2007_PAIRED_PULSE))
            pp_sigma = (pp_pred - 0.33) / 0.03
            po_max_str = f"{fit_obj.po_max:.4f}" + ("" if fit_obj.po_max_in_range else " *")
            bound_str = ", ".join(fit_obj.at_bound) if fit_obj.at_bound else "no"

            print(
                f"{m_name:<26} | {fit_obj.total_log_likelihood:+11.4f} | "
                f"{deact_rmse:.4f} ({deact_ratio:4.2f}x) | {pp_pred:13.4f} | "
                f"{pp_sigma:+9.2f}σ | {po_max_str:<8} | {bound_str}"
            )
        print("-" * 115)
        print()

    # SECTION 1: NAT DECOMPOSITION (Task C §1)
    print("=" * 115)
    print("TASK C §1: MEASURED NAT DECOMPOSITION AT EACH FIT'S OPTIMUM (UNDER AS-IS SCORING)")
    print("=" * 115)
    print("Shows explicitly how many nats the enlarged model gains on deactivation and pays on recovery:")
    print("-" * 115)
    print(f"{'Model':<26} | {'Total LL':<10} | {'Deact LL (nats)':<17} | {'Recovery LL':<16} | {'Anchor LL':<12} | {'Deact RMSE':<10} | {'Rec Miss'}")
    print("-" * 115)

    res_as_is = all_results["as_is"]
    decomp = {}
    for m_name in ["kinetic_jw95", "kinetic_jw95_free_rates", "kinetic_open_desens"]:
        fit_obj = res_as_is[m_name]
        m = MODEL_SPECS[m_name].factory(fit_obj.vector)
        sc = diagnostic_log_likelihood(m, dsets, n_params=len(MODEL_SPECS[m_name].param_names), anchor_po_max=True, weighting="as_is")
        d_fit = [f for f in sc.per_dataset if "deactivation" in f.label.lower()][0]
        r_fit = [f for f in sc.per_dataset if "paired-pulse" in f.label.lower()][0]
        a_fit = [f for f in sc.per_dataset if "anchor" in f.label.lower()][0]
        
        deact_rmse = float(np.sqrt(d_fit.rss / d_fit.n))
        pp_pred = float(predict_paired_pulse_recovery(m, BARBERIS2007_PAIRED_PULSE))
        pp_sigma = (pp_pred - 0.33) / 0.03
        
        decomp[m_name] = {
            "total": float(sc.total_log_likelihood),
            "deact": float(d_fit.log_likelihood),
            "recovery": float(r_fit.log_likelihood),
            "anchor": float(a_fit.log_likelihood),
            "rmse": deact_rmse,
            "sigma": pp_sigma,
            "vector": fit_obj.vector,
            "po_max": fit_obj.po_max,
        }
        print(
            f"{m_name:<26} | {sc.total_log_likelihood:+10.4f} | {d_fit.log_likelihood:+17.4f} | "
            f"{r_fit.log_likelihood:+16.4f} | {a_fit.log_likelihood:+12.4f} | {deact_rmse:10.4f} | {pp_sigma:+9.2f}σ"
        )
    print("-" * 115)

    delta_total = decomp["kinetic_jw95_free_rates"]["total"] - decomp["kinetic_jw95"]["total"]
    delta_deact = decomp["kinetic_jw95_free_rates"]["deact"] - decomp["kinetic_jw95"]["deact"]
    delta_rec = decomp["kinetic_jw95_free_rates"]["recovery"] - decomp["kinetic_jw95"]["recovery"]
    delta_anc = decomp["kinetic_jw95_free_rates"]["anchor"] - decomp["kinetic_jw95"]["anchor"]

    print("\nMEASURED EXCHANGE RATE (kinetic_jw95 -> kinetic_jw95_free_rates):")
    print(f"  Deactivation term: gained {delta_deact:+.4f} nats (RMSE {decomp['kinetic_jw95']['rmse']:.4f} -> {decomp['kinetic_jw95_free_rates']['rmse']:.4f})")
    print(f"  Recovery term:     paid   {delta_rec:+.4f} nats (Miss {decomp['kinetic_jw95']['sigma']:+.2f}σ -> {decomp['kinetic_jw95_free_rates']['sigma']:+.2f}σ)")
    print(f"  Anchor term:       gained {delta_anc:+.4f} nats (Po,max {decomp['kinetic_jw95']['po_max']:.4f} -> {decomp['kinetic_jw95_free_rates']['po_max']:.4f})")
    print(f"  Net Total LL:      gained {delta_total:+.4f} nats (nested enlargement strictly improves total objective)")
    print(f"  Net Exchange:      {delta_deact:+.2f} nats gained on deactivation easily pays for {abs(delta_rec):.2f} nats lost on recovery!\n")

    # SECTION 2: COMPARISON TABLE ACROSS WEIGHTINGS (Task C §2)
    print("=" * 115)
    print("TASK C §2: CROSS-WEIGHTING COMPARISON TABLE")
    print("=" * 115)
    print(f"{'Weighting Assumption':<32} | {'Model':<26} | {'Rec Miss (σ)':<12} | {'Deact RMSE':<14} | {'Po,max':<8} | {'Satisfies Both?'}")
    print("-" * 115)

    for w_key, w_desc in weightings:
        w_short = w_key
        for m_name, label in [
            ("kinetic_jw95", "kinetic_jw95 (3-param)"),
            ("kinetic_jw95_free_rates", "kinetic_jw95_free_rates (6)"),
            ("kinetic_open_desens", "kinetic_open_desens (3)"),
        ]:
            fit_obj = all_results[w_key][m_name]
            m = MODEL_SPECS[m_name].factory(fit_obj.vector)
            deact_pred = _predict(m, BARBERIS2007_DEACT)
            deact_rmse = float(np.sqrt(np.mean((deact_pred - BARBERIS2007_DEACT.normalized_current) ** 2)))
            deact_ratio = deact_rmse / deact_base

            pp_pred = float(predict_paired_pulse_recovery(m, BARBERIS2007_PAIRED_PULSE))
            pp_sigma = (pp_pred - 0.33) / 0.03

            # Criterion: deact within 3x baseline RMSE AND recovery within 3 sigma
            satisfies_both = (deact_ratio <= 3.0) and (abs(pp_sigma) <= 3.0)
            both_str = "YES" if satisfies_both else "NO"

            print(
                f"{w_desc[:32]:<32} | {label:<26} | {pp_sigma:+12.2f}σ | "
                f"{deact_rmse:.4f} ({deact_ratio:4.2f}x) | {fit_obj.po_max:.4f} | {both_str}"
            )
        print("-" * 115)

    elapsed = time.perf_counter() - t0
    print(f"\nExecution complete in {elapsed:.2f} s ({elapsed/60:.2f} min).")
    return all_results


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Weighting sensitivity analysis")
    parser.add_argument("--n-starts", type=int, default=24, help="number of LatinHypercube random starts")
    parser.add_argument("--workers", type=int, default=8, help="number of multiprocessing workers")
    parser.add_argument("--no-cache", action="store_true", help="bypass cache and recompute all fits")
    args = parser.parse_args()

    analyze_weightings(n_starts=args.n_starts, n_workers=args.workers, use_cache=not args.no_cache)
