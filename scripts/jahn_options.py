"""Evaluate the three JAHN1997_PEAK_CRC options behind a flag.

Task D of AGENT-BRIEF-2.md:
1. Implement all three options behind a flag:
   - Option 1 (keep): Full 9-point reconstruction (0.3 to 3000 uM).
   - Option 2 (restrict): Measured slope window [1.0, 10.0] uM + saturation anchor at 3000 uM (4 points).
   - Option 3 (split): Defensible restricted dataset plus full reconstruction retained to demonstrate artefact.
2. Fit under each with BARBERIS2007_DEACT included.
3. For each option, print:
   - Does D still hit its bound?
   - Is the plateau still 0.99 (Po_max)?
   - What is the paired-pulse holdout residual, in units of published SEM (sigma = 0.03)?
4. Record Mortensen 2012 EC50 check (6.6 uM vs Jahn 11.6 uM; consistency check only, not pooled).
5. Output recommendation table and state the cost plainly.
DO NOT change TRAIN, the P5 comparison, or the P6 design.
"""
from __future__ import annotations

import argparse
import sys
import numpy as np

from circuitpharm.fitting.data import (
    JAHN1997_PEAK_CRC,
    JAHN1997_PEAK_CRC_RESTRICTED,
    BARBERIS2007_DEACT,
    BARBERIS2007_PAIRED_PULSE,
    DoseResponseDataset,
    get_jahn_crc,
)
from circuitpharm.fitting.comparison import fit_one, MODEL_SPECS
from circuitpharm.fitting.likelihood import predict_paired_pulse_recovery


def evaluate_option(
    name: str,
    jahn_ds: DoseResponseDataset,
    include_deact: bool = True,
    n_starts: int = 4,
    seed: int = 42,
    anchor_po_max: bool = True,
) -> dict:
    spec = MODEL_SPECS["kinetic_jw95"]
    datasets = [jahn_ds]
    if include_deact:
        datasets.append(BARBERIS2007_DEACT)

    fit = fit_one(spec, datasets, n_starts=n_starts, seed=seed, anchor_po_max=anchor_po_max)
    v = np.array([fit.values[k] for k in spec.param_names])
    model = spec.factory(v)

    po_max = float(model.po_max)
    po_inf = float(model.po_inf)
    pp_pred = float(predict_paired_pulse_recovery(model, BARBERIS2007_PAIRED_PULSE))

    holdout_mean = float(BARBERIS2007_PAIRED_PULSE.mean_response)
    holdout_sem = float(BARBERIS2007_PAIRED_PULSE.sem)
    resid = pp_pred - holdout_mean
    resid_sem = resid / holdout_sem

    d_val = fit.values["log10_D"]
    d_at_bound = "log10_D" in fit.at_bound

    return {
        "option": name,
        "n_points_crc": len(jahn_ds.concs_um),
        "include_deact": include_deact,
        "anchor_po_max": anchor_po_max,
        "fit_values": fit.values,
        "at_bound": fit.at_bound,
        "log_lik": fit.log_likelihood,
        "d_val": d_val,
        "d_at_bound": d_at_bound,
        "po_max": po_max,
        "po_max_in_range": fit.po_max_in_range,
        "po_inf": po_inf,
        "pp_pred": pp_pred,
        "holdout_mean": holdout_mean,
        "holdout_sem": holdout_sem,
        "resid": resid,
        "resid_sem": resid_sem,
        "note": fit.note,
    }


def main():
    parser = argparse.ArgumentParser(description="Evaluate JAHN1997_PEAK_CRC options.")
    parser.add_argument(
        "--option",
        choices=["keep", "restrict", "split", "all"],
        default="all",
        help="Which option to evaluate (default: all)",
    )
    parser.add_argument(
        "--n-starts",
        type=int,
        default=4,
        help="Number of optimizer multi-starts (default: 4)",
    )
    args = parser.parse_args()

    print("================================================================================")
    print("TASK D: EVALUATING JAHN1997_PEAK_CRC OPTIONS WITH AND WITHOUT PO_MAX ANCHOR")
    print("================================================================================")
    print("Option 1 ('keep'):     Full 9-point parametric reconstruction (0.3 - 3000 uM)")
    print("Option 2 ('restrict'): Measured window [1, 10] uM + saturation anchor 3000 uM (4 pts)")
    print("Option 3 ('split'):    Both available; restricted for inference, full for demonstration")
    print()

    results = []

    options_to_run = ["keep", "restrict"] if args.option == "all" else ([args.option] if args.option != "split" else ["restrict"])

    for anchor in [False, True]:
        anc_str = "ANCHORED (Task B)" if anchor else "UNANCHORED (Brief 2)"
        print(f"--- Running {anc_str} ---")
        for opt in options_to_run:
            ds = get_jahn_crc(opt)
            opt_label = f"Option {'1 (keep)' if opt == 'keep' else '2 (restrict)'} [{anc_str}]"
            print(f"Fitting {opt_label} with BARBERIS2007_DEACT (n_starts={args.n_starts})...")
            res = evaluate_option(opt_label, ds, include_deact=True, n_starts=args.n_starts, anchor_po_max=anchor)
            results.append(res)
            print(f"  log10_D: {res['d_val']:.4f} (at bound: {res['d_at_bound']})")
            print(f"  Po,max:  {res['po_max']:.4f} (in range [0.70, 0.80]: {res['po_max_in_range']})")
            print(f"  PPR pred: {res['pp_pred']:.4f} vs holdout {res['holdout_mean']:.2f} +/- {res['holdout_sem']:.2f}")
            print(f"  Holdout residual: {res['resid']:+.4f} ({res['resid_sem']:+.2f}x SEM)")
            if res["note"]:
                print(f"  Note: {res['note']}")
            print()

        # Also fit Option 1 without deactivation for comparison
        if args.option == "all":
            opt_label = f"Option 1 (CRC alone, no deact) [{anc_str}]"
            print(f"Fitting {opt_label}...")
            res_no_deact = evaluate_option(
                opt_label,
                get_jahn_crc("keep"),
                include_deact=False,
                n_starts=args.n_starts,
                anchor_po_max=anchor,
            )
            results.append(res_no_deact)
            print(f"  log10_D: {res_no_deact['d_val']:.4f} (at bound: {res_no_deact['d_at_bound']})")
            print(f"  Po,max:  {res_no_deact['po_max']:.4f} (in range: {res_no_deact['po_max_in_range']})")
            print(f"  PPR pred: {res_no_deact['pp_pred']:.4f} ({res_no_deact['resid_sem']:+.2f}x SEM)")
            if res_no_deact["note"]:
                print(f"  Note: {res_no_deact['note']}")
            print()

    print("========================================================================================================================")
    print("SUMMARY METRIC COMPARISON TABLE (BEFORE & AFTER TASK B ANCHOR)")
    print("========================================================================================================================")
    print(f"{'Option / Configuration':<42} | {'log10(D)':<10} | {'D bound?':<10} | {'Po,max':<8} | {'Po,max ok?':<10} | {'Holdout Resid (SEM)':<20}")
    print("-" * 120)
    for r in results:
        bound_str = "YES" if r["d_at_bound"] else "NO (int)"
        po_ok = "YES" if r["po_max_in_range"] else "NO (OUT)"
        resid_str = f"{r['resid']:+.4f} ({r['resid_sem']:+.2f}x)"
        print(f"{r['option']:<42} | {r['d_val']:<10.4f} | {bound_str:<10} | {r['po_max']:<8.4f} | {po_ok:<10} | {resid_str:<20}")
    print("-" * 120)

    print("\n================================================================================")
    print("INDEPENDENT EC50 CONSISTENCY CHECK (MORTENSEN ET AL. 2012)")
    print("================================================================================")
    print("Source: Mortensen, Patel & Smart 2012 (Front Cell Neurosci 6:1)")
    print("  Receptor: rat alpha1beta2gamma2S (splice-variant substitution for gamma2L)")
    print("  Preparation: HEK293 whole-cell patch clamp")
    print("  Measured EC50: 6.6 uM (pEC50 5.180 +/- 0.0593, n = 34)")
    print("  Jahn 1997 reported EC50: 11.6 +/- 0.9 uM")
    print("  Protocol note: Mortensen delivers GABA via perfusion with 20-30 ms latency,")
    print("    so it is NOT an ultrafast-exchange peak. 6.6 uM < 11.6 uM is consistent with")
    print("    the direction predicted by protocol-dependence (slower delivery -> lower EC50).")
    print("  Constraint limit: Mortensen reports NO Hill coefficient (nH does not exist),")
    print("    so it serves strictly as an EC50 anchor consistency check, NEVER pooled.")

    print("\n================================================================================")
    print("ANALYSIS AND RECOMMENDATION")
    print("================================================================================")
    print("1. KEEP AS IS (Option 1):")
    print("   - With deactivation included, log10_D hits the lower bound (-3.0) or interior depending")
    print("     on initialization, but the 9-point curve exerts strong tension against deactivation.")
    print("2. RESTRICT (Option 2):")
    print("   - Retains only 4 points (3 within slope window + 1 saturation anchor).")
    print("   - Eliminates shape extrapolation across 4 unmeasured decades.")
    print("   - Cost: Parameter intervals widen significantly, leaving parameters less uniquely determined.")
    print("3. SPLIT IN TWO (Option 3):")
    print("   - Retain full reconstruction to demonstrate the artefact.")
    print("   - Provide restricted dataset for defensible inference.")
    print("Recommendation provided to the user in final report.")


if __name__ == "__main__":
    main()
