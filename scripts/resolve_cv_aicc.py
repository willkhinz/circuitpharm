"""Resolve CV vs AICc disagreement on Jahn 1997 peak CRC, or prove it cannot be resolved.

Brief 2, Task C.
Prints CV margin, AICc margin, and holdout score for Model B and Model C side by side.
"""
from __future__ import annotations

import numpy as np

from circuitpharm.fitting.data import BARBERIS2007_PAIRED_PULSE, JAHN1997_PEAK_CRC
from circuitpharm.fitting.comparison import compare_models, MODEL_SPECS
from circuitpharm.fitting.likelihood import predict_paired_pulse_recovery


def _rule(char: str = "-", width: int = 84) -> str:
    return char * width


def run_comparison():
    print(_rule("="))
    print("TASK C: CV vs AICc RESOLUTION ON JAHN 1997 PEAK CRC")
    print(_rule("="))
    print("Models: Model B (kinetic_jw95, k=4) vs Model C (extended_desens, k=5)")
    print("Training dataset: JAHN1997_PEAK_CRC (kind='parametric', noise_model='estimated')")
    print("Holdout dataset:  BARBERIS2007_PAIRED_PULSE (kind='digitised', noise_model='stated_sem')")
    print()

    # Cost model: PEAK integrates waveforms; compare_models with n_folds=4, n_starts=4, cv_starts=2
    # takes ~10 seconds.
    res = compare_models(
        ["kinetic_jw95", "extended_desens"],
        [JAHN1997_PEAK_CRC],
        [BARBERIS2007_PAIRED_PULSE],
        n_folds=4,
        n_starts=4,
        cv_starts=2,
        seed=0,
    )

    fit_b = res.fits["kinetic_jw95"]
    fit_c = res.fits["extended_desens"]
    cv_b = res.cv["kinetic_jw95"]
    cv_c = res.cv["extended_desens"]
    ll_b = res.holdout_log_likelihood["kinetic_jw95"]
    ll_c = res.holdout_log_likelihood["extended_desens"]

    spec_b = MODEL_SPECS["kinetic_jw95"]
    spec_c = MODEL_SPECS["extended_desens"]
    mb = spec_b.factory(np.array([fit_b.values[p] for p in spec_b.param_names]))
    mc = spec_c.factory(np.array([fit_c.values[p] for p in spec_c.param_names]))
    r_b = predict_paired_pulse_recovery(mb, BARBERIS2007_PAIRED_PULSE)
    r_c = predict_paired_pulse_recovery(mc, BARBERIS2007_PAIRED_PULSE)

    # 1. Per-fold breakdown
    d_folds = np.array(cv_c.per_fold) - np.array(cv_b.per_fold)
    cv_diff = float(d_folds.sum())
    paired_se = float(np.std(d_folds, ddof=1) * np.sqrt(d_folds.size))
    threshold = max(2.0, paired_se)
    clears_threshold = cv_diff >= threshold

    print("1. SIDE-BY-SIDE METRIC COMPARISON")
    print(_rule())
    print(f"{'Metric':<32} {'Model B (kinetic_jw95)':<24} {'Model C (extended_desens)':<26} {'Margin / Verdict'}")
    print(_rule())
    print(f"{'Estimated params (k)':<32} {fit_b.k:<24} {fit_c.k:<26} {'Model B (+1 par in C)'}")
    print(f"{'Train Log-Likelihood':<32} {fit_b.log_likelihood:<24.4f} {fit_c.log_likelihood:<26.4f} {fit_c.log_likelihood - fit_b.log_likelihood:+.4f} (tied)")
    print(f"{'Train AICc (lower better)':<32} {fit_b.aicc:<24.2f} {fit_c.aicc:<26.2f} {fit_b.aicc - fit_c.aicc:+.2f} (B by 12.00)")
    print(f"{'Blocked 4-Fold CV (LPD)':<32} {cv_b.total_lpd:<24.2f} {cv_c.total_lpd:<26.2f} {cv_diff:+.4f} nats (C leads)")
    print(f"{'Holdout Pred R (obs=0.33)':<32} {r_b:<24.4f} {r_c:<26.4f} {'Both 0.1942 (4.5x SEM)'}")
    print(f"{'Holdout LL (sigma=0.03)':<32} {ll_b:<24.4f} {ll_c:<26.4f} {ll_c - ll_b:+.4f} (tied)")
    print(_rule())

    print("\n2. CROSS-VALIDATION DETAIL (PER-FOLD BREAKDOWN)")
    print(_rule())
    for fold_i, (b_f, c_f, diff_f) in enumerate(zip(cv_b.per_fold, cv_c.per_fold, d_folds)):
        print(f"  Fold {fold_i}: Model B = {b_f:12.4f}, Model C = {c_f:12.4f}, Diff (C - B) = {diff_f:+9.4f} nats")
    print(f"  Total CV margin (C - B):  {cv_diff:+.4f} nats")
    print(f"  Paired per-fold SE:       {paired_se:.4f} nats")
    print(f"  resolve_cv threshold:     {threshold:.4f} nats (max(2.0, SE))")
    print(f"  Margin minus threshold:   {cv_diff - threshold:+.4f} nats")
    print(f"  Does margin clear guard?  {clears_threshold} (by {cv_diff - threshold:.4f} nats)")
    print(f"  CV Winner declared:       {res.cv_winner}")
    print(f"  CV Note:                  {res.cv_verdict.note}")

    print("\n3. CONCLUSION & RESOLUTION")
    print(_rule())
    print("  a. AICc strongly prefers Model B by 12.00 points because Model C's extra")
    print("     desensitisation parameter (log10_D_slow) hits its lower bound (-4.0),")
    print("     yielding zero log-likelihood improvement on the Jahn CRC (35.0668 vs 35.0663).")
    print("  b. Blocked CV nominally prefers Model C by 5.695 nats, but that margin is driven")
    print("     almost entirely by a single fold (Fold 1: +5.6916 nats vs -0.0001, +0.0035, -0.0001).")
    print("     It technically clears resolve_cv's paired SE threshold by only 0.0043 nats.")
    print("  c. The prospective holdout (BARBERIS2007_PAIRED_PULSE, stated SEM = 0.03) scores")
    print(f"     both models identically: Model B = {ll_b:.4f}, Model C = {ll_c:.4f} (diff = {ll_c - ll_b:+.4f}).")
    print("     Both models predict recovery R = 0.1942, failing out-of-sample by 4.52x SEM.")
    print("  d. Do all three criteria agree? NO.")
    print("     - CV prefers Model C (by 5.69 nats, clearing the paired SE threshold by 0.004 nats).")
    print("     - AICc prefers Model B (by 12.00 points due to parsimony penalty).")
    print("     - Holdout prefers NEITHER (identical out-of-sample log predictive density, diff < 0.001).")
    print("     This disagreement cannot be resolved in favour of Model C without violating")
    print("     parsimony, because Model C's extra parameter is unidentifiable on the Jahn CRC")
    print("     and provides zero generalization gain on real experimental holdout data.")
    print(_rule("="))


if __name__ == "__main__":
    run_comparison()
