"""Machine-check embedding and nesting between kinetic_jw95 and kinetic_jw95_free_rates.

Evaluates Model B at (log10_kd, log10_E, log10_D) = (-0.428, 0.104, 1.9721) on Deact + Recovery
under both the 3-parameter spec and the 6-parameter enlarged spec with conventions embedded.
"""
from __future__ import annotations

import numpy as np

from circuitpharm.fitting.comparison import FIT_CONVENTIONS, MODEL_SPECS
from circuitpharm.fitting.data import BARBERIS2007_DEACT, BARBERIS2007_PAIRED_PULSE
from circuitpharm.fitting.diagnostic import diagnostic_log_likelihood


def run_nesting_check():
    v3 = np.array([-0.428, 0.104, 1.9721])
    spec3 = MODEL_SPECS["kinetic_jw95"]
    spec6 = MODEL_SPECS["kinetic_jw95_free_rates"]

    m3 = spec3.factory(v3)

    log10_r = float(np.log10(FIT_CONVENTIONS["r"][0]))
    log10_koff = float(np.log10(FIT_CONVENTIONS["koff"][0]))
    log10_alpha = float(np.log10(FIT_CONVENTIONS["alpha"][0]))

    v6 = np.array([v3[0], v3[1], v3[2], log10_r, log10_koff, log10_alpha])
    m6 = spec6.factory(v6)

    dsets = [BARBERIS2007_DEACT, BARBERIS2007_PAIRED_PULSE]

    ll3 = diagnostic_log_likelihood(m3, dsets, n_params=3, anchor_po_max=True)
    ll6 = diagnostic_log_likelihood(m6, dsets, n_params=6, anchor_po_max=True)

    diff_ll = abs(float(ll3.total_log_likelihood) - float(ll6.total_log_likelihood))

    q3 = m3.q_matrix(100.0)
    q6 = m6.q_matrix(100.0)
    max_q_diff = float(np.max(np.abs(q3 - q6)))

    print("=" * 80)
    print("NESTING AND EMBEDDING CHECK: kinetic_jw95 vs kinetic_jw95_free_rates")
    print("=" * 80)
    print(f"3-parameter point (log10): kd={v3[0]:.4f}, E={v3[1]:.4f}, D={v3[2]:.4f}")
    print(f"Convention values (log10): r={log10_r:.6f}, koff={log10_koff:.6f}, alpha={log10_alpha:.6f}")
    print(f"6-parameter vector:        [{', '.join(f'{x:.6f}' for x in v6)}]")
    print("-" * 80)
    print(f"{'Evaluated as':<40} | {'Total Log-Likelihood':<22}")
    print("-" * 80)
    print(f"{'kinetic_jw95 (3 params)':<40} | {ll3.total_log_likelihood:.6f}")
    print(f"{'embedded in kinetic_jw95_free_rates (6 params)':<40} | {ll6.total_log_likelihood:.6f}")
    print("-" * 80)
    print(f"Difference in log-likelihood: {diff_ll:.2e} nats")
    print(f"Max abs difference in Q matrix: {max_q_diff:.2e}")
    print(f"Nesting holds exactly: {diff_ll < 1e-4 and max_q_diff < 1e-12}")
    print("=" * 80)

    return ll3.total_log_likelihood, ll6.total_log_likelihood, diff_ll, max_q_diff


if __name__ == "__main__":
    run_nesting_check()
