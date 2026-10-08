"""Multi-seed stretch-reflex drug panel with statistics.

Single runs are noisy (30 Mn neurons, stochastic drive), so every condition is run over
N_SEED independent seeds and reported as mean +/- sd. Reflex gain is the dynamic-phase
Mn response per unit Ia drive; %ctrl is the paired ratio to the control at the same seed.
"""
import os, itertools
import numpy as np
from multiprocessing import Pool
from circuitpharm.circuit import SpinalCircuit
from circuitpharm.cpg import Drug
import scripts.reflex as R

N_SEED = 6
SCALES = (0.22, 0.30, 0.40)

CASES = [
    ("control",                        Drug()),
    ("GABA-A PAM 1.5x",                Drug(gaba_a_gain=1.5, gaba_a_tau=1.3)),
    ("GABA-A PAM 2.0x (benzo-like)",    Drug(gaba_a_gain=2.0, gaba_a_tau=1.6)),
    ("GlyR PAM 1.6x (ethanol-like)",    Drug(glyr_gain=1.6)),
    ("GlyR block 0.4x (strychnine)",    Drug(glyr_gain=0.4)),
    ("NMDA 60% non-selective",          Drug(nmda_block=0.6, glun2b_selectivity=0.0)),
    ("NMDA 60% GluN2B-sel (spinal)",    Drug(nmda_block=0.6, glun2b_selectivity=1.0,
                                             glun2b_fraction=0.15)),
    ("PAM 2x + GluN2B-sel 60%",         Drug(gaba_a_gain=2.0, gaba_a_tau=1.6, nmda_block=0.6,
                                             glun2b_selectivity=1.0, glun2b_fraction=0.15)),
    ("PAM 2x + non-selective 60%",      Drug(gaba_a_gain=2.0, gaba_a_tau=1.6, nmda_block=0.6,
                                             glun2b_selectivity=0.0)),
]

def job(args):
    """Sweep the Ia->Mn weight scale.

    TWO BUGS FIXED HERE (2026-10-07):

    1. DOUBLE SCALING. This mutated SpinalCircuit.W by `scale` and then called R.run,
       which is `stretch_reflex` -- and that applies its OWN `ia_scale`, defaulting to
       config.IA_SCALE = 0.30. The weight actually reaching the plant was therefore
       scale * 0.30, i.e. over 3x weaker than the swept value, so every point on this
       sweep was mislabelled. Now the scale is passed through as `ia_scale`, the single
       place it is applied.

    2. GLOBAL CLASS-STATE MUTATION. Rewriting SpinalCircuit.W corrupts weights for every
       other instance in the process and leaks permanently if an exception lands before
       the `finally`. SpinalCircuit now takes per-instance weights, so nothing is mutated.

    The reflex gain is taken from `stretch_reflex` rather than recomputed, so it inherits
    the guard that returns NaN when the Ia response is non-positive instead of dividing by
    1e-6 and reporting a ~5,000,000 gain.
    """
    scale, ci, seed = args
    r = R.run(CASES[ci][1], seed=seed, ia_scale=scale)
    return (scale, ci, seed, r["gain"], r["mn_dyn"], r["mn_sta"])

if __name__ == "__main__":
    jobs = list(itertools.product(SCALES, range(len(CASES)), range(N_SEED)))
    with Pool(os.cpu_count()) as p:
        out = p.map(job, jobs)
    D = {}
    for sc, ci, sd, g, dyn, sta in out:
        D.setdefault((sc,ci), {})[sd] = (g, dyn, sta)
    for sc in SCALES:
        ctrl = np.array([D[(sc,0)][s][0] for s in range(N_SEED)])
        print(f"\n=== Ia->Mn scale {sc}  (control Mn dyn "
              f"{np.mean([D[(sc,0)][s][1] for s in range(N_SEED)]):.1f} Hz, ceiling 125) ===")
        print(f"{'condition':<32} {'gain':>14} {'%ctrl':>14} {'Mn dyn':>12} {'Mn sta':>12}")
        for ci,(lab,_) in enumerate(CASES):
            g   = np.array([D[(sc,ci)][s][0] for s in range(N_SEED)])
            dyn = np.array([D[(sc,ci)][s][1] for s in range(N_SEED)])
            sta = np.array([D[(sc,ci)][s][2] for s in range(N_SEED)])
            pct = 100*g/np.maximum(1e-9, ctrl)          # paired per seed
            print(f"{lab:<32} {g.mean():6.3f}±{g.std():5.3f} "
                  f"{pct.mean():7.0f}±{pct.std():4.0f}%  {dyn.mean():6.1f}±{dyn.std():4.1f} "
                  f"{sta.mean():6.2f}±{sta.std():4.2f}")
