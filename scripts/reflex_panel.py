"""Multi-seed stretch-reflex drug panel with statistics.

Single runs are noisy (30 Mn neurons, stochastic drive), so every condition is run over
N_SEED independent seeds and reported as mean +/- sd. Reflex gain is the dynamic-phase
Mn response per unit Ia drive; %ctrl is the paired ratio to the control at the same seed.
"""
import sys, os, itertools; sys.path.insert(0,".")
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
    scale, ci, seed = args
    orig = dict(SpinalCircuit.W)
    try:
        SpinalCircuit.W = dict(orig)
        SpinalCircuit.W[("Ia","Mn","ampa")] = orig[("Ia","Mn","ampa")]*scale
        SpinalCircuit.W[("Ia","Mn","nmda")] = orig[("Ia","Mn","nmda")]*scale
        r = R.run(CASES[ci][1], seed=seed)
        gain = (r['mn_dyn']-r['mn_base'])/max(1e-6, r['ia_dyn']-r['ia_base'])
        return (scale, ci, seed, gain, r['mn_dyn'], r['mn_sta'])
    finally:
        SpinalCircuit.W = orig

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
