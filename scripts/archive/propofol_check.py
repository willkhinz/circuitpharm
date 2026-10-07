"""Consistency check: with gaba_sens=0.15 fitted on midazolam, does an UNCAPPED
high-efficacy modulator (propofol-like) reach the measured ~-60% tidal volume?

If yes, the calibration is anchored at two points with the correct mechanism distinction:
a ceiling-limited BZ-site PAM gives the midazolam numbers, and a non-ceiling-limited
modulator gives the propofol numbers. Propofol is NOT a BZ-site PAM (it also hits glycine
receptors, HCN and sodium channels and has direct agonist activity at higher
concentrations), so it must NOT be modelled as a capped PAM.
"""
import os
import numpy as np
from multiprocessing import Pool
from circuitpharm.resp import PreBotC, resp_metrics
from circuitpharm.cpg import Drug

OP = dict(drive=170.0, g_adapt=2.5, tau_adapt=400.0, w=dict(ee_ampa=0.45, ee_nmda=0.2475))
SENS, N_SEED = 0.15, 4

def job(a):
    kind, x, seed = a
    cap = 2.5 if kind == "PAM" else 1e9
    d = Drug(gaba_a_gain=x, gaba_a_tau=1.0+0.6*(x-1.0), gaba_a_efficacy_cap=cap)
    b = PreBotC(drug=d, gaba_sens=SENS, seed=seed, **OP)
    for i in range(int(14000/0.1)):
        b.step(0.1)
        if i % 10 == 0: b.record()
    A = b.arrays(); m = A["t"] > 4000
    r = resp_metrics(A["t"][m], A["Out"][m], A["Exc"][m])
    return kind, x, r["mean"], r["amp"], float(r["alive"])

DOSES = (1.0, 2.0, 4.0, 7.0, 11.0, 16.0, 24.0)
if __name__ == "__main__":
    with Pool(os.cpu_count()) as p:
        out = p.map(job, [(k,x,s) for k in ("PAM","uncapped") for x in DOSES for s in range(N_SEED)])
    agg = {}
    for k,x,mn,amp,al in out: agg.setdefault((k,x),[]).append((mn,amp,al))
    cv = np.mean([a[0] for a in agg[("PAM",1.0)]]); ca = np.mean([a[1] for a in agg[("PAM",1.0)]])
    print(f"gaba_sens={SENS}; control vent {cv:.1f}, amp {ca:.1f}\n")
    print(f"{'dose':>5} | {'PAM vent':>8} {'PAM amp':>8} {'st':>4} | "
          f"{'uncap vent':>10} {'uncap amp':>9} {'st':>4}")
    for x in DOSES:
        row = f"{x:5.1f} |"
        for k in ("PAM","uncapped"):
            v = agg[(k,x)]
            mn = np.mean([a[0] for a in v]); am = np.mean([a[1] for a in v])
            al = np.mean([float(a[2]) for a in v])
            st = "APN" if al < 0.5 else "ok"
            row += f" {100*mn/cv-100:7.1f}% {100*am/ca-100:7.1f}% {st:>4} |"
        print(row)
