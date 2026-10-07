"""Subtype-resolved comparison across SUBJECTIVE INTENSITY, normalised to drug-free control.

Two things the previous run obscured:
  1. normalising to the benzodiazepine arm hid the absolute burden;
  2. at a low subjective target every compound is cheap, so the selectivity advantage only
     appears once the required intensity forces the non-selective compound to a high dose.

Also surfaces a real constraint: a highly selective compound with an efficacy ceiling has a
LOWER MAXIMUM achievable subjective effect, because it carries the whole subjective load on
fewer subtypes. max_subj = subjective_index_per_unit_gain * (ceiling - 1).
"""
import sys, os, itertools; sys.path.insert(0, ".")
import numpy as np
from multiprocessing import Pool
from circuitpharm.resp import PreBotC, resp_metrics
from circuitpharm.circuit import SpinalCircuit
from circuitpharm.plant import JointPlant
from circuitpharm.cpg import Drug
from circuitpharm.subtypes import PROFILES
import scripts.reflex as R

RESP_OP = dict(drive=170.0, g_adapt=2.5, tau_adapt=400.0,
               w=dict(ee_ampa=0.45, ee_nmda=0.2475))
IA_SCALE, N_SEED = 0.30, 3
TARGETS = (0.20, 0.35, 0.50, 0.65, 0.80)
ARMS = ["nonselective_bz", "neurosteroid", "hz_166", "l_838417",
        "mp_iii_022", "ideal_a5"]


def gain_for(key, target):
    p = PROFILES[key]; si = p.subjective_index()
    if si <= 1e-9:
        return None
    g = 1.0 + target / si
    return g if g <= p.ceiling else None        # None = unreachable within the ceiling


def max_subj(key):
    p = PROFILES[key]
    return p.subjective_index() * (p.ceiling - 1.0)


def job(a):
    key, target, kind, seed = a
    if key == "CONTROL":
        d, sens_r, sens_s = Drug(), 0.15, 1.0
    else:
        p = PROFILES[key]
        g = gain_for(key, target)
        if g is None:
            return None
        d = Drug(gaba_a_gain=g, gaba_a_tau=1.0 + 0.6 * (g - 1.0),
                 gaba_a_efficacy_cap=p.ceiling)
        sens_r, sens_s = p.regional_sens("prebotc"), p.regional_sens("spinal")
    if kind == "resp":
        b = PreBotC(drug=d, gaba_sens=sens_r, seed=seed, **RESP_OP)
        for i in range(int(14000 / 0.1)):
            b.step(0.1)
            if i % 10 == 0:
                b.record()
        A = b.arrays(); m = A["t"] > 4000
        r = resp_metrics(A["t"][m], A["Out"][m], A["Exc"][m])
        return key, target, kind, r["mean"]
    orig = dict(SpinalCircuit.W)
    try:
        SpinalCircuit.W = dict(orig)
        for rec in ("ampa", "nmda"):
            SpinalCircuit.W[("Ia", "Mn", rec)] = orig[("Ia", "Mn", rec)] * IA_SCALE
        pl = JointPlant(R.XML, drug=d, rg_gain=0.0, seed=seed, gaba_sens=sens_s)
        for (t, q, v) in R.trajectory(0.1, 0.45):
            pl.step(R.DT, impose=(q, v))
        A = pl.arrays(); tt = A["t"] / 1000.0
        base = (tt > 0.15) & (tt < R.PRE)
        dyn = (tt > R.PRE) & (tt < R.PRE + R.RAMP + 0.03)
        g = ((A["mn_E"][dyn].max() - A["mn_E"][base].mean())
             / max(1e-6, A["ia_E"][dyn].max() - A["ia_E"][base].mean()))
        return key, target, kind, g
    finally:
        SpinalCircuit.W = orig


if __name__ == "__main__":
    jobs = [("CONTROL", 0.0, k, s) for k in ("resp", "reflex") for s in range(N_SEED)]
    jobs += [(k, t, kind, s) for k, t in itertools.product(ARMS, TARGETS)
             for kind in ("resp", "reflex") for s in range(N_SEED)]
    print(f"{len(jobs)} sims on {os.cpu_count()} cores")
    with Pool(os.cpu_count()) as p:
        out = [r for r in p.map(job, jobs) if r]
    A = {}
    for key, t, kind, v in out:
        A.setdefault((key, t, kind), []).append(v)
    cv = np.mean(A[("CONTROL", 0.0, "resp")])
    cr = np.mean(A[("CONTROL", 0.0, "reflex")])
    print(f"drug-free control: ventilation {cv:.2f}, reflex gain {cr:.3f}\n")
    print("max achievable subjective effect within each compound's ceiling:")
    for k in ARMS:
        print(f"  {PROFILES[k].name:<48} {max_subj(k):.2f}")
    print()
    hdr = f"{'GABA arm':<48}" + "".join(f"{t:>14.2f}" for t in TARGETS)
    print(hdr); print("-" * len(hdr))
    print("VENTILATION (% of drug-free control)")
    for k in ARMS:
        row = f"{PROFILES[k].name:<48}"
        for t in TARGETS:
            key = (k, t, "resp")
            row += f"{100*np.mean(A[key])/cv:13.0f}%" if key in A else f"{'--':>14}"
        print(row)
    print("\nSTRETCH REFLEX (% of drug-free control)")
    for k in ARMS:
        row = f"{PROFILES[k].name:<48}"
        for t in TARGETS:
            key = (k, t, "reflex")
            row += f"{100*np.mean(A[key])/cr:13.0f}%" if key in A else f"{'--':>14}"
        print(row)
    print("\n'--' = subjective target unreachable within that compound's efficacy ceiling")
