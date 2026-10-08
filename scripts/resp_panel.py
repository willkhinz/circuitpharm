"""Respiratory dose-response and apnoea threshold -- the safety axis of the margin.

Brainstem is GluN2D-dominant, so glun2b_fraction=0.15 here. The forebrain (subjective)
comparison uses 0.7. Multi-seed with statistics; apnoea is declared on collapse of burst
modulation or mean output, not on crossing count (see resp_metrics docstring).
"""
import os
import numpy as np
from multiprocessing import Pool
from circuitpharm.resp import PreBotC, resp_metrics
from circuitpharm.cpg import Drug
from circuitpharm.config import RESP_OP

# Operating point from circuitpharm.config -- the single source. Duplicated
# literals are HOW recurring error E12 happened twice.
OP = dict(RESP_OP); OP["w"] = dict(RESP_OP["w"])
N_SEED, T, DT = 4, 14000.0, 0.1
BRAINSTEM_2B = 0.15


def one(args):
    label, drug, seed = args
    b = PreBotC(drug=drug, seed=seed, **OP)
    for i in range(int(T / DT)):
        b.step(DT)
        if i % 10 == 0:
            b.record()
    A = b.arrays(); m = A["t"] > 4000
    r = resp_metrics(A["t"][m], A["Out"][m], A["Exc"][m])
    return (label, r["freq"], r["amp"], r["mean"], r["mod"], r["alive"])


def pam(x):  return Drug(gaba_a_gain=x, gaba_a_tau=1.0 + 0.6 * (x - 1.0))
def nsel(x): return Drug(nmda_block=x, glun2b_selectivity=0.0)
def sel(x):  return Drug(nmda_block=x, glun2b_selectivity=1.0,
                         glun2b_fraction=BRAINSTEM_2B)
def combo(g, b, s):
    return Drug(gaba_a_gain=g, gaba_a_tau=1.0 + 0.6 * (g - 1.0), nmda_block=b,
                glun2b_selectivity=s,
                glun2b_fraction=BRAINSTEM_2B if s else 0.7)

CASES = ([("control", Drug())]
         + [(f"PAM {x:.1f}x", pam(x)) for x in (1.5, 2.0, 2.5, 3.0, 4.0)]
         + [(f"NMDA non-sel {x:.0%}", nsel(x)) for x in (0.3, 0.5, 0.7, 0.9)]
         + [(f"NMDA 2B-sel {x:.0%}", sel(x)) for x in (0.3, 0.5, 0.7, 0.9)]
         + [("PAM 2.0x + 2B-sel 60%",  combo(2.0, 0.6, 1.0)),
            ("PAM 2.0x + non-sel 60%", combo(2.0, 0.6, 0.0)),
            ("PAM 2.5x + 2B-sel 90%",  combo(2.5, 0.9, 1.0)),
            ("PAM 2.5x + non-sel 70%", combo(2.5, 0.7, 0.0))])

if __name__ == "__main__":
    jobs = [(l, d, s) for (l, d) in CASES for s in range(N_SEED)]
    with Pool(os.cpu_count()) as p:
        out = p.map(one, jobs)
    agg = {}
    for lab, f, a, mn, md, alive in out:
        agg.setdefault(lab, []).append((f, a, mn, md, alive))
    cm = np.mean([x[2] for x in agg["control"]])
    cf = np.mean([x[0] for x in agg["control"]])
    print(f"control: {cf:.2f} Hz = {cf*60:.0f} breaths/min, mean output {cm:.1f} Hz\n")
    print(f"{'condition':<24} {'freq Hz':>13} {'%ctrl':>6} {'mean out':>13} "
          f"{'%ctrl':>6} {'mod':>12} {'status':>9}")
    print("-" * 92)
    for lab, _ in CASES:
        v = agg[lab]
        fr = np.array([x[0] for x in v]); mo = np.array([x[2] for x in v])
        md = np.array([x[3] for x in v]); al = np.array([x[4] for x in v], float)
        st = "APNOEA" if al.mean() < 0.5 else ("marginal" if al.mean() < 1.0 else "ok")
        print(f"{lab:<24} {fr.mean():6.2f}±{fr.std():5.2f} {100*fr.mean()/cf:6.0f} "
              f"{mo.mean():6.1f}±{mo.std():5.1f} {100*mo.mean()/cm:6.0f} "
              f"{md.mean():6.2f}±{md.std():4.2f} {st:>9}")
