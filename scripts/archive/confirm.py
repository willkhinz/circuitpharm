"""Final confirmation of the recommended combination, 8 seeds, against reference points.

References included so the candidate can be read against things with known human profiles:
  sedative benzodiazepine alone  -> should sit near -17% ventilation (the calibration anchor)
  GluN2B-selective antagonist alone at 40% forebrain block -> the non-alcohol-like option
  non-selective NMDA version of the candidate -> the cost of getting selectivity wrong
"""
import os
import numpy as np
from multiprocessing import Pool
from circuitpharm.resp import PreBotC, resp_metrics
from circuitpharm.circuit import SpinalCircuit
from circuitpharm.cpg import Drug
import scripts.reflex as R

FB_2B, PERIPH_2B, GABA_SENS, IA_SCALE = 0.70, 0.15, 0.15, 0.30
RESP_OP = dict(drive=170.0, g_adapt=2.5, tau_adapt=400.0,
               w=dict(ee_ampa=0.45, ee_nmda=0.2475))
N_SEED = 8


def drug(pam=1.0, fb=0.0, sel=True, cap=2.5):
    kw = dict(gaba_a_gain=pam, gaba_a_tau=1.0 + 0.6 * (pam - 1.0),
              gaba_a_efficacy_cap=cap)
    if fb > 0:
        kw.update(nmda_block=min(fb / (FB_2B if sel else 1.0), 1.0),
                  glun2b_selectivity=1.0 if sel else 0.0,
                  glun2b_fraction=PERIPH_2B)
    return Drug(**kw)


CASES = [
    ("control",                                        dict()),
    ("ref: sedative benzodiazepine alone (PAM 2.0x)",  dict(pam=2.0)),
    ("CANDIDATE A: BZ-PAM 2.2x + 20%fb GluN2B-sel",    dict(pam=2.2, fb=0.20, sel=True)),
    ("  same but NON-selective NMDA",                  dict(pam=2.2, fb=0.20, sel=False)),
    ("CANDIDATE B: neurosteroid 1.8x + 20%fb 2B-sel",  dict(pam=1.8, fb=0.20, sel=True, cap=6.0)),
    ("OPTION C: GluN2B-sel alone, 40% forebrain",      dict(pam=1.0, fb=0.40, sel=True)),
]


def job(a):
    i, kind, seed = a
    d = drug(**CASES[i][1])
    if kind == "resp":
        b = PreBotC(drug=d, gaba_sens=GABA_SENS, seed=seed, **RESP_OP)
        for k in range(int(14000 / 0.1)):
            b.step(0.1)
            if k % 10 == 0:
                b.record()
        A = b.arrays(); m = A["t"] > 4000
        r = resp_metrics(A["t"][m], A["Out"][m], A["Exc"][m])
        return i, kind, r["mean"], r["amp"]
    orig = dict(SpinalCircuit.W)
    try:
        SpinalCircuit.W = dict(orig)
        for rec in ("ampa", "nmda"):
            SpinalCircuit.W[("Ia", "Mn", rec)] = orig[("Ia", "Mn", rec)] * IA_SCALE
        r = R.run(d, seed=seed)
        return i, kind, (r["mn_dyn"] - r["mn_base"]) / max(1e-6, r["ia_dyn"] - r["ia_base"]), 0.0
    finally:
        SpinalCircuit.W = orig


if __name__ == "__main__":
    jobs = [(i, k, s) for i in range(len(CASES))
            for k in ("resp", "reflex") for s in range(N_SEED)]
    with Pool(os.cpu_count()) as p:
        out = p.map(job, jobs)
    A = {}
    for i, k, v, amp in out:
        A.setdefault((i, k), []).append((v, amp))
    cv = np.mean([x[0] for x in A[(0, "resp")]])
    ca = np.mean([x[1] for x in A[(0, "resp")]])
    cr = np.mean([x[0] for x in A[(0, "reflex")]])
    print(f"controls (n={N_SEED}): ventilation {cv:.2f}, amplitude {ca:.1f}, "
          f"reflex gain {cr:.3f}\n")
    print(f"{'condition':<48} {'ventilation':>16} {'amplitude':>14} {'reflex':>16}")
    print("-" * 96)
    for i, (lab, _) in enumerate(CASES):
        rv = np.array([x[0] for x in A[(i, "resp")]])
        ra = np.array([x[1] for x in A[(i, "resp")]])
        rf = np.array([x[0] for x in A[(i, "reflex")]])
        pv = 100 * rv / cv; pa = 100 * ra / ca; pf = 100 * rf / cr
        print(f"{lab:<48} {pv.mean():8.0f}±{pv.std():3.0f}% "
              f"{pa.mean():8.0f}±{pa.std():3.0f}% {pf.mean():10.0f}±{pf.std():3.0f}%")
