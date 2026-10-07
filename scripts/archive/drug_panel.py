"""First scientific output: dose-response of the spinal circuit to each mechanism,
and to the GABA-A PAM + NMDA antagonist combination.

Readouts are the 'data from the simulated brain': motoneuron pool rates, locomotor
period, duty cycle, flexor/extensor alternation, and whether the rhythm survives at all.
"""
import os
import numpy as np
from multiprocessing import Pool
from circuitpharm.circuit import SpinalCircuit
from circuitpharm.cpg import Drug, burst_metrics

GAIN, T, DT = 900.0, 9000.0, 0.1
SMOOTH = np.ones(60)/60.0

def metrics(label, drug):
    c = SpinalCircuit(drug=drug, rg_gain=GAIN, seed=1)
    for i in range(int(T/DT)):
        c.step(DT)
        if i % 10 == 0: c.record()
    A = c.arrays(); m = A["t"] > 3000
    sm = lambda v: np.convolve(v, SMOOTH, mode="same")
    mf, me = sm(A["Mn_F"][m]), sm(A["Mn_E"][m])
    rg = A["RG_F"][m]
    alive = mf.std() > 1e-6 and me.std() > 1e-6
    bm = burst_metrics(A["t"][m], mf) if alive else dict(period=np.nan, duty=np.nan, n_bursts=0)
    return dict(label=label,
                mn=float((mf.mean()+me.mean())/2), mnpeak=float(mf.max()),
                per=bm["period"], duty=bm["duty"],
                corr=float(np.corrcoef(mf, me)[0,1]) if alive else np.nan,
                rg_amp=float(rg.max()), rhythm=bool(bm["n_bursts"] >= 3))

CASES = [
    ("control",                         Drug()),
    # --- GABA-A PAM alone (gain + decay prolongation, with efficacy ceiling) ---
    ("GABA-A PAM  1.5x",                Drug(gaba_a_gain=1.5, gaba_a_tau=1.3)),
    ("GABA-A PAM  2.0x",                Drug(gaba_a_gain=2.0, gaba_a_tau=1.6)),
    ("GABA-A PAM  3.0x (past ceiling)", Drug(gaba_a_gain=3.0, gaba_a_tau=2.0)),
    # --- NMDA antagonist alone, non-selective (MK-801 / ketamine-like) ---
    ("NMDA block 30% non-selective",    Drug(nmda_block=0.30, glun2b_selectivity=0.0)),
    ("NMDA block 60% non-selective",    Drug(nmda_block=0.60, glun2b_selectivity=0.0)),
    # --- NMDA antagonist, GluN2B-selective; spinal cord is GluN2D-dominant ---
    ("NMDA 60% GluN2B-sel (spinal, 2B=0.15)", Drug(nmda_block=0.60, glun2b_selectivity=1.0,
                                                   glun2b_fraction=0.15)),
    ("NMDA 60% GluN2B-sel (forebrain-like 2B=0.7)", Drug(nmda_block=0.60, glun2b_selectivity=1.0,
                                                          glun2b_fraction=0.7)),
    # --- ethanol-like: GABA-A + NMDA + glycine potentiation together ---
    ("ethanol-like (GABA+NMDA+GlyR)",   Drug(gaba_a_gain=1.6, gaba_a_tau=1.4,
                                             nmda_block=0.15, glun2b_selectivity=0.0,
                                             glyr_gain=1.3)),
    # --- the candidate: GABA-A PAM + GluN2B-selective antagonist ---
    ("candidate: PAM 2x + 2B-sel 60%",  Drug(gaba_a_gain=2.0, gaba_a_tau=1.6,
                                             nmda_block=0.60, glun2b_selectivity=1.0,
                                             glun2b_fraction=0.15)),
    ("candidate: PAM 2x + non-sel 60%", Drug(gaba_a_gain=2.0, gaba_a_tau=1.6,
                                             nmda_block=0.60, glun2b_selectivity=0.0)),
]

if __name__ == "__main__":
    with Pool(os.cpu_count()) as p:
        res = p.starmap(metrics, CASES)
    ctrl = res[0]
    print(f"{'condition':<44} {'Mn Hz':>6} {'%ctrl':>6} {'peak':>6} {'per ms':>7} "
          f"{'duty':>5} {'corr':>6} {'rhythm':>7}")
    print("-"*96)
    for r in res:
        f = lambda v: v if v == v else float('nan')
        pct = 100*r['mn']/ctrl['mn'] if ctrl['mn'] > 0 else float('nan')
        print(f"{r['label']:<44} {r['mn']:6.1f} {pct:6.0f} {r['mnpeak']:6.0f} "
              f"{f(r['per']):7.0f} {f(r['duty']):5.2f} {f(r['corr']):+6.2f} "
              f"{'yes' if r['rhythm'] else 'LOST':>7}")
