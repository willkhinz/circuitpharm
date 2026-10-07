"""Calibrate the NMDA share of preBotC recurrent excitation against clinical fact.

ANCHOR: ketamine-class NMDA channel blockers do NOT depress ventilation. Clinically they
preserve respiratory drive and airway reflexes, and s-ketamine STIMULATES breathing and
attenuates propofol-/opioid-induced hypoventilation. Mechanistically this is expected:
preBotC rhythmogenesis is AMPA/kainate-dependent -- AMPA antagonists abolish the rhythm,
NMDA antagonists do not.

The model was built with ee_nmda/ee_ampa = 0.55, i.e. NMDA carrying 55% of recurrent
excitation, which forces the model to predict large respiratory depression from NMDA block.
That is a STRUCTURAL error, not a parameter to fudge with a 'respiratory boost' term.

This sweep finds the NMDA share at which a full non-selective NMDA block leaves ventilation
near control, as clinically observed. Total excitatory drive is held constant so the
drug-free rhythm is unchanged.
"""
import os, itertools
import numpy as np
from multiprocessing import Pool
from circuitpharm.resp import PreBotC, resp_metrics
from circuitpharm.cpg import Drug

BASE_TOTAL = 0.45 + 0.2475     # keep total excitatory weight constant
SHARES = (0.55, 0.35, 0.20, 0.12, 0.06)
BLOCKS = (0.0, 0.3, 0.6, 0.9)
N_SEED = 3

def job(a):
    share, blk, seed = a
    ampa = BASE_TOTAL / (1.0 + share)
    nmda = ampa * share
    d = Drug(nmda_block=blk, glun2b_selectivity=0.0)   # non-selective, full peripheral hit
    b = PreBotC(drug=d, gaba_sens=0.15, seed=seed,
                drive=170.0, g_adapt=2.5, tau_adapt=400.0,
                w=dict(ee_ampa=ampa, ee_nmda=nmda))
    for i in range(int(14000/0.1)):
        b.step(0.1)
        if i%10==0: b.record()
    A=b.arrays(); m=A["t"]>4000
    r=resp_metrics(A["t"][m],A["Out"][m],A["Exc"][m])
    return share, blk, r["mean"], r["mod"]

if __name__=="__main__":
    with Pool(os.cpu_count()) as p:
        out=p.map(job,[(s,b,k) for s,b in itertools.product(SHARES,BLOCKS) for k in range(N_SEED)])
    A={}
    for s,b,mn,md in out: A.setdefault((s,b),[]).append((mn,md))
    print("target: full non-selective NMDA block leaves ventilation NEAR CONTROL\n")
    print(f"{'NMDA share':>11} {'ctrl mod':>9} |" + "".join(f"{f'blk {b:.0%}':>11}" for b in BLOCKS[1:]))
    print("-"*58)
    for s in SHARES:
        c=np.mean([x[0] for x in A[(s,0.0)]]); cm=np.mean([x[1] for x in A[(s,0.0)]])
        row=f"{s:11.2f} {cm:9.2f} |"
        for b in BLOCKS[1:]:
            v=np.mean([x[0] for x in A[(s,b)]])
            row+=f"{100*v/c:10.0f}%"
        print(row)
    print("\n(percentages are ventilation vs the drug-free control AT THAT SHARE;")
    print(" 'ctrl mod' confirms the drug-free rhythm is still bursting at each share)")
