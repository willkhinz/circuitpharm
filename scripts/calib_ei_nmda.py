"""Does an NMDA component on Exc->Inh let NMDA block spare ventilation?

WHY THIS SWEEP EXISTS. The earlier conclusion -- "this preBotC cannot represent the
clinical fact that NMDA antagonists spare ventilation, and no parameter fixes it" -- was
reached by sweeping the NMDA share of RECURRENT EXCITATION (ee_nmda) only. With no NMDA on
the excitatory->inhibitory projection, blocking NMDA could not reduce inhibition by any
amount, so the DISINHIBITION arm of ketamine's action was structurally absent from the
circuit. That makes the earlier conclusion unsupported rather than wrong: the mechanism was
never in the model to begin with.

preBotC inhibitory interneurons do express NMDA receptors, so a nonzero ei_nmda is the
biologically correct structure. This sweep asks whether it changes the SIGN of the
ventilation response to NMDA block.

Two opposing effects are now both represented:
  - block reduces recurrent excitation  -> weaker/shorter bursts  (depressant)
  - block reduces drive to Inh          -> less inhibition on Exc (stimulant)

TOTAL exc->inh weight is held constant across shares so the DRUG-FREE rhythm is unchanged
and every column is comparable to its own control. Non-selective block (glun2b_sel=0) is
used because it is the harshest test and the one with the clearest clinical anchor.

Run:  python scripts/calib_ei_nmda.py
"""
import sys, os, itertools; sys.path.insert(0, ".")
import numpy as np
from multiprocessing import Pool
from spinal.resp import PreBotC, resp_metrics
from spinal.cpg import Drug

# operating point used everywhere else in the project (simulator.RESP_OP)
EE_AMPA, EE_NMDA = 0.45, 0.2475
EI_TOTAL = 0.55                    # held constant; only its AMPA/NMDA split varies
SHARES = (0.0, 0.2, 0.5, 1.0, 2.0)  # ei_nmda / ei_ampa
BLOCKS = (0.0, 0.3, 0.6, 0.9)
N_SEED = 3
DUR_MS, DT = 14000.0, 0.1


def job(a):
    share, blk, seed = a
    ei_ampa = EI_TOTAL / (1.0 + share)
    ei_nmda = ei_ampa * share
    d = Drug(nmda_block=blk, glun2b_selectivity=0.0)   # non-selective = harshest
    b = PreBotC(drug=d, gaba_sens=0.15, seed=seed,
                drive=170.0, g_adapt=2.5, tau_adapt=400.0,
                w=dict(ee_ampa=EE_AMPA, ee_nmda=EE_NMDA,
                       ei_ampa=ei_ampa, ei_nmda=ei_nmda))
    for i in range(int(DUR_MS / DT)):
        b.step(DT)
        if i % 10 == 0:
            b.record()
    A = b.arrays(); m = A["t"] > 4000
    r = resp_metrics(A["t"][m], A["Out"][m], A["Exc"][m])
    return share, blk, r["mean"], r["mod"], r["freq"], float(r["alive"])


if __name__ == "__main__":
    jobs = [(s, b, k) for s, b in itertools.product(SHARES, BLOCKS)
            for k in range(N_SEED)]
    print(f"{len(jobs)} sims on {os.cpu_count()} cores "
          f"({len(SHARES)} ei_nmda shares x {len(BLOCKS)} blocks x {N_SEED} seeds)")
    with Pool(os.cpu_count()) as p:
        out = p.map(job, jobs)
    A = {}
    for s, b, mn, md, fq, al in out:
        A.setdefault((s, b), []).append((mn, md, fq, al))

    print("\nCLINICAL ANCHOR: non-selective NMDA block should leave ventilation "
          "NEAR OR ABOVE control.\n")
    print(f"{'ei_nmda share':>13} {'ctrl mod':>9} {'ctrl Hz':>8} |"
          + "".join(f"{f'blk {b:.0%}':>11}" for b in BLOCKS[1:]))
    print("-" * 70)
    for s in SHARES:
        ctrl = A[(s, 0.0)]
        c = np.mean([x[0] for x in ctrl])
        cm = np.mean([x[1] for x in ctrl])
        cf = np.mean([x[2] for x in ctrl])
        row = f"{s:13.2f} {cm:9.2f} {cf:8.2f} |"
        for b in BLOCKS[1:]:
            v = np.mean([x[0] for x in A[(s, b)]])
            alive = np.mean([x[3] for x in A[(s, b)]])
            cell = f"{100*v/c:.0f}%" + ("" if alive > 0.5 else "*")
            row += f"{cell:>11}"
        print(row)
    print("\npercentages are ventilation vs the drug-free control AT THAT SHARE")
    print("'ctrl mod' and 'ctrl Hz' confirm the drug-free rhythm is intact at each share")
    print("'*' = rhythm declared dead (modulation or mean collapsed)")
    print("\nshare 0.00 is the model AS IT STOOD when the limitation was declared.")
