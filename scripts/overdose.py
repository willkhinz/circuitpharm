"""THE deciding safety test: overdose behaviour of the candidate combination.

At therapeutic dose both GABA-arm options look similar. What separates them is what
happens when someone takes 5-20x the intended amount -- which is the failure mode that
actually kills people with alcohol, and the harm this whole project exists to remove.

A dose multiplier D scales the WHOLE product (both components), as taking more of a
mixture does:
    PAM gain   -> 1 + (gain - 1) * D
    NMDA block -> min(block * D, 1)

Arms:
    BZ-PAM        efficacy ceiling 2.5  (requires endogenous GABA; saturates)
    neurosteroid  efficacy ceiling 6.0  (high-efficacy, partial direct agonism)
    uncapped      no ceiling            (barbiturate-like reference; should die early)

Reported: ventilation vs D, and the apnoea threshold = the overdose index. For reference,
ethanol's own lethal/intoxicating ratio is roughly 5x (20 mM strong intoxication,
~100 mM lethal), so an overdose index above ~5 is an improvement on alcohol and below it
is not.
"""
import sys, os, itertools; sys.path.insert(0, ".")
import numpy as np
from multiprocessing import Pool
from spinal.resp import PreBotC, resp_metrics
from spinal.cpg import Drug

FOREBRAIN_2B, PERIPHERAL_2B, GABA_SENS = 0.70, 0.15, 0.15
RESP_OP = dict(drive=170.0, g_adapt=2.5, tau_adapt=400.0,
               w=dict(ee_ampa=0.45, ee_nmda=0.2475))
N_SEED = 4
DOSES = (1.0, 2.0, 3.0, 5.0, 8.0, 12.0, 18.0, 26.0)
# therapeutic points chosen by scripts/final_combo.py
ARMS = [
    ("BZ-PAM        PAM2.2 + 20%fb", 2.2, 0.20, 2.5),
    ("neurosteroid  PAM1.8 + 20%fb", 1.8, 0.20, 6.0),
    ("uncapped ref  PAM2.2 + 20%fb", 2.2, 0.20, 1e9),
]


def job(a):
    lab, pam0, fb0, cap, D, seed = a
    pam = 1.0 + (pam0 - 1.0) * D
    fb = min(fb0 * D, 1.0)
    b = fb / FOREBRAIN_2B                      # GluN2B-selective arm
    d = Drug(gaba_a_gain=pam, gaba_a_tau=1.0 + 0.6 * (pam - 1.0),
             gaba_a_efficacy_cap=cap, nmda_block=min(b, 1.0),
             glun2b_selectivity=1.0, glun2b_fraction=PERIPHERAL_2B)
    bot = PreBotC(drug=d, gaba_sens=GABA_SENS, seed=seed, **RESP_OP)
    for i in range(int(14000 / 0.1)):
        bot.step(0.1)
        if i % 10 == 0:
            bot.record()
    A = bot.arrays(); m = A["t"] > 4000
    r = resp_metrics(A["t"][m], A["Out"][m], A["Exc"][m])
    return lab, D, r["mean"], r["mod"], float(r["alive"])


if __name__ == "__main__":
    jobs = [(lab, p, f, c, D, s) for (lab, p, f, c) in ARMS
            for D in DOSES for s in range(N_SEED)]
    print(f"{len(jobs)} sims on {os.cpu_count()} cores\n")
    with Pool(os.cpu_count()) as pool:
        out = pool.map(job, jobs)
    agg = {}
    for lab, D, mn, mod, al in out:
        agg.setdefault((lab, D), []).append((mn, mod, al))
    ctrl = {lab: np.mean([x[0] for x in agg[(lab, 1.0)]]) for (lab, _, _, _) in ARMS}

    print(f"{'dose x':>7} | " + " | ".join(f"{lab.split()[0]:>22}" for lab, _, _, _ in ARMS))
    print("-" * 80)
    for D in DOSES:
        row = f"{D:7.1f} |"
        for (lab, _, _, _) in ARMS:
            v = agg[(lab, D)]
            mn = np.mean([x[0] for x in v])
            al = np.mean([float(x[2]) for x in v])
            st = "APNOEA" if al < 0.5 else ("marg" if al < 1.0 else "")
            row += f" {100*mn/ctrl[lab]:16.0f}% {st:>5} |"
        print(row)

    print("\noverdose index (dose multiple at which the rhythm fails):")
    for (lab, _, _, _) in ARMS:
        thr = None
        for D in DOSES:
            if np.mean([float(x[2]) for x in agg[(lab, D)]]) < 0.5:
                thr = D
                break
        ref = "  <-- worse than alcohol (~5x)" if (thr and thr <= 5) else ""
        print(f"  {lab:<32} {thr if thr else f'>{DOSES[-1]:.0f}x (not reached)'}{ref}")
