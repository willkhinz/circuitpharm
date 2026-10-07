"""THE key respiratory safety validation: the benzodiazepine / barbiturate ceiling
dissociation.

A positive allosteric modulator requires endogenous GABA and therefore saturates, so its
respiratory depression plateaus. A direct agonist does not require GABA and has no
ceiling, so it proceeds to apnoea. This is the textbook reason benzodiazepines are safer
in overdose than barbiturates. If the model reproduces *why* one plateaus and the other
does not, it can be trusted on a novel compound's ceiling -- which is the single property
the product's overdose safety depends on.

Implementation: a PAM is `gaba_a_efficacy_cap = 2.5`; a direct agonist is the same
mechanism with the cap removed (set very high). Nothing else differs, so this is a clean
paired comparison.

Primary readout is MEAN OUTPUT (minute-ventilation proxy = rate x amplitude). Frequency
is NOT reported: the FFT peak follows fast fluctuations once the burst pattern degrades,
so it rises during depression and is unusable in this regime.
"""
import os
import numpy as np
from multiprocessing import Pool
from circuitpharm.resp import PreBotC, resp_metrics
from circuitpharm.cpg import Drug

OP = dict(drive=170.0, g_adapt=2.5, tau_adapt=400.0,
          w=dict(ee_ampa=0.45, ee_nmda=0.2475))
N_SEED, T, DT = 4, 14000.0, 0.1
DOSES = (1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 6.0, 9.0, 14.0)


def one(args):
    kind, x, seed = args
    cap = 2.5 if kind == "PAM" else 1e9          # direct agonist: no ceiling
    d = Drug(gaba_a_gain=x, gaba_a_tau=1.0 + 0.6 * (x - 1.0),
             gaba_a_efficacy_cap=cap)
    b = PreBotC(drug=d, seed=seed, **OP)
    for i in range(int(T / DT)):
        b.step(DT)
        if i % 10 == 0:
            b.record()
    A = b.arrays(); m = A["t"] > 4000
    r = resp_metrics(A["t"][m], A["Out"][m], A["Exc"][m])
    return (kind, x, r["mean"], r["mod"], r["alive"])


if __name__ == "__main__":
    jobs = [(k, x, s) for k in ("PAM", "agonist") for x in DOSES
            for s in range(N_SEED)]
    with Pool(os.cpu_count()) as p:
        out = p.map(one, jobs)
    agg = {}
    for k, x, mn, md, al in out:
        agg.setdefault((k, x), []).append((mn, md, al))
    ctrl = np.mean([v[0] for v in agg[("PAM", 1.0)]])
    print(f"control mean output = {ctrl:.1f} Hz (minute-ventilation proxy)\n")
    print(f"{'dose':>6} | {'PAM (cap 2.5)':>28} | {'direct agonist (no cap)':>28}")
    print(f"{'':>6} | {'vent':>8} {'%ctrl':>6} {'mod':>5} {'st':>5} | "
          f"{'vent':>8} {'%ctrl':>6} {'mod':>5} {'st':>5}")
    print("-" * 76)
    for x in DOSES:
        row = f"{x:6.1f} |"
        for k in ("PAM", "agonist"):
            v = agg[(k, x)]
            mn = np.mean([a[0] for a in v]); md = np.mean([a[1] for a in v])
            al = np.mean([float(a[2]) for a in v])
            st = "APN" if al < 0.5 else ("marg" if al < 1.0 else "ok")
            row += f" {mn:8.1f} {100*mn/ctrl:6.0f} {md:5.2f} {st:>5} |"
        print(row)

    # apnoea threshold per arm
    print()
    for k in ("PAM", "agonist"):
        thr = None
        for x in DOSES:
            al = np.mean([float(a[2]) for a in agg[(k, x)]])
            if al < 0.5:
                thr = x; break
        print(f"{k:>16}: apnoea threshold = "
              f"{thr if thr else 'NOT REACHED up to %.0fx' % DOSES[-1]}")
