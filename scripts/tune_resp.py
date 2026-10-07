"""Tune the preBotC to a bursting eupnoeic regime.

Key constraint learned the hard way (twice -- same error as the locomotor RG): the tonic
drive MUST be sub-rheobase (~200 pA with tonic GABA). A group pacemaker needs a silent
interburst phase so that recurrent excitation can recruit the population regeneratively.
With suprathreshold drive the population fires tonically and only ripples -- which the
old relative-threshold metric then mis-reported as a fast rhythm.
"""
import sys, os, itertools; sys.path.insert(0, ".")
import numpy as np
from multiprocessing import Pool
from spinal.resp import PreBotC, resp_metrics


def probe(a):
    drive, gadp, tau_a, wee = a
    b = PreBotC(drive=drive, g_adapt=gadp, tau_adapt=tau_a,
                w=dict(ee_ampa=wee, ee_nmda=wee * 0.55), seed=1)
    dt = 0.1
    for i in range(int(14000 / dt)):
        b.step(dt)
        if i % 10 == 0:
            b.record()
    A = b.arrays(); m = A["t"] > 4000
    r = resp_metrics(A["t"][m], A["Out"][m], A["Exc"][m])
    return dict(drive=drive, gadp=gadp, tau_a=tau_a, wee=wee,
                exc=float(A["Exc"][m].mean()), **r)


GRID = list(itertools.product(
    (90., 130., 170.),            # sub-rheobase drive (pA)
    (1.2, 2.5),                   # g_adapt
    (400., 700.),                 # tau_adapt (ms)
    (0.25, 0.45, 0.70),           # recurrent excitation weight
))

if __name__ == "__main__":
    with Pool(os.cpu_count()) as p:
        res = p.map(probe, GRID)
    ok = [r for r in res if r["alive"] and 0.6 < r["freq"] < 3.0 and r["mod"] > 1.2]
    ok.sort(key=lambda r: -r["mod"])
    print(f"{len(ok)}/{len(GRID)} bursting in eupnoea range (0.6-3 Hz) with mod>1.2")
    hdr = (f"{'drive':>6} {'gadp':>5} {'tauA':>5} {'wee':>5} | {'freq':>6} {'mod':>5} "
           f"{'amp':>6} {'mean':>6} {'Exc':>6} {'n':>3}")
    print(hdr)
    for r in (ok[:10] if ok else sorted(res, key=lambda r: -r["mod"])[:10]):
        print(f"{r['drive']:6.0f} {r['gadp']:5.1f} {r['tau_a']:5.0f} {r['wee']:5.2f} | "
              f"{r['freq']:6.2f} {r['mod']:5.2f} {r['amp']:6.1f} {r['mean']:6.1f} "
              f"{r['exc']:6.1f} {r['n']:3d}")
    if not ok:
        print("\n(no hits -- rows above are best-by-modulation)")
