"""Parallel grid search for a physiological locomotor CPG regime.

Targets (rat hindlimb locomotion):
  period 300-1500 ms, flexor/extensor Mn correlation < -0.5,
  Mn firing 10-50 Hz, duty cycle 0.3-0.7, both half-centres active.
"""
import sys, os, itertools, time; sys.path.insert(0, ".")
import numpy as np
from multiprocessing import Pool
from circuitpharm.cpg import HalfCentreCPG, burst_metrics

T, DT, WARM = 7000.0, 0.1, 2000.0

def run(args):
    drive, gadp, tau_a, w_gly, drive_in = args
    try:
        c = HalfCentreCPG(drive=drive, tau_adapt=tau_a, drive_in=drive_in, seed=1)
        c.WIRING = [(a,b,r,(w_gly if (r=="gly" and b=="RG*") else w))
                    for a,b,r,w in HalfCentreCPG.WIRING]
        for p in c.pops.values(): p.g_adapt = gadp
        for i in range(int(T/DT)):
            c.step(DT)
            if i % 10 == 0: c.record()
        A = c.arrays(); m = A["t"] > WARM
        rf, re_ = A["RG_F"][m], A["RG_E"][m]
        mf, me = A["Mn_F"][m], A["Mn_E"][m]
        if mf.std() < 1e-6 or me.std() < 1e-6: return None
        cm = float(np.corrcoef(mf, me)[0,1])
        bm = burst_metrics(A["t"][m], mf)
        per, duty = bm["period"], bm["duty"]
        if not (per == per) or bm["n_bursts"] < 4: return None
        # score: distance from the physiological target box
        pen = 0.0
        pen += max(0.0, 300-per)/300 + max(0.0, per-1500)/1500
        pen += max(0.0, 10-mf.mean())/10 + max(0.0, mf.mean()-50)/50
        if duty == duty: pen += max(0.0, 0.3-duty)/0.3 + max(0.0, duty-0.7)/0.7
        else: pen += 0.5
        pen += max(0.0, cm + 0.5)              # want cm <= -0.5
        pen += abs(rf.mean()-re_.mean())/max(1.0, rf.mean())   # symmetry
        return dict(score=pen, drive=drive, gadp=gadp, tau_a=tau_a, w_gly=w_gly,
                    drive_in=drive_in, per=per, duty=duty, cm=cm,
                    mnf=float(mf.mean()), mne=float(me.mean()), rgf=float(rf.mean()))
    except Exception:
        return None

GRID = list(itertools.product(
    (200., 280., 360., 450.),      # drive (pA) -- must rise with adaptation
    (1.0, 1.8, 2.8, 4.0),          # g_adapt (nS/spike)
    (400., 700., 1100.),           # tau_adapt (ms)
    (8., 20., 40.),                # w_gly RG reciprocal
    (75., 130.),                   # interneuron drive
))

if __name__ == "__main__":
    n = os.cpu_count()
    print(f"{len(GRID)} combos on {n} cores, {T/1000:.0f}s sim each")
    t0 = time.time()
    with Pool(n) as pool:
        res = [r for r in pool.imap_unordered(run, GRID, chunksize=2) if r]
    el = time.time() - t0
    print(f"{len(res)}/{len(GRID)} produced a rhythm in {el:.0f}s "
          f"({len(GRID)/el:.1f} combos/s, {len(GRID)*T/1000/el:.0f}x realtime aggregate)\n")
    res.sort(key=lambda r: r["score"])
    print(f"{'score':>6} {'drive':>6} {'gadp':>5} {'tauA':>5} {'wgly':>5} {'drvIn':>6} | "
          f"{'per':>6} {'duty':>5} {'Mncor':>6} {'Mnf':>6} {'Mne':>6} {'RGf':>5}")
    for r in res[:12]:
        d = r['duty'] if r['duty']==r['duty'] else float('nan')
        print(f"{r['score']:6.3f} {r['drive']:6.0f} {r['gadp']:5.2f} {r['tau_a']:5.0f} "
              f"{r['w_gly']:5.0f} {r['drive_in']:6.0f} | {r['per']:6.0f} {d:5.2f} "
              f"{r['cm']:+6.2f} {r['mnf']:6.1f} {r['mne']:6.1f} {r['rgf']:5.1f}")
