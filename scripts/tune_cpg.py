import sys; sys.path.insert(0,".")
import numpy as np
from spinal.cpg import HalfCentreCPG, burst_metrics

def probe(drive, w_gly, gadp, tau_a, T=10000., dt=0.1):
    c = HalfCentreCPG(drive=drive, tau_adapt=tau_a, seed=1)
    c.WIRING=[(a,b,r,(w_gly if (r=="gly" and b=="RG*") else w)) for a,b,r,w in HalfCentreCPG.WIRING]
    for p in c.pops.values(): p.g_adapt = gadp
    for i in range(int(T/dt)):
        c.step(dt)
        if i%10==0: c.record()
    A=c.arrays(); m=A["t"]>2500
    rf,re_,mf,me=A["RG_F"][m],A["RG_E"][m],A["Mn_F"][m],A["Mn_E"][m]
    cr = np.corrcoef(rf,re_)[0,1] if rf.std()>1e-6 and re_.std()>1e-6 else np.nan
    cm = np.corrcoef(mf,me)[0,1] if mf.std()>1e-6 and me.std()>1e-6 else np.nan
    bm = burst_metrics(A["t"][m], rf)
    return rf.mean(), A["InRG_F"][m].mean(), mf.mean(), cr, cm, bm["period"], bm["n_bursts"], bm["duty"]

print(f"{'drv':>4} {'wgly':>4} {'gadp':>5} {'tauA':>5} | {'RGf':>5} {'InRG':>5} {'Mnf':>5} "
      f"{'RGcor':>6} {'Mncor':>6} {'per':>6} {'nb':>3} {'duty':>5}")
best=None
for drive in (195.,):
    for w_gly in (4., 12., 28.):
        for gadp in (0.45,):
            for tau_a in (500., 900., 1600.):
                r=probe(drive,w_gly,gadp,tau_a)
                per=r[5] if r[5]==r[5] else float('nan'); dty=r[7] if r[7]==r[7] else float('nan')
                print(f"{drive:4.0f} {w_gly:4.0f} {gadp:5.2f} {tau_a:5.0f} | {r[0]:5.1f} {r[1]:5.1f} "
                      f"{r[2]:5.1f} {r[3]:+6.2f} {r[4]:+6.2f} {per:6.0f} {r[6]:3d} {dty:5.2f}")
                score = (r[4] if r[4]==r[4] else 0)
                if r[1]>2 and r[2]>2 and score<-0.3 and 250<per<2000:
                    if best is None or score<best[0]: best=(score,drive,w_gly,gadp,tau_a,per)
print("\nbest (Mncorr, drive, wgly, gadp, tauA, period):", best)
